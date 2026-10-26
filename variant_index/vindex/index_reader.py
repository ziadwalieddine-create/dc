"""Read access to a built index (a directory of Parquet tables named as in schema.TABLES, or in-memory tables).

An index on disk is never loaded whole: evidence runs to tens of millions of rows and the difficult-region stratifications to millions of
intervals. Reads go through Arrow datasets with filters (a chromosome, a set of subject ids, a kind), and overlap queries use a per-chromosome
interval structure that is built on first use and dropped when the caller moves on to other chromosomes (calls arrive in VCF order). Tables passed
in memory (tests, small indexes) are filtered the same way, so both give the same answers.
"""
import os
from collections import OrderedDict

import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.dataset as pds
import pyarrow.parquet as pq

from . import schema
from .intervals import IntervalStore

KEEP_CHROMS = 3         # per-chromosome structures kept per table


class _Overlaps:
    """Interval structure over the rows of one chromosome of one table; `kind_field` splits it into one store per kind (regions)."""

    def __init__(self, table, kind_field=None):
        self.table = table
        starts, ends = table['start'].to_pylist(), table['end'].to_pylist()
        kinds = table[kind_field].to_pylist() if kind_field else [None] * table.num_rows
        groups = {}
        for i, (s, e, k) in enumerate(zip(starts, ends, kinds)):
            if s is None or e is None:
                continue
            groups.setdefault(k, []).append(i)
        self.stores = {k: (IntervalStore([starts[i] for i in idx], [max(ends[i], starts[i]) for i in idx]), idx) for k, idx in groups.items()}

    def rows(self, start, end, kinds=None):
        rows = []
        for k, (store, idx) in self.stores.items():
            if kinds is not None and k not in kinds:
                continue
            rows.extend(idx[j] for j in store.overlapping(start, end))
        rows.sort()
        return [self.table.slice(i, 1).to_pylist()[0] for i in rows]


class Index:
    def __init__(self, tables=None, root=None):
        self.root = root
        self._t = dict(tables or {})
        self._c = OrderedDict()         # (table, chrom) -> _Overlaps
        self._defs = None
        self._mnv = None

    @classmethod
    def open(cls, root):
        return cls(root=root)

    @classmethod
    def from_tables(cls, **tables):
        return cls(tables=tables)

    def table(self, name):
        if name not in self._t:
            if self.root is None:
                self._t[name] = schema.empty_table(name)
            else:
                p = os.path.join(self.root, name)
                if os.path.isdir(p):
                    self._t[name] = pds.dataset(p, format='parquet', ignore_prefixes=[]).to_table()
                elif os.path.exists(p + '.parquet'):
                    self._t[name] = pq.read_table(p + '.parquet')
                else:
                    self._t[name] = schema.empty_table(name)
        return self._t[name]

    def _on_disk(self, name, sub=''):
        """Dataset over (a part of) a table directory when the index is on disk and the table has not been loaded; else None."""
        if self.root is None or name in self._t:
            return None
        p = os.path.join(self.root, name, sub)
        return pds.dataset(p, format='parquet', ignore_prefixes=[]) if os.path.isdir(p) else None

    def filtered(self, name, expr, sub='', columns=None):
        """Rows of table `name` for which the Arrow expression holds, without loading the whole table when the index is on disk. `columns` limits what is read."""
        if self.root is not None and name not in self._t:
            ds = self._on_disk(name, sub)
            if ds is not None:
                return ds.to_table(filter=expr, columns=columns)
            if os.path.isdir(os.path.join(self.root, name)):          # a partition that does not exist holds no rows
                return schema.empty_table(name).select(columns) if columns else schema.empty_table(name)
        t = self.table(name)
        if not t.num_rows:
            return t.select(columns) if columns else t
        return pds.dataset(t).to_table(filter=expr, columns=columns)

    # ---- evidence ---------------------------------------------------------------------------------------------------------------------
    def evidence_for_variants(self, keys):
        """{variant_key: [evidence row dicts]} for mapped variant-level rows. Reads one chromosome partition at a time."""
        by_part = {}
        for k in set(keys):
            by_part.setdefault(schema.partition_of(k.split(':', 1)[0]), []).append(k)
        out = {}
        for part, ks in sorted(by_part.items()):
            expr = pc.field('subject_id').isin(ks)
            if self.root is not None and 'evidence' not in self._t:
                t = self.filtered('evidence', expr, sub=os.path.join('variant', part))
            else:
                t = self.filtered('evidence', expr & (pc.field('level') == 'variant'))
            for r in t.to_pylist():
                out.setdefault(r['subject_id'], []).append(r)
        for rows in out.values():
            rows.sort(key=lambda r: r['evidence_id'])
        return out

    def evidence_for(self, level, ids):
        ids = list(set(ids))
        if not ids:
            return {}
        if self.root is not None and 'evidence' not in self._t:
            t = self.filtered('evidence', pc.field('subject_id').isin(ids), sub=level)
        else:
            t = self.filtered('evidence', pc.field('subject_id').isin(ids) & (pc.field('level') == level))
        out = {}
        for r in t.to_pylist():
            out.setdefault(r['subject_id'], []).append(r)
        for rows in out.values():
            rows.sort(key=lambda r: r['evidence_id'])
        return out

    # ---- overlap queries ----------------------------------------------------------------------------------------------------------------
    def _overlaps(self, name, chrom, kind_field=None):
        key = (name, chrom)
        if key in self._c:
            self._c.move_to_end(key)
            return self._c[key]
        o = _Overlaps(self.filtered(name, pc.field('chrom') == chrom), kind_field)
        self._c[key] = o
        same = [k for k in self._c if k[0] == name]
        for k in same[:-KEEP_CHROMS]:
            del self._c[k]
        return o

    def events_overlapping(self, chrom, start, end):
        return self._overlaps('entity_event', chrom).rows(start, end)

    def regions_overlapping(self, chrom, start, end, kinds=None):
        return self._overlaps('entity_region', chrom, 'kind').rows(start, end, None if kinds is None else set(kinds))

    def genes_overlapping(self, chrom, start, end):
        return self._overlaps('entity_gene', chrom).rows(start, end)

    def regions_by_kind(self, kinds):
        """Region rows of the given kinds without loading the whole region table (the difficult-region stratifications alone are millions of intervals)."""
        return self.filtered('entity_region', pc.field('kind').isin(list(kinds))).to_pylist()

    # ---- haplotypes, genes, variants --------------------------------------------------------------------------------------------------
    def haplotypes(self, gene_symbol=None, system=None, levels=None):
        t = self.table('entity_haplotype')
        rows = t.to_pylist()
        return [r for r in rows if (gene_symbol is None or r['gene_symbol'] == gene_symbol) and (system is None or r['system'] == system)
                and (levels is None or r.get('allele_level') in levels)]

    def haplotype_definitions(self, haplotype_ids):
        ids = set(haplotype_ids)
        if self._defs is None:
            self._defs = {}
            for r in self.table('haplotype_definition').to_pylist():
                self._defs.setdefault(r['haplotype_id'], []).append(r)
        return {h: self._defs[h] for h in ids if h in self._defs}

    def variant_partitions(self):
        """Partitions that hold mapped variant-level evidence: contigs of the primary assembly and/or 'NA' (see schema.partition_of)."""
        if self.root is not None and 'evidence' not in self._t:
            d = os.path.join(self.root, 'evidence', 'variant')
            return sorted(os.listdir(d)) if os.path.isdir(d) else []
        t = self.table('evidence')
        if not t.num_rows:
            return []
        t = t.filter((pc.field('level') == 'variant') & (pc.field('mapping_status') == 'mapped'))
        return sorted({schema.partition_of(c) for c in pc.unique(t['chrom']).to_pylist() if c})

    def variant_genes(self, keys):
        ks = list(set(keys))
        out = {}
        if ks:
            for r in self.filtered('variant_gene', pc.field('variant_key').isin(ks)).to_pylist():
                out.setdefault(r['variant_key'], []).append(r)
        return out

    def mnv_components(self):
        """{component SNV key: [index MNV keys]} for MNV entities in the index (used for phase-aware component matching)."""
        if self._mnv is None:
            out = {}
            for r in self.filtered('entity_variant', pc.field('vtype') == 'mnv').to_pylist():
                for i, (a, b) in enumerate(zip(r['ref'], r['alt'])):
                    if a != b:
                        out.setdefault(f"{r['chrom']}:{r['pos'] + i}:{a}:{b}", []).append(r['variant_key'])
            self._mnv = out
        return self._mnv
