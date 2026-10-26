"""Read access to a built index (a directory of Parquet tables named as in schema.TABLES, or in-memory tables)."""
import os

import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.dataset as pds
import pyarrow.parquet as pq

from . import schema


class Index:
    def __init__(self, tables=None, root=None):
        self.root = root
        self._t = dict(tables or {})

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

    def evidence_for_variants(self, keys):
        """{variant_key: [evidence row dicts]} for mapped variant-level rows."""
        ev = self.table('evidence')
        keys = list(set(keys))
        if not keys or ev.num_rows == 0:
            return {}
        mask = pc.and_(pc.equal(ev['level'], 'variant'), pc.is_in(ev['subject_id'], value_set=pa.array(keys, pa.string())))
        out = {}
        for r in ev.filter(mask).to_pylist():
            out.setdefault(r['subject_id'], []).append(r)
        return out

    def evidence_for(self, level, ids):
        ev = self.table('evidence')
        ids = list(set(ids))
        if not ids or ev.num_rows == 0:
            return {}
        mask = pc.and_(pc.equal(ev['level'], level), pc.is_in(ev['subject_id'], value_set=pa.array(ids, pa.string())))
        out = {}
        for r in ev.filter(mask).to_pylist():
            out.setdefault(r['subject_id'], []).append(r)
        return out

    def events_overlapping(self, chrom, start, end):
        t = self.table('entity_event')
        if t.num_rows == 0:
            return []
        m = pc.and_(pc.equal(t['chrom'], chrom), pc.and_(pc.less_equal(t['start'], end), pc.greater_equal(t['end'], start)))
        return t.filter(m).to_pylist()

    def regions_overlapping(self, chrom, start, end, kinds=None):
        t = self.table('entity_region')
        if t.num_rows == 0:
            return []
        m = pc.and_(pc.equal(t['chrom'], chrom), pc.and_(pc.less_equal(t['start'], end), pc.greater_equal(t['end'], start)))
        rows = t.filter(m).to_pylist()
        return [r for r in rows if kinds is None or r['kind'] in kinds]

    def genes_overlapping(self, chrom, start, end):
        t = self.table('entity_gene')
        if t.num_rows == 0:
            return []
        m = pc.and_(pc.equal(t['chrom'], chrom), pc.and_(pc.less_equal(t['start'], end), pc.greater_equal(t['end'], start)))
        return t.filter(m).to_pylist()

    def haplotypes(self, gene_symbol=None, system=None):
        t = self.table('entity_haplotype')
        rows = t.to_pylist()
        return [r for r in rows if (gene_symbol is None or r['gene_symbol'] == gene_symbol) and (system is None or r['system'] == system)]

    def haplotype_definitions(self, haplotype_ids):
        t = self.table('haplotype_definition')
        ids = set(haplotype_ids)
        out = {}
        for r in t.to_pylist():
            if r['haplotype_id'] in ids:
                out.setdefault(r['haplotype_id'], []).append(r)
        return out

    def variant_genes(self, keys):
        t = self.table('variant_gene')
        ks = set(keys)
        out = {}
        for r in t.to_pylist():
            if r['variant_key'] in ks:
                out.setdefault(r['variant_key'], []).append(r)
        return out

    def mnv_components(self):
        """{component SNV key: [index MNV keys]} for MNV entities in the index (used for phase-aware component matching)."""
        t = self.table('entity_variant')
        out = {}
        if t.num_rows == 0:
            return out
        for r in t.filter(pc.equal(t['vtype'], 'mnv')).to_pylist():
            for i, (a, b) in enumerate(zip(r['ref'], r['alt'])):
                if a != b:
                    out.setdefault(f"{r['chrom']}:{r['pos'] + i}:{a}:{b}", []).append(r['variant_key'])
        return out
