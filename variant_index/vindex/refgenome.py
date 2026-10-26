"""Reference sequence access, contig naming and assembly validation.

The canonical contig name is the one used in the mirrored Ensembl FASTA ('1'..'22', 'X', 'Y', 'MT', and accession-style
names for scaffolds). VCFs name contigs in other styles (UCSC 'chr1', 'chrM', 'chr1_KI270706v1_random'; GenBank and RefSeq
accessions). Primary chromosomes are mapped by rule; every other name is mapped through the NCBI assembly report, never
by guessing. A contig that cannot be mapped is reported, not dropped.
"""
import csv
import os
import re

import pysam

PRIMARY = [str(i) for i in range(1, 23)] + ['X', 'Y', 'MT']
CHR1_LENGTH_BUILD = {248956422: 'GRCh38', 249250621: 'GRCh37', 247249719: 'NCBI36', 248387328: 'T2T-CHM13v2.0'}


def read_assembly_report(path):
    """Rows of an NCBI assembly report as dicts (the header is the last '#' line)."""
    header, rows = None, []
    with open(path, encoding='utf-8', newline='') as f:
        for line in f:
            line = line.rstrip('\r\n')
            if line.startswith('#'):
                if line.startswith('# Sequence-Name'):
                    header = line[2:].split('\t')
                continue
            if line and header:
                rows.append(dict(zip(header, line.split('\t'))))
    return rows


class Reference:
    """A bgzip/plain FASTA with .fai, plus contig-name resolution."""

    def __init__(self, fasta_path, assembly_report=None, name='GRCh38', builds=None):
        self.path = fasta_path
        self.name = name
        self.builds = dict(CHR1_LENGTH_BUILD if builds is None else builds)
        self.fa = pysam.FastaFile(fasta_path)
        self.lengths = dict(zip(self.fa.references, self.fa.lengths))
        self._alias = {}
        for c in self.lengths:
            self._alias[c] = c
        if assembly_report:
            self._add_report(read_assembly_report(assembly_report))
        self._add_rules()

    def _add_report(self, rows):
        for r in rows:
            target = None
            for k in ('Sequence-Name', 'GenBank-Accn', 'RefSeq-Accn', 'UCSC-style-name'):
                v = r.get(k, '')
                if v in self.lengths:
                    target = v
                    break
            if target is None:
                continue
            L = r.get('Sequence-Length', '')
            if L.isdigit() and int(L) != self.lengths[target]:
                continue  # a name whose length disagrees is never an alias
            for k in ('Sequence-Name', 'GenBank-Accn', 'RefSeq-Accn', 'UCSC-style-name'):
                v = r.get(k, '')
                if v and v != 'na':
                    self._alias.setdefault(v, target)

    def _add_rules(self):
        for c in PRIMARY:
            if c in self.lengths:
                self._alias.setdefault('chr' + c, c)
        if 'MT' in self.lengths:
            for v in ('chrM', 'M', 'chrMT'):
                self._alias.setdefault(v, 'MT')

    def contig(self, name):
        """Canonical FASTA name for a VCF contig name, or None."""
        return self._alias.get(name)

    def length(self, canon):
        return self.lengths[canon]

    def is_primary(self, canon):
        return canon in PRIMARY

    def fetch(self, canon, start1, end1):
        """Bases start1..end1 (1-based, inclusive), upper-case; raises ValueError outside the contig."""
        if start1 < 1 or end1 > self.lengths[canon] or end1 < start1 - 1:
            raise ValueError(f'{canon}:{start1}-{end1} outside contig of length {self.lengths.get(canon)}')
        if end1 < start1:
            return ''
        return self.fa.fetch(canon, start1 - 1, end1).upper()

    def check_ref(self, canon, pos1, ref):
        """('match' | 'mismatch' | 'n' | 'out_of_contig', reference_bases): 'n' when either side has an N (not verifiable)."""
        end = pos1 + len(ref) - 1
        if pos1 < 1 or end > self.lengths[canon]:
            return 'out_of_contig', ''
        seq = self.fetch(canon, pos1, end)
        if 'N' in seq or 'N' in ref.upper():
            return 'n', seq
        return ('match' if seq == ref.upper() else 'mismatch'), seq

    def assembly_check(self, header_contigs):
        """Compare VCF header contig lengths with this reference.

        header_contigs: {name: length or None}. Returns dict(ok, build_guess, matched, mismatched, unmapped, no_length).
        ok is False when any mapped contig has a different length, or when chr1's length names another build."""
        matched, mismatched, unmapped, no_length = [], [], [], []
        build = None
        for n, L in header_contigs.items():
            c = self.contig(n)
            if c is None:
                unmapped.append(n)
                continue
            if L is None:
                no_length.append(n)
                continue
            (matched if L == self.lengths[c] else mismatched).append(n if L == self.lengths[c] else (n, L, self.lengths[c]))
            if c == '1':
                build = self.builds.get(L, f'unknown(chr1={L})')
        return dict(ok=not mismatched and build in (None, self.name), build_guess=build, matched=len(matched), mismatched=mismatched,
                    unmapped=unmapped, no_length=no_length)
