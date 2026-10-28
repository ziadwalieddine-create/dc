"""Synthetic alignments for the read-evidence tests: reads are cut from the synthetic reference with the edits a test states, and written as a coordinate-sorted, indexed BAM.

Nothing here is human data. A read is described by its operations on the reference, so that the same allele can be written with the gap in different places:

    ('M', n)     copy n reference bases            ('X', 'G')   one base that differs from the reference (aligned, CIGAR M)
    ('D', n)     skip n reference bases            ('I', 'AC')  inserted bases              ('S', 'AC')  soft-clipped bases        ('N', n) reference skip
"""
import os
from dataclasses import dataclass, field

import pysam

FLAG_PAIRED, FLAG_PROPER, FLAG_UNMAPPED, FLAG_REVERSE, FLAG_MATE_REVERSE, FLAG_READ1, FLAG_READ2 = 1, 2, 4, 16, 32, 64, 128
FLAG_SECONDARY, FLAG_QCFAIL, FLAG_DUPLICATE, FLAG_SUPPLEMENTARY = 256, 512, 1024, 2048


@dataclass
class Read:
    name: str
    chrom: str
    start1: int
    seq: str
    cigar: str
    flag: int = 0
    mapq: int = 60
    quals: object = 35                     # a base quality for every base, a list, or None for none stored
    tags: dict = field(default_factory=dict)


def make_read(fa, chrom, start1, ops, name, flag=0, mapq=60, quals=35):
    seq, cig, r = [], [], start1
    for op, arg in ops:
        if op == 'M':
            seq.append(fa.fetch(chrom, r - 1, r - 1 + arg))
            r += arg
            n = arg
        elif op == 'X':
            seq.append(arg)
            r += 1
            op, n = 'M', 1
        elif op in ('D', 'N'):
            r += arg
            n = arg
        elif op in ('I', 'S'):
            seq.append(arg)
            n = len(arg)
        else:
            raise ValueError(op)
        if cig and cig[-1][0] == op:
            cig[-1][1] += n
        else:
            cig.append([op, n])
    return Read(name, chrom, start1, ''.join(seq), ''.join(f'{n}{op}' for op, n in cig), flag, mapq, quals)


def span(fa, chrom, start1, end1, name, edits=(), **kw):
    """A read that copies the reference over start1..end1 (inclusive) except where `edits` {pos1: base} substitute a base."""
    ops, r = [], start1
    for pos, base in sorted(dict(edits).items()):
        if pos < start1 or pos > end1:
            continue
        if pos > r:
            ops.append(('M', pos - r))
        ops.append(('X', base))
        r = pos + 1
    if r <= end1:
        ops.append(('M', end1 - r + 1))
    return make_read(fa, chrom, start1, ops, name, **kw)


def write_bam(path, reads, contigs, rename=None):
    """A coordinate-sorted, indexed BAM of `reads` on `contigs` [(name, length)]; `rename` {contig: name in the header} makes the header name differ from the reference's."""
    rename = rename or {}
    order = {c: i for i, (c, _) in enumerate(contigs)}
    header = {'HD': {'VN': '1.6', 'SO': 'coordinate'}, 'SQ': [{'SN': rename.get(c, c), 'LN': n} for c, n in contigs]}
    with pysam.AlignmentFile(path, 'wb', header=header) as out:
        for r in sorted(reads, key=lambda r: (order[r.chrom], r.start1)):
            a = pysam.AlignedSegment(out.header)
            a.query_name = r.name
            a.flag = r.flag
            a.reference_id = order[r.chrom]
            a.reference_start = r.start1 - 1
            a.mapping_quality = r.mapq
            if not r.flag & FLAG_UNMAPPED:
                a.cigarstring = r.cigar
            a.query_sequence = r.seq
            if r.quals is None:
                a.query_qualities = None
            elif isinstance(r.quals, int):
                a.query_qualities = pysam.qualitystring_to_array(chr(33 + r.quals) * len(r.seq))
            else:
                a.query_qualities = pysam.qualitystring_to_array(''.join(chr(33 + q) for q in r.quals))
            for k, v in r.tags.items():
                a.set_tag(k, v)
            out.write(a)
    pysam.index(path)
    return path
