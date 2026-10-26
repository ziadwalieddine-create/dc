"""Deterministic synthetic reference ('SYN38') for the tests: small contigs with engineered repeats, so that left-alignment,
equivalent representations, strand-ambiguous sites and star-allele loci can be tested without any real genome.

The sequence is pseudo-random (fixed seed) apart from the engineered features listed in FEATURES. Nothing here is human data."""
import hashlib
import os
import random
import sys

CONTIGS = {'1': 3000, '2': 1500, 'X': 1500, 'MT': 400}
BUILD = {3000: 'SYN38'}

# (contig, 1-based start, sequence) overrides
FEATURES = [
    ('1', 100, 'GCAAAAAAAAGT'),            # A8 homopolymer, anchored by C (101) and G (110..)
    ('1', 200, 'TCAGCAGCAGCAGCAGTT'),       # (CAG)5 tandem repeat
    ('1', 300, 'CGTGTGTGTGTGAC'),           # (GT)5 dinucleotide repeat
    ('1', 400, 'ACGTACGTAC'),               # plain mixed sequence for SNVs
    ('1', 500, 'GATTACAGATTACA'),           # adjacent bases for MNV / phase tests
    ('1', 600, 'AATTCCGG'),                 # A/T and C/G (palindromic) SNV sites
    ('1', 700, 'CCATGGCCTTAGGCATCC'),       # 'gene' with defining variants at 705, 710, 715 (pharmacogene-like)
    ('2', 100, 'TTTTGAAAAC'),
]


def build(seed=20261026):
    rng = random.Random(seed)
    seqs = {}
    for c, L in CONTIGS.items():
        seqs[c] = [rng.choice('ACGT') for _ in range(L)]
    for c, start, s in FEATURES:
        for i, b in enumerate(s):
            seqs[c][start - 1 + i] = b
    return {c: ''.join(v) for c, v in seqs.items()}


def write(path):
    seqs = build()
    with open(path, 'w') as f:
        for c, s in seqs.items():
            f.write(f'>{c}\n')
            for i in range(0, len(s), 60):
                f.write(s[i:i + 60] + '\n')
    return hashlib.sha256(open(path, 'rb').read()).hexdigest()


if __name__ == '__main__':
    out = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.abspath(__file__)), 'syn38.fa')
    print(write(out), out)
