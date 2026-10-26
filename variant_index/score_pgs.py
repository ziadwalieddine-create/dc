#!/usr/bin/env python3
"""Score a sample with PGS Catalog models, one chromosome at a time, with the registry's gates beside every result.

    python3 variant_index/score_pgs.py --vcf SAMPLE.vcf.gz --out data/results/SAMPLE/pgs (--pgs PGS000001,PGS000002 | --trait EFO_0004339)
                                       [--sample-ancestry European] [--reference-distributions dist.json] [--callable-bed BED] [--per-weight]

--trait selects the acquired scores whose mapped EFO ids include the given id. --reference-distributions is {pgs_id: {ancestry group: {"mean": m, "sd": s}}} computed
from a reference panel; without it every score stays uncalibrated. Results go under data/ (git-ignored) or outside the repository.
"""
import argparse
import json
import os
import sys

import pyarrow.parquet as pq

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from vindex import pgs_run                                         # noqa: E402
from vindex.build.common import DATA, REPO                          # noqa: E402
from vindex.refgenome import Reference                              # noqa: E402
from vindex.vcfio import QCPolicy                                    # noqa: E402

FASTA = os.path.join(DATA, 'reference_sequence', 'Homo_sapiens.GRCh38.dna.toplevel.fa.gz')
REPORT = os.path.join(DATA, 'refseq', 'GCF_000001405.40_GRCh38.p14_assembly_report.txt')


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--vcf', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--registry', default=os.path.join(REPO, 'data', 'pgs_registry'))
    ap.add_argument('--scoring-dir', default=os.path.join(DATA, 'pgs_catalog', 'scores_grch38'))
    ap.add_argument('--pgs', help='comma-separated PGS ids')
    ap.add_argument('--trait', help='EFO/MONDO id: every acquired score mapped to it')
    ap.add_argument('--sample')
    ap.add_argument('--sample-ancestry', help='genetically inferred ancestry group of the sample, as the catalogue names groups (for example European); unassessed when omitted')
    ap.add_argument('--reference-distributions')
    ap.add_argument('--callable-bed')
    ap.add_argument('--per-weight', action='store_true')
    ap.add_argument('--weights-budget', type=int, default=4_000_000, help='most scoring-file weights held in memory at once; more scores are scored in further passes over the callset')
    ap.add_argument('--min-gq', type=int, default=20)
    ap.add_argument('--min-dp', type=int, default=10)
    ap.add_argument('--fasta', default=FASTA)
    ap.add_argument('--assembly-report', default=REPORT)
    a = ap.parse_args(argv)
    models = {r['pgs_id']: r for r in pq.read_table(os.path.join(a.registry, 'pgs_model.parquet')).to_pylist()}
    ids = [x for x in (a.pgs or '').split(',') if x]
    if a.trait:
        ids += sorted(i for i, m in models.items() if a.trait in (m['efo_ids'] or []) and m['scoring_file'] == 'acquired')
    ids = sorted(set(ids))
    if not ids:
        raise SystemExit('no scores selected: give --pgs or a --trait that matches acquired scores')
    dist = json.load(open(a.reference_distributions)) if a.reference_distributions else None
    man = pgs_run.score_models(a.vcf, Reference(a.fasta, a.assembly_report), a.scoring_dir, models, ids, a.out, sample=a.sample, qc=QCPolicy(min_gq=a.min_gq, min_dp=a.min_dp),
                               sample_ancestry=a.sample_ancestry, reference_distributions=dist, callable_bed=a.callable_bed, per_weight=a.per_weight, weights_budget=a.weights_budget)
    print(json.dumps(man, indent=1))
    return 0


if __name__ == '__main__':
    sys.exit(main())
