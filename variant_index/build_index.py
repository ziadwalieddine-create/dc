#!/usr/bin/env python3
"""Build the High-Value Human Variant Reference Index from the acquired sources.

    python3 variant_index/build_index.py [--out data/index] [--sources clinvar,clingen,...] [--limit N]

Every build starts from an empty output directory and runs the selected adapters in dependency order, so the same sources and code give the
same table digests (table_digests.json). The index stays under data/ and is never committed: it contains only public-source records, but it is
large and reproducible. Personal calls are matched against it by vindex.match_* and are never written into it.
"""
import argparse
import json
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from vindex.build import runner                                   # noqa: E402
from vindex.build.common import DATA, OUT_DEFAULT                 # noqa: E402
from vindex.refgenome import Reference                            # noqa: E402

FASTA = os.path.join(DATA, 'reference_sequence', 'Homo_sapiens.GRCh38.dna.toplevel.fa.gz')
REPORT = os.path.join(DATA, 'refseq', 'GCF_000001405.40_GRCh38.p14_assembly_report.txt')

ORDER = ['genes', 'clinvar', 'clingen_erepo', 'clingen_gene_region', 'gwas', 'cpic', 'pharmvar', 'clinpgx']


def make_adapters(names, ref):
    from vindex.build.genes import GeneMap, GeneAdapter
    from vindex.build.clinvar_vcf import ClinVarVcfAdapter
    from vindex.build.clingen import ClinGenErepoAdapter, ClinGenGeneRegionAdapter
    from vindex.build.gwas import GwasCatalogAdapter
    from vindex.build.alleles import OtVariantIndex
    from vindex.build.pgx import CpicAdapter, PharmVarAdapter
    from vindex.build.clinpgx import ClinPgxAdapter
    genes = GeneMap(ref)
    factories = {
        'genes': lambda: GeneAdapter(genes),
        'clinvar': lambda: ClinVarVcfAdapter(genes),
        'clingen_erepo': lambda: ClinGenErepoAdapter(),
        'clingen_gene_region': lambda: ClinGenGeneRegionAdapter(genes),
        'gwas': lambda: GwasCatalogAdapter(genes, OtVariantIndex(os.path.join(DATA, 'opentargets', '26.09', 'variant'))),
        'cpic': lambda: CpicAdapter(genes),
        'pharmvar': lambda: PharmVarAdapter(),
        'clinpgx': lambda: ClinPgxAdapter(genes),
    }
    unknown = [n for n in names if n not in factories]
    if unknown:
        raise SystemExit(f'unknown source(s) {unknown}; available: {", ".join(ORDER)}')
    return [factories[n]() for n in ORDER if n in names]


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--out', default=OUT_DEFAULT)
    ap.add_argument('--sources', default=','.join(ORDER), help='comma-separated subset of: ' + ', '.join(ORDER))
    ap.add_argument('--limit', type=int, default=None, help='stop each adapter after this many source records (for trials)')
    ap.add_argument('--fasta', default=FASTA)
    ap.add_argument('--assembly-report', default=REPORT)
    a = ap.parse_args(argv)
    t0 = time.time()
    ref = Reference(a.fasta, a.assembly_report)
    names = [s for s in a.sources.split(',') if s]
    man = runner.build(a.out, make_adapters(names, ref), ref, limit=a.limit)
    print(json.dumps(dict(out=a.out, seconds=round(time.time() - t0, 1), tables=man['tables'], adapters={x['adapter']: x['stats'] for x in man['adapters']}, problems=man['problems']), indent=1))
    return 0


if __name__ == '__main__':
    sys.exit(main())
