"""Gene identity helper: MANE summary -> Ensembl gene id, NCBI Gene id, HGNC id, symbol, MANE Select transcript."""
import gzip
import os

from .common import DATA


class GeneMap:
    def __init__(self, reference=None, path=None):
        path = path or os.path.join(DATA, 'mane', 'MANE.GRCh38.v1.5.summary.txt.gz')
        self.by_ncbi, self.by_symbol, self.by_ensg, self.rows = {}, {}, {}, []
        if not os.path.exists(path):
            return
        with gzip.open(path, 'rt', encoding='utf-8') as f:
            hdr = f.readline().rstrip('\n').lstrip('#').split('\t')
            for line in f:
                r = dict(zip(hdr, line.rstrip('\n').split('\t')))
                ensg = r['Ensembl_Gene'].split('.')[0]
                g = dict(gene_id=ensg, symbol=r['symbol'], hgnc_id=r['HGNC_ID'], ncbi=r['NCBI_GeneID'], chrom=(reference.contig(r['GRCh38_chr']) if reference else None) or r['GRCh38_chr'],
                         start=int(r['chr_start']), end=int(r['chr_end']), strand=r['chr_strand'], mane_select=r['Ensembl_nuc'] if r['MANE_status'] == 'MANE Select' else '')
                if r['MANE_status'] == 'MANE Select' or g['gene_id'] not in self.by_ensg:
                    self.by_ensg[g['gene_id']] = g
                self.by_ncbi.setdefault(g['ncbi'].replace('GeneID:', ''), g)
                self.by_symbol.setdefault(g['symbol'], g)
        self.rows = list(self.by_ensg.values())

    def from_ncbi(self, ncbi_id):
        g = self.by_ncbi.get(str(ncbi_id))
        return g['gene_id'] if g else None

    def from_symbol(self, symbol):
        g = self.by_symbol.get(symbol)
        return g['gene_id'] if g else None
