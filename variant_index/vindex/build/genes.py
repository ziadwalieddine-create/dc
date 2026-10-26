"""Gene identity: Ensembl gene id (version stripped) is the entity id; HGNC id, NCBI Gene id and symbols are cross-references.

The gene table is the Open Targets 26.09 target dataset (78,733 genes: approved symbol, GRCh38 location, HGNC id, previous symbols) with the
MANE Select transcript and NCBI Gene id from the MANE summary. Symbol lookups never guess: a previous symbol that belongs to more than one
current gene is treated as unresolvable."""
import glob
import gzip
import os

from .common import DATA


class GeneMap:
    def __init__(self, reference=None, ot_dir=None, mane_path=None):
        self.by_ensg, self.by_symbol, self.by_hgnc, self.by_ncbi, self.prev = {}, {}, {}, {}, {}
        ot_dir = ot_dir or os.path.join(DATA, 'opentargets', '26.09', 'target')
        mane_path = mane_path or os.path.join(DATA, 'mane', 'MANE.GRCh38.v1.5.summary.txt.gz')
        files = sorted(glob.glob(os.path.join(ot_dir, '*.parquet')))
        if files:
            import pyarrow.parquet as pq
            for f in files:
                for r in pq.read_table(f, columns=['id', 'approvedSymbol', 'biotype', 'genomicLocation', 'dbXrefs', 'obsoleteSymbols']).to_pylist():
                    loc = r['genomicLocation'] or {}
                    hg = next((x['id'] for x in (r['dbXrefs'] or []) if x['source'] == 'HGNC'), None)
                    g = dict(gene_id=r['id'], symbol=r['approvedSymbol'], hgnc_id=f'HGNC:{hg}' if hg else '', chrom=str(loc.get('chromosome') or ''), start=loc.get('start'), end=loc.get('end'),
                             strand={1: '+', -1: '-'}.get(loc.get('strand'), ''), mane_select='', biotype=r['biotype'] or '', ncbi='')
                    self.by_ensg[g['gene_id']] = g
                    if g['symbol']:
                        self.by_symbol.setdefault(g['symbol'], g)
                    if hg:
                        self.by_hgnc[g['hgnc_id']] = g
                    for ob in (r['obsoleteSymbols'] or []):
                        s = ob['label'] if isinstance(ob, dict) else ob
                        self.prev.setdefault(s, set()).add(g['gene_id'])
        if os.path.exists(mane_path):
            with gzip.open(mane_path, 'rt', encoding='utf-8') as f:
                hdr = f.readline().rstrip('\n').lstrip('#').split('\t')
                for line in f:
                    r = dict(zip(hdr, line.rstrip('\n').split('\t')))
                    ensg = r['Ensembl_Gene'].split('.')[0]
                    g = self.by_ensg.get(ensg)
                    if g is None:
                        g = dict(gene_id=ensg, symbol=r['symbol'], hgnc_id=r['HGNC_ID'], chrom=(reference.contig(r['GRCh38_chr']) if reference else '') or r['GRCh38_chr'], start=int(r['chr_start']),
                                 end=int(r['chr_end']), strand=r['chr_strand'], mane_select='', biotype='', ncbi='')
                        self.by_ensg[ensg] = g
                        self.by_symbol.setdefault(g['symbol'], g)
                    g['ncbi'] = r['NCBI_GeneID'].replace('GeneID:', '')
                    self.by_ncbi[g['ncbi']] = g
                    if r['MANE_status'] == 'MANE Select':
                        g['mane_select'] = r['Ensembl_nuc']
                    if not self.by_hgnc.get(r['HGNC_ID']):
                        self.by_hgnc[r['HGNC_ID']] = g
        self.rows = list(self.by_ensg.values())

    def from_ncbi(self, ncbi_id):
        g = self.by_ncbi.get(str(ncbi_id))
        return g['gene_id'] if g else None

    def from_hgnc(self, hgnc):
        g = self.by_hgnc.get(hgnc if str(hgnc).startswith('HGNC:') else f'HGNC:{hgnc}')
        return g['gene_id'] if g else None

    def from_symbol(self, symbol):
        """Approved symbol, else a previous symbol that identifies exactly one current gene; None otherwise."""
        g = self.by_symbol.get(symbol)
        if g:
            return g['gene_id']
        p = self.prev.get(symbol)
        return next(iter(p)) if p and len(p) == 1 else None

    def resolve(self, hgnc=None, symbol=None, ncbi=None):
        """(gene_id or None, how). HGNC id first, then NCBI Gene id, then symbol."""
        for key, fn in (('hgnc', self.from_hgnc), ('ncbi', self.from_ncbi), ('symbol', self.from_symbol)):
            v = {'hgnc': hgnc, 'ncbi': ncbi, 'symbol': symbol}[key]
            if v:
                g = fn(v)
                if g:
                    return g, key
        return None, 'unresolved'


class GeneAdapter:
    """Writes the gene entity table (all genes of the Open Targets target dataset) and symbol/HGNC cross-references."""
    source_id = 'opentargets'
    targets = ['opentargets_26_09_target', 'mane_summary_v1_5']

    def __init__(self, genes):
        self.genes = genes

    def run(self, ctx, limit=None):
        ctx.adapter = self
        rel = ctx.register_release('opentargets', 'opentargets_26_09_target', '26.09')
        st = ctx.st()
        for g in sorted(self.genes.rows, key=lambda g: g['gene_id']):
            ctx.emit('entity_gene', dict(gene_id=g['gene_id'], symbol=g['symbol'], hgnc_id=g['hgnc_id'], chrom=g['chrom'], start=g['start'], end=g['end'], strand=g['strand'],
                                         mane_select=g['mane_select'], biotype=g['biotype'], release_id=rel), name='genes')
            if g['hgnc_id']:
                ctx.emit('xref', dict(namespace='hgnc', identifier=g['hgnc_id'], entity_level='gene', entity_id=g['gene_id'], release_id=rel, relation='open_targets_dbXrefs'), name='genes')
            if g['ncbi']:
                ctx.emit('xref', dict(namespace='ncbi_gene', identifier=g['ncbi'], entity_level='gene', entity_id=g['gene_id'], release_id=rel, relation='mane_summary'), name='genes')
            st.inc('genes')
        return st.c
