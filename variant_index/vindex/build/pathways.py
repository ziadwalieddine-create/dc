"""Reactome pathway membership -> gene_set_member rows.

Ensembl2Reactome maps an Ensembl gene to the lowest-level Reactome pathways it takes part in (CC0). Membership is an annotation used to group matched
variants in a pathway view; it is not evidence about any variant, and a gene's membership says nothing about the allele a sample carries. Only human
rows with an Ensembl gene id are kept; the other species and the transcript/protein identifiers in the file are counted, not used.
"""
import os

from .common import DATA, retrieved

REACTOME = os.path.join(DATA, 'reactome')


class ReactomeAdapter:
    source_id = 'reactome'
    targets = ['reactome_ensembl2reactome']

    def __init__(self, genes=None, root=None):
        self.genes, self.root = genes, root or REACTOME

    def run(self, ctx, limit=None):
        ctx.adapter = self
        st = ctx.st()
        rel = ctx.register_release('reactome', 'reactome_ensembl2reactome', retrieved(ctx, 'reactome_ensembl2reactome'))
        known = getattr(self.genes, 'by_ensg', None)
        seen = set()
        with open(os.path.join(self.root, 'Ensembl2Reactome.txt'), encoding='utf-8', newline='') as f:
            for line in f:
                c = line.rstrip('\n').split('\t')
                st.inc('rows')
                if len(c) < 6 or c[5] != 'Homo sapiens':
                    st.inc('rows_of_other_species')
                    continue
                if not c[0].startswith('ENSG'):
                    st.inc('human_rows_not_keyed_by_gene')
                    continue
                if limit and st.c.get('memberships', 0) >= limit:
                    break
                gid = c[0].split('.')[0]
                if (c[1], gid) in seen:
                    st.inc('duplicate_memberships')
                    continue
                seen.add((c[1], gid))
                if known is not None and gid not in known:
                    st.inc('memberships_of_genes_not_in_gene_map')
                ctx.emit('gene_set_member', dict(set_id=c[1], source_id='reactome', kind='reactome_pathway', name=c[3], gene_id=gid, release_id=rel), name='reactome')
                st.inc('memberships')
        return st.c
