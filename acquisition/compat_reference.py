"""Tests that tie the mirrored GRCh38 reference sequence and the Ensembl release-116 gene annotation to the other coordinate-bearing
files (called from compat_tests.py; needs pysam, the FASTA and the GTF).

  * the remote range reads (fetch_reference_regional.py) equal the mirrored FASTA
  * the REF allele of every ClinVar record and of every AlphaMissense position equals the reference sequence
  * the start and stop codons of the annotation read ATG and TAA/TAG/TGA in the reference sequence
  * MANE transcripts, Open Targets gene locations and the MaveDB gene set resolve in the annotation
  * the amino-acid change implied by the reference sequence and the annotation's coding exons equals AlphaMissense's protein_variant

These are consistency tests on files, not an annotation engine. The small codon helper below exists only to run the last test.
A test is PASS or FAIL against a stated criterion, or INFO where the result describes a relationship; a test that raises is
recorded as FAIL with the error.
"""
import collections
import csv
import glob
import gzip
import json
import os
import re
import zipfile

from mavedb_snvs import NC2CHR

FASTA_ID, GTF_ID, LOF_ID, PV_ID = 'ensembl_grch38_dna_toplevel', 'ensembl_grch38_gtf_116', 'gnomad_v4_1_lof_curation_incomplete_penetrance', 'pharmvar_6_2_29'
_BASES = 'TCAG'
_AA = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
CODON = {a + b + c: _AA[16 * i + 4 * j + k] for i, a in enumerate(_BASES) for j, b in enumerate(_BASES) for k, c in enumerate(_BASES)}
_COMP = str.maketrans('ACGTN', 'TGCAN')
GID, TID, TVER, TAG, TBIO, GNAME = (re.compile(x) for x in (r'gene_id "([^"]+)"', r'transcript_id "([^"]+)"', r'transcript_version "([^"]+)"', r'tag "([^"]+)"', r'transcript_biotype "([^"]+)"', r'gene_name "([^"]+)"'))
CLNSIG = re.compile(r'CLNSIG=([^;]+)')
PV = re.compile(r'^([A-Z*])(\d+)([A-Z*])$')


def rc(s):
    return s.translate(_COMP)[::-1]


def parse_gtf(p):
    """Genes (with names), transcripts, CDS pieces and codon pieces of the annotation; exon and UTR lines are not needed here."""
    genes, names, tx = {}, {}, {}
    cds = collections.defaultdict(list)
    codons = {'start_codon': collections.defaultdict(list), 'stop_codon': collections.defaultdict(list)}
    with gzip.open(p, 'rt') as f:
        for line in f:
            if line[0] == '#':
                continue
            c = line.split('\t')
            feat = c[2]
            if feat in ('exon', 'five_prime_utr', 'three_prime_utr', 'Selenocysteine'):
                continue
            attr = c[8]
            if feat == 'gene':
                gid = GID.search(attr).group(1)
                genes[gid] = (c[0], int(c[3]), int(c[4]), c[6])
                gn = GNAME.search(attr)
                names[gid] = gn.group(1) if gn else None
            elif feat == 'transcript':
                v, b = TVER.search(attr), TBIO.search(attr)
                tx[TID.search(attr).group(1)] = dict(chrom=c[0], start=int(c[3]), end=int(c[4]), strand=c[6], version=v.group(1) if v else None, tags=set(TAG.findall(attr)), biotype=b.group(1) if b else None)
            elif feat == 'CDS':
                cds[TID.search(attr).group(1)].append((int(c[3]), int(c[4])))
            elif feat in codons:
                codons[feat][TID.search(attr).group(1)].append((c[0], int(c[3]), int(c[4])))
    return genes, names, tx, cds, codons


def spliced(fa, chrom, pieces, strand):
    """Sequence of genomic pieces [(start, end) 1-based inclusive] read on the transcript strand."""
    seq = ''.join(fa.fetch(chrom, s - 1, e) for s, e in sorted(pieces))
    return seq if strand == '+' else rc(seq)


def layout(pieces, strand):
    """CDS pieces in transcription order with the number of coding bases before each."""
    out, cum = [], 0
    for s, e in sorted(pieces, reverse=strand == '-'):
        out.append((s, e, cum))
        cum += e - s + 1
    return out, cum


def implied_change(fa, chrom, strand, pieces, pos, ref, alt, cache, key):
    """Amino-acid change a single-nucleotide change at genomic `pos` causes in a coding sequence: (position, ref aa, alt aa, reference base agrees) or None if outside the CDS."""
    if key not in cache:
        cache[key] = (layout(pieces, strand)[0], spliced(fa, chrom, pieces, strand))
    lay, seq = cache[key]
    for s, e, cum in lay:
        if s <= pos <= e:
            cds_pos = cum + (pos - s if strand == '+' else e - pos)
            break
    else:
        return None
    idx, off = divmod(cds_pos, 3)
    codon = seq[3 * idx:3 * idx + 3]
    if len(codon) < 3:
        return None
    b_ref, b_alt = (ref, alt) if strand == '+' else (ref.translate(_COMP), alt.translate(_COMP))
    new = codon[:off] + b_alt + codon[off + 1:]
    return idx + 1, CODON.get(codon, '?'), CODON.get(new, '?'), codon[off] == b_ref


def guarded(src, sid, name, fn):
    """Run one group of tests; an exception is a recorded FAIL, never a silent skip."""
    try:
        fn()
    except Exception as e:  # noqa: BLE001
        src(sid, name + ' (raised an error)', False, error=repr(e)[:300])


def run(ctx):
    src, log, D, MAN = ctx['src'], ctx['log'], ctx['D'], ctx['MAN']
    fa_path = os.path.join(D, 'reference_sequence', 'Homo_sapiens.GRCh38.dna.toplevel.fa.gz')
    gtf_path = os.path.join(D, 'ensembl_annotation', 'Homo_sapiens.GRCh38.116.gtf.gz')
    if not (MAN.get(FASTA_ID, {}).get('ok') and MAN.get(GTF_ID, {}).get('ok') and os.path.exists(fa_path + '.gzi') and os.path.exists(gtf_path)):
        src('ensembl_genomes', 'Reference FASTA and annotation tests not run: the FASTA or the GTF is not acquired', None)
        return
    try:
        import pysam
    except ImportError:
        src('ensembl_genomes', 'Reference FASTA and annotation tests not run: pysam is not installed', None, note='pip install -r acquisition/requirements.txt')
        return
    fa = pysam.FastaFile(fa_path)
    refs = set(fa.references)
    state = {}
    lof_fp = os.path.join(D, 'gnomad', 'lof_curation', 'incomplete_penetrance_curation_results.csv')
    if MAN.get(LOF_ID, {}).get('ok') and os.path.exists(lof_fp):
        with open(lof_fp, newline='', encoding='utf-8') as f:
            state['lof_rows'] = list(csv.DictReader(f))
        state['lof_keys'] = {}
        for r in state['lof_rows']:
            c_, p_, ref_, alt_ = r['Variant ID'].split('-')
            state['lof_keys'][(c_, int(p_), ref_, alt_)] = r['Verdict']
        state['lof_clinvar'] = {}

    def t_slice():
        """The remote range reads equal the mirrored file."""
        slice_fp = os.path.join(D, 'reference_sequence', 'ensembl_grch38_toplevel_at_mavedb_regions.tsv.gz')
        if not (MAN.get('ensembl_grch38_reference_regional', {}).get('ok') and os.path.exists(slice_fp)):
            return
        n = equal = 0
        with gzip.open(slice_fp, 'rt', newline='') as f:
            for r in csv.DictReader(f, delimiter='\t'):
                n += 1
                equal += fa.fetch(r['chrom'], int(r['start']) - 1, int(r['end'])) == r['sequence']
        src('ensembl_genomes', 'Remote range reads of the reference equal the mirrored FASTA', n > 1000 and equal == n, targets=[FASTA_ID, 'ensembl_grch38_reference_regional'], regions=n, equal=equal,
            note='two independent routes to the same release-116 file: HTTP range reads through htslib on the server copy, and the checksum-verified local copy')

    def t_clinvar():
        """REF allele of every ClinVar record."""
        log('reference: ClinVar REF alleles ...')
        lof_keys = state.get('lof_keys', {})
        checked, match, unknown, bad = collections.Counter(), collections.Counter(), 0, []
        with gzip.open(os.path.join(D, 'clinvar', 'clinvar.vcf.gz'), 'rt', errors='replace') as f:
            for line in f:
                if line[0] == '#':
                    continue
                c = line.split('\t', 8)
                chrom, pos, ref = c[0], int(c[1]), c[3]
                if chrom not in refs:
                    unknown += 1
                    continue
                if (chrom, pos, ref, c[4]) in lof_keys:
                    m_ = CLNSIG.search(c[7])
                    state['lof_clinvar'][(chrom, pos, ref, c[4])] = m_.group(1) if m_ else 'NA'
                kind = 'snv' if len(ref) == 1 and len(c[4]) == 1 else 'other'
                checked[kind] += 1
                if fa.fetch(chrom, pos - 1, pos - 1 + len(ref)).upper() == ref.upper():
                    match[kind] += 1
                elif len(bad) < 5:
                    bad.append((chrom, pos, ref, c[4]))
        tot, ok_ = sum(checked.values()), sum(match.values())
        src('clinvar', 'REF allele of every ClinVar record equals the GRCh38 reference sequence', tot > 100000 and ok_ / tot >= 0.9999, targets=['clinvar_vcf_grch38', FASTA_ID], records=tot, equal=ok_,
            share=round(ok_ / max(1, tot), 6), snv=dict(checked=checked['snv'], equal=match['snv']), other=dict(checked=checked['other'], equal=match['other']), not_tested_contig_named_differently_in_fasta=unknown,
            mismatch_examples=bad, note='all records are tested, indels and MNVs included; the untested records lie on RefSeq NT_/NW_ contigs, which the Ensembl FASTA names differently')

    def t_annotation():
        """Codons, MANE, Open Targets and MaveDB genes against the annotation."""
        log('reference: parsing the GTF ...')
        genes, names, tx, cds, codons = parse_gtf(gtf_path)
        state.update(tx=tx, cds=cds, genes=genes, names=names)
        for feat, want, name in (('start_codon', {'ATG'}, 'start'), ('stop_codon', {'TAA', 'TAG', 'TGA'}, 'stop')):
            seen = collections.Counter()
            for t, pieces in codons[feat].items():
                info = tx[t]
                seen[spliced(fa, info['chrom'], [(s, e) for _, s, e in pieces], info['strand']).upper()] += 1
            tot = sum(seen.values())
            ok_ = sum(v for k, v in seen.items() if k in want)
            src('ensembl_genomes', f'{name.capitalize()} codons of the Ensembl annotation read {"/".join(sorted(want))} in the reference sequence', tot > 50000 and ok_ / tot >= 0.99, targets=[GTF_ID, FASTA_ID],
                transcripts=tot, as_expected=ok_, share=round(ok_ / max(1, tot), 5), other_codons=dict(seen.most_common(8)) if len(seen) > len(want) else {k: v for k, v in seen.items() if k not in want},
                note='a coordinate, strand or assembly mismatch between annotation and sequence would make these shares collapse')
        mane = []
        with gzip.open(os.path.join(D, 'mane', 'MANE.GRCh38.v1.5.summary.txt.gz'), 'rt') as f:
            mane = list(csv.DictReader(f, delimiter='\t'))
        present = same_version = tagged = same_span = n_sel = off_chromosome = 0
        for r in mane:
            t, _, ver = r['Ensembl_nuc'].partition('.')
            info = tx.get(t)
            if info is None:
                off_chromosome += r['GRCh38_chr'][:9] not in NC2CHR       # a fix-patch or alternate-locus sequence rather than a chromosome
                continue
            present += 1
            same_version += info['version'] == ver
            if r['MANE_status'] == 'MANE Select':
                n_sel += 1
                tagged += 'MANE_Select' in info['tags']
            same_span += (info['chrom'], info['start'], info['end']) == (NC2CHR.get(r['GRCh38_chr'][:9]), int(r['chr_start']), int(r['chr_end']))
        src('gencode_mane', 'MANE v1.5 transcripts resolve in the Ensembl release-116 annotation', len(mane) > 18000 and present / len(mane) >= 0.95, targets=['mane_summary_v1_5', GTF_ID], mane_transcripts=len(mane),
            present_by_stable_id=present, same_version=same_version, mane_select=n_sel, mane_select_tagged_in_gtf=tagged, same_genomic_span=same_span, not_found=len(mane) - present, not_found_on_patch_or_alternate_sequences=off_chromosome,
            note='the MANE transcripts not found lie on fix-patch and alternate-locus sequences, which the default Ensembl GTF does not annotate; for the others the version, the span and the MANE_Select tag are compared')
        try:
            import pyarrow.parquet as pq
        except ImportError:
            pq = None
        if pq is not None:
            n_ot = found = same = 0
            for fp in sorted(glob.glob(os.path.join(D, 'opentargets', '26.09', 'target', '*.parquet'))):
                for r in pq.read_table(fp, columns=['id', 'genomicLocation']).to_pylist():
                    g = r['genomicLocation']
                    if not g:
                        continue
                    n_ot += 1
                    if r['id'] in genes:
                        found += 1
                        same += genes[r['id']] == (g['chromosome'], g['start'], g['end'], '+' if g['strand'] == 1 else '-')
            src('opentargets', 'Open Targets target locations agree with the Ensembl release-116 annotation', n_ot > 50000 and found / n_ot >= 0.9, targets=['opentargets_26_09_target', GTF_ID], targets_with_location=n_ot,
                found_by_gene_id=found, same_chromosome_start_end_strand=same, share_of_found=round(same / max(1, found), 4))
        M = ctx['M_ensg_all']
        hit = sum(1 for g in M if g in genes)
        src('mavedb', 'MaveDB target genes (Ensembl IDs) resolve in the Ensembl release-116 annotation', len(M) > 400 and hit / len(M) >= 0.9, targets=['mavedb_zenodo_dump_v4', GTF_ID], mavedb_genes=len(M), found=hit,
            share=round(hit / max(1, len(M)), 4))

    def t_alphamissense():
        """REF alleles at every AlphaMissense position, and the amino-acid change implied by reference + annotation for a fixed sample of rows."""
        log('reference: AlphaMissense positions ...')
        checked = match = unknown = n_pos = 0
        bad, samples = [], []
        cur, seq, last = None, None, None
        with gzip.open(os.path.join(D, 'alphamissense', 'AlphaMissense_hg38.tsv.gz'), 'rt') as f:
            for line in f:
                if line[0] == '#':
                    continue
                i = line.find('\t')
                j = line.find('\t', i + 1)
                key = line[:j]
                if key == last:                                     # up to three ALT rows share a position and its REF
                    continue
                last = key
                name = 'MT' if line[:i] == 'chrM' else line[3:i]
                if name != cur:
                    cur = name
                    seq = fa.fetch(name) if name in refs else None
                n_pos += 1
                if seq is None:
                    unknown += 1
                    continue
                pos = int(line[i + 1:j])
                ref = line[j + 1:line.find('\t', j + 1)]
                checked += 1
                if seq[pos - 1] == ref:
                    match += 1
                elif len(bad) < 5:
                    bad.append((name, pos, ref, seq[pos - 1]))
                if n_pos % 400 == 0:
                    samples.append(line.rstrip('\n').split('\t'))
        seq = None
        src('alphamissense', 'REF allele at every AlphaMissense position equals the GRCh38 reference sequence', checked > 1000000 and match == checked, targets=['alphamissense_hg38', FASTA_ID], positions=checked,
            equal=match, share=round(match / max(1, checked), 6), not_tested_contig_not_in_fasta=unknown, mismatch_examples=bad,
            note='one check per distinct position (up to three ALT rows share a position and its REF); chrM is compared with the Ensembl MT sequence')
        tx, cds = state['tx'], state['cds']
        cache = {}
        st = collections.Counter()
        miss = []
        for chrom, pos, ref, alt, _g, _u, trans, pv, _s, _c in samples:
            t, _, ver = trans.partition('.')
            info = tx.get(t)
            if info is None:
                st['transcript_not_in_annotation'] += 1
            elif info['version'] != ver:
                st['transcript_version_differs'] += 1
            elif info['chrom'] == 'MT':
                st['mitochondrial_code_not_tested'] += 1
            elif t not in cds or 'cds_start_NF' in info['tags']:
                st['no_complete_cds'] += 1
            else:
                res = implied_change(fa, info['chrom'], info['strand'], cds[t], int(pos), ref, alt, cache, t)
                st['tested'] += 1
                if res is None:
                    st['position_outside_cds'] += 1
                    continue
                m = PV.match(pv)
                p_, a1, a2, base_ok = res
                agree = bool(m) and (int(m.group(2)), m.group(1), m.group(3)) == (p_, a1, a2)
                st['reference_base_agrees'] += base_ok
                st['change_agrees'] += agree
                if not agree and len(miss) < 8:
                    miss.append((trans, chrom, pos, ref, alt, pv, f'{a1}{p_}{a2}'))
        tested = st['tested']
        src('alphamissense', 'Amino-acid change implied by reference sequence and annotation equals AlphaMissense protein_variant', tested > 10000 and st['change_agrees'] / tested >= 0.999,
            targets=['alphamissense_hg38', GTF_ID, FASTA_ID], rows_sampled=len(samples), **dict(st), share_agreeing=round(st['change_agrees'] / max(1, tested), 5), disagreement_examples=miss,
            note='fixed sample (every 400th AlphaMissense position); only rows whose transcript has the same version in release 116 are tested; the codon helper applies the standard genetic code to the spliced CDS')

    def t_lof():
        """gnomAD LoF curation: REF alleles, genes and gene spans in the annotation, and ClinVar classes by curation verdict."""
        rows, keys = state['lof_rows'], state['lof_keys']
        genes, names = state['genes'], state['names']
        ref_ok = known = in_span = same_name = 0
        for r in rows:
            chrom, pos, ref, alt = r['Variant ID'].split('-')
            pos = int(pos)
            ref_ok += fa.fetch(chrom, pos - 1, pos - 1 + len(ref)).upper() == ref
            gid, _, sym = r['Gene'].partition(':')
            g = genes.get(gid)
            if g:
                known += 1
                in_span += g[0] == chrom and g[1] <= pos <= g[2]
                same_name += names.get(gid) == sym
        by_verdict = collections.defaultdict(collections.Counter)
        for k, verdict in keys.items():
            sig = state['lof_clinvar'].get(k)
            cls = 'not in ClinVar' if sig is None else ('pathogenic or likely pathogenic' if 'athogenic' in sig and 'onflicting' not in sig else ('benign or likely benign' if 'enign' in sig and 'onflicting' not in sig else 'other classification'))
            by_verdict[verdict][cls] += 1
        n = len(rows)
        src('gnomad_constraint', 'gnomAD v4.1 LoF curation: REF alleles equal the reference and genes resolve in the Ensembl annotation', n > 500 and ref_ok == n and known / n >= 0.99 and in_span / max(1, known) >= 0.99,
            targets=[LOF_ID, FASTA_ID, GTF_ID], curated_variants=n, ref_equals_reference=ref_ok, genes_found_in_annotation=known, variant_inside_gene_span=in_span, gene_symbol_equals_annotation_name=same_name,
            clinvar_class_by_curation_verdict={v: dict(c) for v, c in by_verdict.items()},
            note='the curation covers 69 genes with incomplete penetrance; ClinVar classes are shown only as a relationship between two independent curations')

    def t_pharmvar():
        """PharmVar GRCh38 haplotype definitions against the reference sequence and the annotation, and RefSeqGene accessions against CPIC."""
        z = zipfile.ZipFile(os.path.join(D, 'pharmvar', 'pharmvar-6.2.29.zip'))
        genes, names = state['genes'], state['names']
        by_symbol = collections.defaultdict(list)
        for gid, nm in names.items():
            if nm:
                by_symbol[nm].append(gid)
        flank = 10000                                         # the RefSeqGene sequences extend about 5 kb upstream and 2 kb downstream of the gene
        checked = equal = inserts = flank_ok = 0
        in_span = in_flank = rows = unresolved = 0
        bad = []
        refseq = {}
        n_vcf = vcf_ok = vcf_rec = 0
        for n in z.namelist():
            p = n.split('/')
            if n.endswith('.haplotypes.tsv') and p[2] == 'RefSeqGene':
                for l in z.read(n).decode().split('\n')[2:]:
                    c = l.split('\t')
                    if len(c) == 9 and c[3] != 'REFERENCE':
                        refseq[p[1]] = c[3].split('.')[0]
                        break
            if n.endswith('.haplotypes.tsv') and p[2] == 'GRCh38':
                gene = p[1]
                cands = [genes[g] for g in by_symbol.get(gene, [])]
                if not cands:
                    unresolved += 1
                for l in z.read(n).decode().split('\n')[2:]:
                    c = l.split('\t')
                    if len(c) != 9 or c[3] == 'REFERENCE':
                        continue
                    chrom, start, stop, ref, typ = NC2CHR[c[3][:9]], int(c[4]), int(c[5]), c[6], c[8]
                    rows += 1
                    if typ in ('substitution', 'deletion'):
                        checked += 1
                        if fa.fetch(chrom, start - 1, stop).upper() == ref:
                            equal += 1
                        elif len(bad) < 5:
                            bad.append((gene, c[0], c[3], start, stop, ref))
                    else:
                        inserts += 1
                        flank_ok += stop == start + 1
                    g = next((x for x in cands if x[0] == chrom), None)
                    if g:
                        in_span += g[1] <= start and stop <= g[2]
                        in_flank += g[1] - flank <= start and stop <= g[2] + flank
            elif n.endswith('.vcf') and p[2] == 'GRCh38':
                n_vcf += 1
                ok = True
                for l in z.read(n).decode().split('\n'):
                    if not l or l[0] == '#':
                        continue
                    c = l.split('\t')
                    vcf_rec += 1
                    ok &= fa.fetch(c[0][3:] if c[0] != 'chrM' else 'MT', int(c[1]) - 1, int(c[1]) - 1 + len(c[3])).upper() == c[3]
                vcf_ok += ok
        src('pharmvar', 'PharmVar GRCh38 haplotype definitions: REF alleles equal the reference sequence and variants lie in their gene', rows > 1000 and checked + inserts == rows and equal == checked and flank_ok == inserts and n_vcf > 0 and vcf_ok == n_vcf and in_flank / max(1, rows) >= 0.99 and unresolved == 0,
            targets=[PV_ID, FASTA_ID, GTF_ID], table_rows=rows, substitutions_and_deletions_checked=checked, ref_equals_reference=equal, insertions=inserts, insertions_with_adjacent_flank_coordinates=flank_ok,
            vcf_files=n_vcf, vcf_files_with_all_ref_equal=vcf_ok, vcf_records=vcf_rec, rows_inside_gene_span=in_span, rows_within_10kb_of_gene=in_flank, genes_not_found_in_annotation=unresolved, mismatch_examples=bad,
            note='the table uses PharmVar-style coordinates (insertions sit between two flanking bases); the VCFs use VCF conventions; both are on GRCh38 chromosomes (NC_0000NN.11 and chrN)')
        cp = {x['symbol']: x for x in json.load(open(os.path.join(D, 'cpic', 'cpic_genes.json')))}
        both = sorted(set(refseq) & set(cp))
        same_ng = [g for g in both if (cp[g].get('genesequenceid') or '').split('.')[0] == refseq[g]]
        same_ens = [g for g in both if cp[g].get('ensemblid') in by_symbol.get(g, [])]
        src('pharmvar', 'PharmVar genes and RefSeqGene accessions agree with the CPIC gene table and the Ensembl annotation', len(both) >= 10 and len(same_ng) / len(both) >= 0.9 and len(same_ens) / len(both) >= 0.9,
            targets=[PV_ID, 'cpic_genes', GTF_ID], pharmvar_genes=len(refseq), also_in_cpic=len(both), same_refseqgene_accession=len(same_ng), same_ensembl_gene_id=len(same_ens),
            different_refseqgene=[(g, refseq[g], cp[g].get('genesequenceid')) for g in both if g not in same_ng], not_in_cpic=sorted(set(refseq) - set(cp)))

    guarded(src, 'ensembl_genomes', 'Remote range reads of the reference equal the mirrored FASTA', t_slice)
    guarded(src, 'clinvar', 'REF allele of every ClinVar record equals the GRCh38 reference sequence', t_clinvar)
    guarded(src, 'ensembl_genomes', 'Tests of the Ensembl release-116 annotation', t_annotation)
    if 'tx' in state:
        guarded(src, 'alphamissense', 'AlphaMissense against the reference sequence and annotation', t_alphamissense)
        if 'lof_rows' in state:
            guarded(src, 'gnomad_constraint', 'gnomAD LoF curation against the reference sequence and annotation', t_lof)
        if MAN.get(PV_ID, {}).get('ok'):
            guarded(src, 'pharmvar', 'PharmVar against the reference sequence and annotation', t_pharmvar)
