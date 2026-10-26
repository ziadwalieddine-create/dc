"""Mini source files in the real formats, built on the synthetic reference, for the source -> index -> matcher tests.
Everything here is synthetic. REF alleles are read from the synthetic reference so records are valid unless a test needs them not to be."""
import json
import os

import pysam

HERE = os.path.dirname(os.path.abspath(__file__))
REF = pysam.FastaFile(os.path.join(HERE, 'syn38.fa'))


def b(c, pos):
    return REF.fetch(c, pos - 1, pos).upper()


def alt(c, pos, skip=()):
    base = b(c, pos)
    comp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}[base]
    return next(x for x in 'ACGT' if x not in (base, comp) and x not in skip)


def mini_clinvar_vcf(path):
    hdr = ['##fileformat=VCFv4.1', '##fileDate=2026-01-01', '##source=ClinVar', '##reference=SYN38',
           '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO']
    L = []

    def rec(chrom, pos, vid, ref, a, **info):
        L.append('\t'.join([chrom, str(pos), str(vid), ref, a, '.', '.', ';'.join(f'{k}={v}' for k, v in info.items())]))
    rec('1', 402, 1001, b('1', 402), alt('1', 402), ALLELEID=2001, CLNDISDB='MedGen:C0000001,OMIM:100001', CLNDN='Test_syndrome', CLNREVSTAT='reviewed_by_expert_panel',
        CLNSIG='Pathogenic', CLNSIGSCV='SCV000000001|SCV000000002', CLNVC='single_nucleotide_variant', CLNVI='ClinGen:CA000001', GENEINFO='GENEA:1001', MC='SO:0001583|missense_variant',
        ORIGIN='1', RS='111')
    rec('1', 410, 1002, b('1', 410), alt('1', 410), ALLELEID=2002, CLNDISDB='MedGen:C0000001', CLNDN='Test_syndrome', CLNREVSTAT='criteria_provided,_multiple_submitters,_no_conflicts',
        CLNSIG='Benign/Likely_benign', CLNSIGSCV='SCV000000003|SCV000000004', CLNVC='single_nucleotide_variant', GENEINFO='GENEA:1001', ORIGIN='1')
    rec('1', 420, 1003, b('1', 420), alt('1', 420), ALLELEID=2003, CLNDISDB='MedGen:CN169374', CLNDN='not_specified', CLNREVSTAT='criteria_provided,_single_submitter',
        CLNSIG='Uncertain_significance', CLNSIGSCV='SCV000000005', CLNVC='single_nucleotide_variant', GENEINFO='GENEA:1001|GENEB:1002', ORIGIN='1')
    rec('1', 430, 1004, b('1', 430), alt('1', 430), ALLELEID=2004, CLNDISDB='MedGen:C0000002|MedGen:C0000003', CLNDN='Cond_A|Cond_B', CLNREVSTAT='criteria_provided,_conflicting_classifications',
        CLNSIG='Conflicting_classifications_of_pathogenicity', CLNSIGCONF='Pathogenic(1)|Likely_benign(1)', CLNSIGSCV='SCV000000006|SCV000000007', CLNVC='single_nucleotide_variant', ORIGIN='1')
    rec('1', 105, 1005, 'AA', 'A', ALLELEID=2005, CLNDISDB='MedGen:C0000004', CLNDN='Homopolymer_disorder', CLNREVSTAT='criteria_provided,_single_submitter', CLNSIG='Likely_pathogenic',
        CLNSIGSCV='SCV000000008', CLNVC='Deletion', ORIGIN='1')                       # same deletion as the next record after normalisation
    rec('1', 101, 1006, 'CA', 'C', ALLELEID=2006, CLNDISDB='MedGen:C0000004', CLNDN='Homopolymer_disorder', CLNREVSTAT='no_assertion_criteria_provided', CLNSIG='Pathogenic',
        CLNSIGSCV='SCV000000009', CLNVC='Deletion', ORIGIN='1')
    s = REF.fetch('1', 499, 501).upper()
    rec('1', 500, 1007, s[:2], ''.join('C' if x != 'C' else 'G' for x in s[:2]), ALLELEID=2007, CLNDISDB='MedGen:C0000005', CLNDN='Mnv_disorder', CLNREVSTAT='criteria_provided,_single_submitter',
        CLNSIG='Pathogenic', CLNSIGSCV='SCV000000010', CLNVC='Microsatellite', ORIGIN='1')
    wrong = 'A' if b('1', 440) != 'A' else 'C'
    rec('1', 440, 1008, wrong, 'T' if wrong != 'T' else 'G', ALLELEID=2008, CLNDISDB='MedGen:C0000006', CLNDN='Wrong_ref', CLNREVSTAT='criteria_provided,_single_submitter', CLNSIG='Pathogenic',
        CLNSIGSCV='SCV000000011', CLNVC='single_nucleotide_variant', ORIGIN='1')
    rec('1', 450, 1009, b('1', 450), '.', ALLELEID=2009, CLNDISDB='MedGen:C0000007', CLNDN='No_alt', CLNREVSTAT='no_assertion_criteria_provided', CLNSIG='not_provided', CLNVC='Variation', ORIGIN='1')
    rec('1', 460, 1010, b('1', 460), alt('1', 460), ALLELEID=2010, ONC='Oncogenic', ONCREVSTAT='criteria_provided,_single_submitter', ONCDN='Tumour_x', ONCDISDB='MedGen:C0000008',
        ONCSCV='SCV000000012', SCI='Tier_I_-_Strong', SCIREVSTAT='criteria_provided,_single_submitter', SCIDN='Tumour_x', SCIDISDB='MedGen:C0000008', SCISCV='SCV000000013', CLNVC='single_nucleotide_variant', ORIGIN='2')
    rec('chrUn_xx', 100, 1011, 'A', 'T', ALLELEID=2011, CLNDISDB='MedGen:C0000009', CLNDN='Unplaced', CLNREVSTAT='criteria_provided,_single_submitter', CLNSIG='Pathogenic', CLNSIGSCV='SCV000000014', CLNVC='single_nucleotide_variant', ORIGIN='1')
    rec('1', 470, 1012, b('1', 470), alt('1', 470), AF_ESP='0.00039', AF_EXAC='0.00038', AF_TGP='0.0014', ALLELEID=2012, CLNDISDB='MedGen:CN169374', CLNDN='not_specified',
        CLNREVSTAT='criteria_provided,_single_submitter', CLNSIG='Likely_benign', CLNSIGSCV='SCV000000015', CLNVC='single_nucleotide_variant', ORIGIN='1')
    rec('1', 480, 1013, b('1', 480), 'N' + alt('1', 480), ALLELEID=2013, CLNDISDB='MedGen:C0000010', CLNDN='N_allele', CLNREVSTAT='criteria_provided,_single_submitter', CLNSIG='Pathogenic', CLNSIGSCV='SCV000000016', CLNVC='Insertion', ORIGIN='1')
    with open(path, 'w') as f:
        f.write('\n'.join(hdr + L) + '\n')


def mini_clingen(out):
    """Mini ERepo export, gene-disease validity CSV, dosage gene and region lists, recurrent CNV BED and actionability TSVs in the real layouts."""
    os.makedirs(os.path.join(out, 'actionability'), exist_ok=True)
    hdr = ['Variation', 'ClinVar Variation Id', 'Allele Registry Id', 'HGVS Expressions', 'HGNC Gene Symbol', 'Disease', 'Mondo Id', 'Mode of Inheritance', 'Assertion',
           'Applied Evidence Codes (Met)', 'Applied Evidence Codes (Not Met)', 'Summary of interpretation', 'PubMed Articles', 'Expert Panel', 'Guideline', 'Approval Date',
           'Published Date', 'Retracted', 'Evidence Repo Link', 'Uuid']
    rows = [
        ['NM_X.1(GENEA):c.1A>G', '1001', 'CA000001', 'NM_X.1:c.1A>G', 'GENEA', 'Test syndrome', 'MONDO:0000001', 'Autosomal dominant inheritance', 'Pathogenic', 'PVS1,PS3_Moderate,PM2_Supporting', 'PS4,BA1', 'summary', '111,222', 'Test VCEP', 'ACMG spec v1', '2025-01-01', '2025-02-01', 'false', 'https://example.org/1', 'uuid-1'],
        ['NM_X.1(GENEA):c.2G>T', '1002', 'CA000002', '', 'GENEA', 'Test syndrome', 'MONDO:0000001', 'Autosomal dominant inheritance', 'Likely Benign', 'BP4,BS1_Supporting', 'PP3', '', '', 'Test VCEP', 'ACMG spec v1', '2025-01-02', '', 'false', '', 'uuid-2'],
        ['NM_X.1(GENEA):c.3C>A', '99999', 'CA999999', '', 'GENEA', 'Test syndrome', 'MONDO:0000001', 'Autosomal dominant inheritance', 'Uncertain Significance', 'PM2_Supporting', 'PM1', '', '', 'Test VCEP', 'ACMG spec v1', '2025-01-03', '', 'false', '', 'uuid-3'],
        ['NM_Y.1(GENEB):c.9del', '', '', 'NM_Y.1:c.9del', 'GENEB', 'Other', '', '', 'Pathogenic', 'PVS1', '', '', '', 'Test VCEP', 'ACMG spec v1', '2025-01-04', '', 'false', '', 'uuid-4'],
        ['NM_X.1(GENEA):c.1A>G', '', 'CA000001', '', 'GENEA', 'Test syndrome', 'MONDO:0000001', 'Autosomal dominant inheritance', 'Pathogenic', 'PS4', '', '', '', 'Second VCEP', 'ACMG spec v1', '2025-01-05', '', 'false', '', 'uuid-5'],
    ]
    with open(os.path.join(out, 'erepo_variant_classifications.tsv'), 'w') as f:
        f.write('\t'.join(hdr) + '\n' + '\n'.join('\t'.join(r) for r in rows) + '\n')
    q = lambda *x: ','.join('"%s"' % v for v in x)
    with open(os.path.join(out, 'clingen_gene_disease_validity.csv'), 'w') as f:
        f.write('\n'.join([q('CLINGEN GENE DISEASE VALIDITY CURATIONS', '', '', '', '', '', '', '', '', ''), q('FILE CREATED: 2026-01-01', *[''] * 9), q('WEBPAGE: x', *[''] * 9), q(*['+++'] * 10),
                           q('GENE SYMBOL', 'GENE ID (HGNC)', 'DISEASE LABEL', 'DISEASE ID (MONDO)', 'MOI', 'SOP', 'CLASSIFICATION', 'ONLINE REPORT', 'CLASSIFICATION DATE', 'GCEP'), q(*['+++'] * 10),
                           q('AARS', 'HGNC:20', 'CMT type 2N', 'MONDO:0013212', 'AD', 'SOP10', 'Definitive', 'https://x/CGGV:assertion_1', '2024-03-14T16:00:00.000Z', 'CMT GCEP'),
                           q('NOSUCHGENE', 'HGNC:99999999', 'Imaginary', 'MONDO:0000002', 'AR', 'SOP10', 'Limited', 'https://x/CGGV:assertion_2', '2024-03-15T16:00:00.000Z', 'Test GCEP'),
                           q('BRCA1', 'HGNC:1100', 'Disputed example', 'MONDO:0000003', 'AD', 'SOP10', 'Disputed', 'https://x/CGGV:assertion_3', '2024-03-16T16:00:00.000Z', 'Test GCEP')]) + '\n')
    gh = '#Gene Symbol\tGene ID\tcytoBand\tGenomic Location\tHaploinsufficiency Score\tHaploinsufficiency Description\tHaploinsufficiency PMID1\tHaploinsufficiency PMID2\tHaploinsufficiency PMID3\tHaploinsufficiency PMID4\tHaploinsufficiency PMID5\tHaploinsufficiency PMID6\tTriplosensitivity Score\tTriplosensitivity Description\tTriplosensitivity PMID1\tTriplosensitivity PMID2\tTriplosensitivity PMID3\tTriplosensitivity PMID4\tTriplosensitivity PMID5\tTriplosensitivity PMID6\tDate Last Evaluated\tHaploinsufficiency Disease ID\tTriplosensitivity Disease ID'
    with open(os.path.join(out, 'ClinGen_gene_curation_list_GRCh38.tsv'), 'w') as f:
        f.write('#ClinGen Gene Curation Results\n' + gh + '\n' + 'BRCA1\t672\t17q21.31\tchr17:43044295-43170245\t3\tSufficient evidence\t1\t2\t\t\t\t\t0\tNo evidence available\t\t\t\t\t\t\t2025-01-01\tMONDO:0011450\t\n'
                'NOSUCHGENE\t0\t1p1\tchr1:1-2\t30\tGene associated with autosomal recessive phenotype\t\t\t\t\t\t\t40\tDosage sensitivity unlikely\t\t\t\t\t\t\t2025-01-02\t\t\n')
    rh = '#ISCA ID\tISCA Region Name\tcytoBand\tGenomic Location\t' + gh.split('\t', 4)[4]
    with open(os.path.join(out, 'ClinGen_region_curation_list_GRCh38.tsv'), 'w') as f:
        f.write('#ClinGen Region Curation Results\n' + rh + '\n' + 'ISCA-1\t1q region\t1q1\tchr1:1000-2000\t3\tSufficient\t\t\t\t\t\t\t1\tLittle\t\t\t\t\t\t\t2025-01-01\tMONDO:0000001\t\n'
                'ISCA-2\tbroken region\t1q2\ttbd\t0\tNo evidence\t\t\t\t\t\t\t0\tNo evidence\t\t\t\t\t\t\t2025-01-01\t\t\n')
    with open(os.path.join(out, 'ClinGen_recurrent_CNV_V2.1-hg38.bed'), 'w') as f:
        f.write('track name="x"\n#chrom\tstart\tend\tname\tscore\tstrand\tthickStart\tthickEnd\titemRgb\nchr1\t999\t2000\tBP1\t0\t.\t999\t2000\t0,0,0\n')
    ah = '# docId\tiri\tcurationType\tlatestSearchDate\tlastUpdated\tlastAuthor\tcontext\tcontextIri\trelease\tgeneOrVariant\tgeneOmim\tdisease\tomim\tmondo\tsuggestedAssertion\tscorer\tpreliminaryAssertion\tconsensusAssertion\tstatus-assertion\tstatus-overall\tstatus-stg1'
    sh = '# docId\tiri\tcurationType\tlatestSearchDate\tlastUpdated\tlastAuthor\tcontext\tcontextIri\trelease\treleaseDate\tgeneOrVariant\tgeneOmim\tdisease\tomim\tstatus-overall\tstatus-stg1\tstatus-stg2\tstatus-scoring\toutcome\toutcomeScoringGroup\tintervention\tinterventionScoringGroup\tscorer\tprelim-severity\tprelim-likelihood\tprelim-natureOfIntervention\tprelim-effectiveness\tfinal-severity\tfinal-likelihood\tfinal-natureOfIntervention\tfinal-effectiveness\tfinal-overall'
    for ctxn, C in (('adult', 'Adult'), ('pediatric', 'Pediatric')):
        with open(os.path.join(out, 'actionability', f'{ctxn}_assertions_all.tsv'), 'w') as f:
            f.write(ah + '\n' + f'AC1\t/AC1\tGene-Condition\tx\tx\tx\t{C}\thttps://x\t1.0.0\tBRCA1\t113705\tBreast cancer\t114480\tMONDO_0007254\t\tA\t\tStrong Actionability\tComplete\tReleased\tComplete\n'
                    f'AC2\t/AC2\tVariant-Condition\tx\tx\tx\t{C}\thttps://x\t1.0.0\tNM_000059.4(BRCA2):c.1A>G\t600185\tHBOC\t604370\t\t\tA\t\tAssertion Pending\tComplete\tReleased\tComplete\n')
        with open(os.path.join(out, 'actionability', f'{ctxn}_scoring.tsv'), 'w') as f:
            f.write(sh + '\n' + f'AC1\t/AC1\tGene-Condition\tx\tx\tx\t{C}\thttps://x\t1.0.0\tx\tBRCA1\t113705\tBreast cancer\t114480\tReleased\tComplete\tComplete\tComplete\tcancer\tGroupA\tsurveillance\tGroupA\tA\t2\t3C\t3\t2B\t2\t3C\t3\t2B\t2N\n')


GWAS_COLS = ['DATE ADDED TO CATALOG', 'PUBMEDID', 'FIRST AUTHOR', 'DATE', 'JOURNAL', 'LINK', 'STUDY', 'DISEASE/TRAIT', 'INITIAL SAMPLE SIZE', 'REPLICATION SAMPLE SIZE', 'REGION', 'CHR_ID', 'CHR_POS',
             'REPORTED GENE(S)', 'MAPPED_GENE', 'UPSTREAM_GENE_ID', 'DOWNSTREAM_GENE_ID', 'SNP_GENE_IDS', 'UPSTREAM_GENE_DISTANCE', 'DOWNSTREAM_GENE_DISTANCE', 'STRONGEST SNP-RISK ALLELE', 'SNPS',
             'MERGED', 'SNP_ID_CURRENT', 'CONTEXT', 'INTERGENIC', 'RISK ALLELE FREQUENCY', 'P-VALUE', 'PVALUE_MLOG', 'P-VALUE (TEXT)', 'OR or BETA', '95% CI (TEXT)', 'PLATFORM [SNPS PASSING QC]',
             'CNV', 'MAPPED_TRAIT', 'MAPPED_TRAIT_URI', 'STUDY ACCESSION', 'GENOTYPING TECHNOLOGY']


def gwas_row(**kw):
    base = {c: '' for c in GWAS_COLS}
    base.update({'PUBMEDID': '1', 'DISEASE/TRAIT': 'Test trait', 'MAPPED_TRAIT': 'test trait', 'MAPPED_TRAIT_URI': 'http://www.ebi.ac.uk/efo/EFO_0000001', 'STUDY ACCESSION': 'GCST000001',
                 'MERGED': '0', 'P-VALUE': '1E-10', 'RISK ALLELE FREQUENCY': 'NR', 'SNP_GENE_IDS': 'ENSG0001', 'CONTEXT': 'intron_variant'})
    base.update(kw)
    return '\t'.join(base[c] for c in GWAS_COLS)


def mini_gwas(out):
    """Associations zip, studies and ancestries in the real layouts. Returns the cases so tests can refer to them by name."""
    import zipfile
    C = lambda pos, n=1: REF.fetch('1', pos - 1, pos - 1 + n).upper()
    comp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    # pick non-palindromic sites for the plain cases and an A/T or C/G site for the palindromic ones
    plain = [p for p in range(402, 700) if C(p) in 'ACGT']
    cases = {}
    rows = []

    def add(name, **kw):
        cases[name] = kw
        rows.append(gwas_row(**kw))
    p1 = 402
    r1 = C(p1)
    a1 = next(x for x in 'ACGT' if x not in (r1, comp[r1]))
    cases['plus_alt'] = dict(pos=p1, ref=r1, alt=a1)
    rows.append(gwas_row(**{'CHR_ID': '1', 'CHR_POS': str(p1), 'STRONGEST SNP-RISK ALLELE': f'rs101-{a1}', 'SNPS': 'rs101', 'RISK ALLELE FREQUENCY': '0.3', 'P-VALUE': '3E-9', 'OR or BETA': '1.25', '95% CI (TEXT)': '[1.1-1.4]'}))
    p2 = 410
    r2 = C(p2)
    a2 = next(x for x in 'ACGT' if x not in (r2, comp[r2]))
    cases['ref_effect'] = dict(pos=p2, ref=r2, alt=a2)
    rows.append(gwas_row(**{'CHR_ID': '1', 'CHR_POS': str(p2), 'STRONGEST SNP-RISK ALLELE': f'rs102-{r2}', 'SNPS': 'rs102', 'RISK ALLELE FREQUENCY': '0.7', 'P-VALUE': '2E-6', 'OR or BETA': '0.04', '95% CI (TEXT)': '[0.02-0.06] unit decrease'}))
    p3 = 420
    r3 = C(p3)
    a3 = next(x for x in 'ACGT' if x not in (r3, comp[r3]))
    cases['minus_strand'] = dict(pos=p3, ref=r3, alt=a3, reported=comp[a3])
    rows.append(gwas_row(**{'CHR_ID': '1', 'CHR_POS': str(p3), 'STRONGEST SNP-RISK ALLELE': f'rs103-{comp[a3]}', 'SNPS': 'rs103', 'OR or BETA': '1.5', '95% CI (TEXT)': '[1.2-1.9]'}))
    pal = [p for p in plain if C(p) in 'AT'][:3] if [p for p in plain if C(p) in 'AT'] else plain[:3]
    pa, pb = pal[0], pal[1]
    ra, rb = C(pa), C(pb)
    cases['palindromic_resolved'] = dict(pos=pa, ref=ra, alt=comp[ra], raf='0.8')
    rows.append(gwas_row(**{'CHR_ID': '1', 'CHR_POS': str(pa), 'STRONGEST SNP-RISK ALLELE': f'rs104-{ra}', 'SNPS': 'rs104', 'RISK ALLELE FREQUENCY': '0.8', 'OR or BETA': '1.1', '95% CI (TEXT)': '[1.05-1.15]'}))
    cases['palindromic_unresolved'] = dict(pos=pb, ref=rb, alt=comp[rb])
    rows.append(gwas_row(**{'CHR_ID': '1', 'CHR_POS': str(pb), 'STRONGEST SNP-RISK ALLELE': f'rs105-{rb}', 'SNPS': 'rs105', 'RISK ALLELE FREQUENCY': 'NR', 'OR or BETA': '1.1', '95% CI (TEXT)': '[1.05-1.15]'}))
    p6 = 430
    r6 = C(p6)
    a6 = next(x for x in 'ACGT' if x not in (r6, comp[r6]))
    cases['unknown_allele'] = dict(pos=p6, ref=r6, alt=a6)
    rows.append(gwas_row(**{'CHR_ID': '1', 'CHR_POS': str(p6), 'STRONGEST SNP-RISK ALLELE': 'rs106-?', 'SNPS': 'rs106', 'OR or BETA': '2.0', '95% CI (TEXT)': '[1.5-2.5]'}))
    rows.append(gwas_row(**{'STRONGEST SNP-RISK ALLELE': 'rs1-A; rs2-G', 'SNPS': 'rs1; rs2', 'CHR_ID': '1;1', 'CHR_POS': '402;410'}))
    rows.append(gwas_row(**{'STRONGEST SNP-RISK ALLELE': 'rs1-? x rs2-?', 'SNPS': 'rs1 x rs2', 'CHR_ID': '1 x 2', 'CHR_POS': '402 x 410'}))
    rows.append(gwas_row(**{'STRONGEST SNP-RISK ALLELE': 'HLA-A*31:01-?', 'SNPS': 'HLA-A*31:01'}))
    p9 = 440
    r9 = C(p9)
    a9 = next(x for x in 'ACGT' if x not in (r9, comp[r9]))
    cases['chr_verified'] = dict(pos=p9, ref=r9, alt=a9)
    rows.append(gwas_row(**{'STRONGEST SNP-RISK ALLELE': f'chr1:{p9}-{a9}', 'SNPS': f'chr1:{p9}', 'OR or BETA': '0.2', '95% CI (TEXT)': '[0.1-0.3] unit increase'}))
    rows.append(gwas_row(**{'STRONGEST SNP-RISK ALLELE': 'chr1:2900-A', 'SNPS': 'chr1:2900'}))
    rows.append(gwas_row(**{'STRONGEST SNP-RISK ALLELE': 'chr1:128781506:D-C', 'SNPS': 'chr1:128781506:D'}))
    rows.append(gwas_row(**{'STRONGEST SNP-RISK ALLELE': 'rs107-CT', 'SNPS': 'rs107', 'CHR_ID': '1', 'CHR_POS': '450'}))
    rows.append(gwas_row(**{'STRONGEST SNP-RISK ALLELE': 'rs108-A', 'SNPS': 'rs108', 'CHR_ID': '', 'CHR_POS': ''}))
    p14 = 460
    r14 = C(p14)
    a14 = next(x for x in 'ACGT' if x not in (r14, comp[r14]))
    cases['merged_rsid_p_zero'] = dict(pos=p14, ref=r14, alt=a14)
    rows.append(gwas_row(**{'CHR_ID': '1', 'CHR_POS': str(p14), 'STRONGEST SNP-RISK ALLELE': f'rs109-{a14}', 'SNPS': 'rs109', 'MERGED': '1', 'SNP_ID_CURRENT': '110', 'P-VALUE': '0', 'OR or BETA': '3.1', '95% CI (TEXT)': 'unit increase'}))
    p15 = 470
    r15 = C(p15)
    cases['no_alt_known'] = dict(pos=p15, ref=r15)
    rows.append(gwas_row(**{'CHR_ID': '1', 'CHR_POS': str(p15), 'STRONGEST SNP-RISK ALLELE': f'rs111-{next(x for x in "ACGT" if x != r15)}', 'SNPS': 'rs111'}))
    with zipfile.ZipFile(os.path.join(out, 'gwas-catalog-associations_ontology-annotated-full.zip'), 'w', zipfile.ZIP_DEFLATED) as z:
        z.writestr('gwas-catalog-download-associations-alt-full.tsv', '\t'.join(GWAS_COLS) + '\n' + '\n'.join(rows) + '\n')
    sh = ['DATE ADDED TO CATALOG', 'PUBMED ID', 'FIRST AUTHOR', 'DATE', 'JOURNAL', 'LINK', 'STUDY', 'DISEASE/TRAIT', 'INITIAL SAMPLE SIZE', 'REPLICATION SAMPLE SIZE', 'PLATFORM [SNPS PASSING QC]',
          'ASSOCIATION COUNT', 'MAPPED_TRAIT', 'MAPPED_TRAIT_URI', 'STUDY ACCESSION', 'GENOTYPING TECHNOLOGY', 'SUBMISSION DATE', 'STATISTICAL MODEL', 'BACKGROUND TRAIT', 'MAPPED BACKGROUND TRAIT',
          'MAPPED BACKGROUND TRAIT URI', 'COHORT', 'FULL SUMMARY STATISTICS', 'SUMMARY STATS LOCATION', 'GXE']
    with open(os.path.join(out, 'gwas-catalog-download-studies-v1.0.3.1.tsv'), 'w') as f:
        d = {c: '' for c in sh}
        d.update({'PUBMED ID': '1', 'DISEASE/TRAIT': 'Test trait', 'STUDY ACCESSION': 'GCST000001', 'COHORT': 'TESTCOHORT', 'FULL SUMMARY STATISTICS': 'yes', 'MAPPED_TRAIT_URI': 'http://www.ebi.ac.uk/efo/EFO_0000001'})
        f.write('\t'.join(sh) + '\n' + '\t'.join(d[c] for c in sh) + '\n')
    ah = ['STUDY ACCESSION', 'PUBMED ID', 'FIRST AUTHOR', 'DATE', 'INITIAL SAMPLE DESCRIPTION', 'REPLICATION SAMPLE DESCRIPTION', 'STAGE', 'NUMBER OF INDIVIDUALS', 'BROAD ANCESTRAL CATEGORY', 'COUNTRY OF ORIGIN',
          'COUNTRY OF RECRUITMENT', 'ADDITIONAL ANCESTRY DESCRIPTION', 'ANCESTRY DESCRIPTOR', 'FOUNDER/GENETICALLY ISOLATED POPULATION', 'NUMBER OF CASES', 'NUMBER OF CONTROLS', 'SAMPLE DESCRIPTION']
    with open(os.path.join(out, 'gwas-catalog-download-ancestries-v1.0.3.1.tsv'), 'w') as f:
        def arow(**kw):
            d = {c: '' for c in ah}
            d.update(kw)
            return '\t'.join(d[c] for c in ah)
        f.write('\t'.join(ah) + '\n' + arow(**{'STUDY ACCESSION': 'GCST000001', 'STAGE': 'initial', 'NUMBER OF INDIVIDUALS': '1000', 'NUMBER OF CASES': '400', 'BROAD ANCESTRAL CATEGORY': 'European'}) + '\n'
                + arow(**{'STUDY ACCESSION': 'GCST000001', 'STAGE': 'initial', 'NUMBER OF INDIVIDUALS': '500', 'BROAD ANCESTRAL CATEGORY': 'East Asian'}) + '\n'
                + arow(**{'STUDY ACCESSION': 'GCST000001', 'STAGE': 'replication', 'NUMBER OF INDIVIDUALS': '300', 'BROAD ANCESTRAL CATEGORY': 'African'}) + '\n')
    return cases


IUPAC2 = {frozenset('AG'): 'R', frozenset('CT'): 'Y', frozenset('AC'): 'M', frozenset('GT'): 'K', frozenset('CG'): 'S', frozenset('AT'): 'W'}


def mini_cpic(out):
    """Mini CPIC API tables in the real layout. TESTG: *1 reference, *2 defined by 1:705 and tolerating 1:710 (IUPAC code), *3 defined by 1:710.
    TESTH: one evaluable allele pair plus alleles the index cannot evaluate (repeat count, IUPAC code that excludes the reference, no location rows).
    NOCHR has no chromosome; MTG sits on the mitochondrial contig. Returns {name: value} of the bases the tests need."""
    os.makedirs(os.path.join(out, 'api'), exist_ok=True)
    g = lambda n: os.path.join(out, 'api', n + '.json')
    locs, defs, vals, alleles = [], [], [], []
    n = [100]

    def loc(chrom, text):
        n[0] += 1
        locs.append(dict(id=n[0], version=1, name=text, chromosomelocation=text, genelocation='', proteinlocation='', genesymbol='', dbsnpid=None, position=0))
        return n[0]

    def allele_def(gene, name, ref=False, structural=False, rows=(), func=None, act=None):
        n[0] += 1
        did = n[0]
        defs.append(dict(id=did, version=1, genesymbol=gene, name=name, pharmvarid=None, matchesreferencesequence=ref, structuralvariation=structural))
        for lid, a in rows:
            vals.append(dict(alleledefinitionid=did, locationid=lid, variantallele=a, version=1))
        n[0] += 1
        alleles.append(dict(id=n[0], version=1, genesymbol=gene, name=name, functionalstatus=func, clinicalfunctionalstatus=None, clinicalfunctionalsubstrate=None, activityvalue=act, definitionid=did,
                            citations=[], strength='Strong' if func else None, functioncomments=None, findings=None, frequency=None, inferredfrequency=False))
        return did

    b705, b710 = b('1', 705), b('1', 710)
    a705, a710 = alt('1', 705), alt('1', 710)
    l705 = loc('1', f'g.705{b705}>{a705}')
    l710 = loc('1', f'g.710{b710}>{a710}')
    code710 = IUPAC2[frozenset((b710, a710))]
    allele_def('TESTG', '*1', ref=True, rows=[(l705, b705), (l710, b710)], func='Normal function', act='1.0')
    allele_def('TESTG', '*2', rows=[(l705, a705), (l710, code710)], func='Decreased function', act='0.5')
    allele_def('TESTG', '*3', rows=[(l705, b705), (l710, a710)], func='No function', act='0.0')
    # TESTH: unevaluable alleles
    b402, a402 = b('1', 402), alt('1', 402)
    lh402 = loc('1', f'g.402{b402}>{a402}')
    lrep = loc('1', 'g.200TCAG[5]')
    lbad = loc('1', f'g.404{b("1", 404)}>{alt("1", 404)}')
    comp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    excl = next(k for k, v in IUPAC2.items() if b('1', 404) not in k)     # a code that does not contain the reference base
    allele_def('TESTH', '*1', ref=True, rows=[(lh402, b402), (lbad, b('1', 404))], func='Normal function')
    allele_def('TESTH', '*2', rows=[(lh402, a402), (lbad, b('1', 404))], func='Decreased function')
    allele_def('TESTH', '*3', rows=[(lh402, b402), (lrep, 'TCAG(7)')], func='Uncertain function')
    allele_def('TESTH', '*4', rows=[], func='Unknown function')
    allele_def('TESTH', '*5', rows=[(lh402, b402), (lbad, IUPAC2[excl])], func=None)
    allele_def('TESTH', '*6', rows=[(loc('1', 'g.406{0}>{1}'.format('C' if b('1', 406) != 'C' else 'G', 'T')), 'T')], func=None)      # reference base conflicts with the FASTA
    # genes without a chromosome or on the mitochondrial contig
    lno = loc('1', f'g.100{b("1", 100)}>{alt("1", 100)}')
    allele_def('NOCHR', '*1', ref=True, rows=[(lno, b('1', 100))])
    allele_def('NOCHR', '*2', rows=[(lno, alt('1', 100))])
    bm, am = b('MT', 50), alt('MT', 50)
    lm = loc('MT', f'm.50{bm}>{am}')
    allele_def('MTG', '*1', ref=True, rows=[(lm, bm)])
    allele_def('MTG', '*2', rows=[(lm, am)], func='No function')
    for name, rows in (('sequence_location', locs), ('allele_definition', defs), ('allele_location_value', vals), ('allele', alleles)):
        with open(g(name), 'w') as f:
            json.dump(rows, f)
    with open(os.path.join(out, 'cpic_genes.json'), 'w') as f:
        json.dump([dict(symbol='TESTG', chr='chr1'), dict(symbol='TESTH', chr='chr1'), dict(symbol='NOCHR', chr=None), dict(symbol='MTG', chr='chrM')], f)
    with open(os.path.join(out, 'cpic_guidelines.json'), 'w') as f:
        json.dump([dict(id=1, name='TESTG, TESTH and testdrug', genes=['TESTG', 'TESTH'], url='', clinpgxid='PA1', notesonusage=None, version=1)], f)
    rec = lambda i, drug, cls, ph, extra=None: dict(id=i, guidelineid=1, drugid=drug, implications={k: 'impl' for k in ph}, drugrecommendation='do x', classification=cls, phenotypes=ph, activityscore={},
                                                    allelestatus={}, lookupkey=dict(ph), population='general', comments='n/a', version=1, dosinginformation=False, alternatedrugavailable=False,
                                                    otherprescribingguidance=False)
    with open(g('recommendation'), 'w') as f:
        json.dump([rec(1, 'RxNorm:1', 'Strong', {'TESTG': 'Poor Metabolizer'}), rec(2, 'RxNorm:1', 'Optional', {'TESTG': 'Normal Metabolizer', 'TESTH': 'Poor Metabolizer'}),
                   rec(3, 'RxNorm:2', 'Moderate', {'UNKNOWNGENE': 'Poor Metabolizer'})], f)
    return dict(b705=b705, a705=a705, b710=b710, a710=a710, b402=b402, a402=a402)


def mini_assembly_report(path):
    hdr = '# Sequence-Name\tSequence-Role\tAssigned-Molecule\tAssigned-Molecule-Location/Type\tGenBank-Accn\tRelationship\tRefSeq-Accn\tAssembly-Unit\tSequence-Length\tUCSC-style-name'
    rows = [('1', 'SYN001.1', 'NC_SYN001.1', 3000, 'chr1'), ('MT', 'SYN002.1', 'NC_SYN002.1', 400, 'chrM')]
    with open(path, 'w') as f:
        f.write('# Assembly name:  SYN38\n' + hdr + '\n' + '\n'.join(f'{a}\tassembled-molecule\t{a}\tChromosome\t{g}\t=\t{r}\tPrimary Assembly\t{L}\t{u}' for a, g, r, L, u in rows) + '\n')


def mini_pharmvar(path):
    """A zip in the PharmVar layout (<release>/<gene>/GRCh38/<gene>.<refseq>.haplotypes.tsv). TESTG: *1 reference; *2 = 1:705 substitution; *2.001 = *2 + 1:712 substitution;
    *3 = 1:710 substitution; *4 = deletion of the two bases at 1:706-707; *5 = insertion between 1:712 and 713; *6 has a deleted sequence that differs from the reference."""
    import zipfile
    hdr = 'Haplotype Name\tGene\trsID\tReferenceSequence\tVariant Start\tVariant Stop\tReference Allele\tVariant Allele\tType'
    ns = 'NC_SYN001.1'
    sub = lambda hap, pos, rs: f'{hap}\tTESTG\t{rs}\t{ns}\t{pos}\t{pos}\t{b("1", pos)}\t{alt("1", pos)}\tsubstitution'
    rows = ['TESTG*1\tTESTG\t\tREFERENCE\t.\t\t\t\t', 'TESTG*1.001\tTESTG\t\tREFERENCE\t.\t\t\t\t',
            sub('TESTG*2', 705, 'rs705'), sub('TESTG*2.001', 705, 'rs705'), sub('TESTG*2.001', 712, 'rs712'), sub('TESTG*3', 710, 'rs710'),
            f'TESTG*4\tTESTG\t\t{ns}\t706\t707\t{REF.fetch("1", 705, 707).upper()}\t-\tdeletion',
            f'TESTG*5\tTESTG\t\t{ns}\t712\t713\t-\tTT\tinsertion',
            f'TESTG*6\tTESTG\t\t{ns}\t706\t707\tAA\t-\tdeletion',
            f'TESTG*7\tTESTG\t\tNC_UNKNOWN.1\t705\t705\tA\tC\tsubstitution']
    with zipfile.ZipFile(path, 'w') as z:
        z.writestr('pharmvar-test/TESTG/GRCh38/TESTG.NC_SYN001.1.haplotypes.tsv', '#version=pharmvar-test\n' + hdr + '\n' + '\n'.join(rows) + '\n')
        z.writestr('pharmvar-test/TESTG/GRCh37/TESTG.NC_SYN001.1.haplotypes.tsv', '#version=pharmvar-test\n' + hdr + '\n' + sub('TESTG*2', 705, 'rs705') + '\n')


def write_all(out):
    """Write every mini source file into `out` (a temporary directory: synthetic VCF-like files are never committed)."""
    os.makedirs(out, exist_ok=True)
    mini_clinvar_vcf(os.path.join(out, 'mini_clinvar.vcf'))
    mini_clingen(out)
    write_all.gwas_cases = mini_gwas(out)
    os.makedirs(os.path.join(out, 'cpic'), exist_ok=True)
    write_all.cpic = mini_cpic(os.path.join(out, 'cpic'))
    mini_assembly_report(os.path.join(out, 'assembly_report.txt'))
    mini_pharmvar(os.path.join(out, 'pharmvar-test.zip'))
    return out


if __name__ == '__main__':
    import sys
    print('written', write_all(sys.argv[1] if len(sys.argv) > 1 else 'sources'))
