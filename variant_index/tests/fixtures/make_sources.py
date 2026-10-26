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
    rec('1', 490, 1014, b('1', 490), alt('1', 490), ALLELEID=2014, CLNVC='single_nucleotide_variant', GENEINFO='GENEA:1001', RS='222|333', ORIGIN='1')        # no classification of any kind
    rec('1', 495, 1015, b('1', 495), alt('1', 495), ALLELEID=2015, CLNVC='single_nucleotide_variant', ORIGIN='1')                                         # no classification, no RS
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
    # TESTK: named (non-star) alleles, one of them with a comma in its name, as in the G6PD names
    b408, a408 = b('1', 408), alt('1', 408)
    lk = loc('1', f'g.408{b408}>{a408}')
    allele_def('TESTK', 'Ref (reference)', ref=True, rows=[(lk, b408)], func='Normal function')
    allele_def('TESTK', 'Odd, name', rows=[(lk, a408)], func='Decreased function')
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
        json.dump([dict(symbol='TESTG', chr='chr1'), dict(symbol='TESTH', chr='chr1'), dict(symbol='TESTK', chr='chr1'), dict(symbol='NOCHR', chr=None), dict(symbol='MTG', chr='chrM')], f)
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


def _tsv(header, rows):
    return '\t'.join(header) + '\n' + '\n'.join('\t'.join(str(r.get(h, '')) for h in header) for r in rows) + '\n'


SUMMARY_HDR = ['Summary Annotation ID', 'Variant/Haplotypes', 'Gene', 'Level of Evidence', 'Level Override', 'Level Modifiers', 'Score', 'Phenotype Category', 'PMID Count', 'Evidence Count',
               'Drug(s)', 'Phenotype(s)', 'Latest History Date (YYYY-MM-DD)', 'URL', 'Specialty Population']
VAR_DRUG_HDR = ['Variant Annotation ID', 'Variant/Haplotypes', 'Gene', 'Drug(s)', 'PMID', 'Phenotype Category', 'Significance', 'Notes', 'Sentence', 'Alleles', 'Specialty Population', 'Metabolizer types',
                'isPlural', 'Is/Is Not associated', 'Direction of effect', 'PD/PK terms', 'Multiple drugs And/or', 'Population types', 'Population Phenotypes or diseases',
                'Multiple phenotypes or diseases And/or', 'Comparison Allele(s) or Genotype(s)', 'Comparison Metabolizer types']
VAR_PHENO_HDR = ['Variant Annotation ID', 'Variant/Haplotypes', 'Gene', 'Drug(s)', 'PMID', 'Phenotype Category', 'Significance', 'Notes', 'Sentence', 'Alleles', 'Specialty Population', 'Metabolizer types',
                 'isPlural', 'Is/Is Not associated', 'Direction of effect', 'Side effect/efficacy/other', 'Phenotype', 'Multiple phenotypes And/or', 'When treated with/exposed to/when assayed with',
                 'Multiple drugs And/or', 'Population types', 'Population Phenotypes or diseases', 'Multiple phenotypes or diseases And/or', 'Comparison Allele(s) or Genotype(s)', 'Comparison Metabolizer types']
VAR_FA_HDR = ['Variant Annotation ID', 'Variant/Haplotypes', 'Gene', 'Drug(s)', 'PMID', 'Phenotype Category', 'Significance', 'Notes', 'Sentence', 'Alleles', 'Specialty Population', 'Assay type',
              'Metabolizer types', 'isPlural', 'Is/Is Not associated', 'Direction of effect', 'Functional terms', 'Gene/gene product', 'When treated with/exposed to/when assayed with',
              'Multiple drugs And/or', 'Cell type', 'Comparison Allele(s) or Genotype(s)', 'Comparison Metabolizer types']
STUDY_HDR = ['Study Parameters ID', 'Variant Annotation ID', 'Study Type', 'Study Cases', 'Study Controls', 'Characteristics', 'Characteristics Type', 'Frequency In Cases', 'Allele Of Frequency In Cases',
             'Frequency In Controls', 'Allele Of Frequency In Controls', 'P Value', 'Ratio Stat Type', 'Ratio Stat', 'Confidence Interval Start', 'Confidence Interval Stop', 'Biogeographical Groups']
LABEL_HDR = ['PharmGKB ID', 'Name', 'Source', 'Biomarker Flag', 'Testing Level', 'Has Prescribing Info', 'Has Dosing Info', 'Has Alternate Drug', 'Has Other Prescribing Guidance', 'Cancer Genome',
             'Prescribing', 'Chemicals', 'Genes', 'Variants/Haplotypes', 'Latest History Date (YYYY-MM-DD)']


def mini_clinpgx(out):
    """Mini ClinPGx archives in the real layout (variants, summary + legacy clinical annotations, variant annotations with study parameters, clinical variants,
    guideline annotations, drug labels) on the synthetic reference. Returns the keys the tests refer to."""
    import zipfile
    os.makedirs(out, exist_ok=True)
    ns = 'NC_SYN001.1'
    b705, a705 = b('1', 705), alt('1', 705)
    b710 = b('1', 710)
    syn = lambda *items: ', '.join(items)
    variants = [
        dict(**{'Variant ID': 'PA1', 'Variant Name': 'rs705', 'Gene IDs': 'PAG1', 'Gene Symbols': 'TESTG', 'Location': f'{ns}:705',
                'Synonyms': syn('705', f'{ns}:g.705{b705}>{a705}', f'{ns}:g.705=', 'PA1', 'rs705', 'rs9705')}),
        {'Variant ID': 'PA2', 'Variant Name': 'rs710', 'Gene IDs': 'PAG1', 'Gene Symbols': 'TESTG', 'Location': f'{ns}:710',
         'Synonyms': syn(f'{ns}:g.710{b710}>C', f'{ns}:g.710{b710}>T', f'{ns}:g.710=', 'rs710')},
        {'Variant ID': 'PA3', 'Variant Name': 'rs900001', 'Gene Symbols': '', 'Location': '', 'Synonyms': 'rs900001'},
        {'Variant ID': 'PA4', 'Variant Name': 'rs900002', 'Gene Symbols': '', 'Location': f'{ns}:402', 'Synonyms': syn(f'{ns}:g.402=', 'rs900002')},
        {'Variant ID': 'PA5', 'Variant Name': 'rs900003', 'Gene Symbols': '', 'Location': 'NC_999999.1:5', 'Synonyms': syn('NC_999999.1:g.5A>C', 'rs900003')},
        {'Variant ID': 'PA6', 'Variant Name': 'rs900004', 'Gene Symbols': 'TESTG', 'Location': f'{ns}:706', 'Synonyms': syn(f'{ns}:g.706_707del', 'rs900004')},
        {'Variant ID': 'PA7', 'Variant Name': 'rs900005', 'Gene Symbols': '', 'Location': f'{ns}:200', 'Synonyms': syn(f'{ns}:g.200_203TCAG[5]', 'rs900005')},
    ]
    vh = ['Variant ID', 'Variant Name', 'Gene IDs', 'Gene Symbols', 'Location', 'Variant Annotation count', 'Clinical Annotation count', 'Level 1/2 Clinical Annotation count',
          'Guideline Annotation count', 'Label Annotation count', 'Synonyms']

    def zip_(name, created, files):
        with zipfile.ZipFile(os.path.join(out, name + '.zip'), 'w') as z:
            z.writestr(f'CREATED_{created}.txt', f'Created on {created}\n')
            for fn, content in files.items():
                z.writestr(fn, content)

    zip_('variants', '2026-10-05', {'variants.tsv': _tsv(vh, variants)})
    S = lambda i, v, g, lvl, score, **kw: dict({'Summary Annotation ID': i, 'Variant/Haplotypes': v, 'Gene': g, 'Level of Evidence': lvl, 'Score': score, 'Phenotype Category': 'Efficacy',
                                                'PMID Count': 1, 'Evidence Count': 1, 'Drug(s)': 'testdrug', 'Phenotype(s)': 'Test disease', 'Latest History Date (YYYY-MM-DD)': '2026-01-01',
                                                'URL': f'https://www.clinpgx.org/summaryAnnotation/{i}'}, **kw)
    summary = [S(1001, 'rs705', 'TESTG', '1A', 12.5), S(1002, 'rs710', 'TESTG', '2B', 3), S(1003, 'rs710', 'TESTG', '3', 1.5), S(1004, 'rs900001', 'TESTG', '3', 1), S(1005, 'rs900002', 'TESTG', '3', 1),
               S(1006, 'TESTG*1, TESTG*2', 'TESTG', '1B', 6), S(1007, 'TESTG poor metabolizer', 'TESTG', '3', 1), S(1008, 'TESTK Ref (reference), TESTK Odd, name', 'TESTK', '2A', 4),
               S(1009, 'TESTG*9', 'TESTG', '3', 1), S(1010, 'rs705, rs710', 'TESTG;TESTH', '3', 2), S(1011, 'rs705', 'TESTG', '4', 0.5), S(1012, 'c.123A>G', 'TESTG', '3', 1),
               S(1013, 'rs900003', 'TESTG', '3', 1), S(1014, 'rs900004', 'TESTG', '3', 1), S(1015, 'rs900005', 'TESTG', '3', 1), S(1016, 'TESTG HTTLPR long form (L allele)', 'TESTG', '3', 1),
               S(1017, 'rs710', 'TESTG', '3', 1), S(1018, 'HLA-B*13:01', 'HLA-B', '3', 1)]
    ah = ['Summary Annotation ID', 'Genotype/Allele', 'Annotation Text', 'Allele Function']
    alleles = [dict(zip(ah, (1001, 'GG', 'Patients with GG may respond.', ''))), dict(zip(ah, (1001, 'AG', 'Patients with AG may respond less.', ''))), dict(zip(ah, (1001, 'AA', 'AA text', ''))),
               dict(zip(ah, (1002, 'CC', 'cc text', ''))), dict(zip(ah, (1002, 'AC', 'ac text', ''))), dict(zip(ah, (1003, 'T', 't text', ''))), dict(zip(ah, (1006, '*1/*2', 'text', 'Normal function'))),
               dict(zip(ah, (1010, 'AA', 'two sites', ''))), dict(zip(ah, (1017, 'CC', 'cc', ''))), dict(zip(ah, (1017, 'TT', 'tt', '')))]
    eh = ['Summary Annotation ID', 'Evidence ID', 'Evidence Type', 'Evidence URL', 'PMID', 'Summary', 'Score']
    evid = [dict(zip(eh, (1001, 2001, 'Variant Drug Annotation', 'https://x/2001', 111, 'summary text', 5.0)))]
    zip_('summaryAnnotations', '2026-10-05', {'summary_annotations.tsv': _tsv(SUMMARY_HDR, summary), 'summary_ann_alleles.tsv': _tsv(ah, alleles), 'summary_ann_evidence.tsv': _tsv(eh, evid),
                                               'summary_ann_history.tsv': 'Summary Annotation ID\tDate (YYYY-MM-DD)\tType\tComment\n'})
    legacy = [r for r in summary if r['Summary Annotation ID'] != 1011] + [S(1099, 'rs705', 'TESTG', '3', 1)]
    legacy = [dict(r, **{'Level of Evidence': '2A'}) if r['Summary Annotation ID'] == 1002 else r for r in legacy]
    lh = ['Clinical Annotation ID'] + SUMMARY_HDR[1:]
    legacy = [{**{('Clinical Annotation ID' if k == 'Summary Annotation ID' else k): v for k, v in r.items()}} for r in legacy]
    zip_('clinicalAnnotations', '2025-07-05', {'clinical_annotations.tsv': _tsv(lh, legacy)})
    cvh = ['variant', 'gene', 'type', 'level of evidence', 'chemicals', 'phenotypes']
    zip_('clinicalVariants', '2026-10-05', {'clinicalVariants.tsv': _tsv(cvh, [dict(zip(cvh, ('rs705', 'TESTG', 'Efficacy', '1A', 'testdrug', 'x'))), dict(zip(cvh, ('rs999999', 'TESTG', 'Toxicity', '2A', 'd', 'y')))])})

    D = lambda i, v, g, pmid, sig, **kw: dict({'Variant Annotation ID': i, 'Variant/Haplotypes': v, 'Gene': g, 'Drug(s)': 'testdrug', 'PMID': pmid, 'Phenotype Category': 'Efficacy', 'Significance': sig,
                                               'Is/Is Not associated': 'Associated with' if sig != 'no' else 'Not associated with', 'Sentence': f'sentence {i}'}, **kw)
    drug = [D(2001, 'rs705', 'TESTG', 111, 'yes', **{'Alleles': 'AA', 'Comparison Allele(s) or Genotype(s)': 'GG', 'Direction of effect': 'increased'}),
            D(2002, 'rs705', 'TESTG', 111, 'no', **{'Alleles': 'AA', 'Comparison Allele(s) or Genotype(s)': 'GG'}),
            D(2003, 'TESTG*1, TESTG*2', 'TESTG', 222, 'not stated', **{'Alleles': '*2/*2', 'Comparison Allele(s) or Genotype(s)': '*1/*1'}),
            D(2004, 'rs710', 'TESTG', 444, 'yes', **{'Alleles': 'GG'}), D(2005, 'TESTG*1', 'TESTG', 555, 'yes', **{'Alleles': '*1/*1'})]
    pheno = [dict(D(2101, 'rs705', 'TESTG', 333, 'yes', **{'Alleles': 'AA'}), Phenotype='Side Effect:Test reaction')]
    fa = [dict(D(2201, 'TESTG*2', 'TESTG', 666, 'yes', **{'Alleles': '*2'}), **{'Assay type': 'in human liver microsomes', 'Cell type': 'hepatocyte', 'Functional terms': 'activity of'})]
    P = lambda i, va, **kw: dict({'Study Parameters ID': i, 'Variant Annotation ID': va, 'Study Type': 'case/control'}, **kw)
    params = [P(3001, 2001, **{'Study Cases': 100, 'Study Controls': 200, 'P Value': '= 0.012', 'Ratio Stat Type': 'OR', 'Ratio Stat': 2.5, 'Confidence Interval Start': 1.1,
                               'Confidence Interval Stop': 5.7, 'Biogeographical Groups': 'European'}),
              P(3002, 2001, **{'Study Cases': 50, 'Study Controls': 60, 'P Value': '< 0.001', 'Ratio Stat Type': 'HR', 'Ratio Stat': 1.8, 'Biogeographical Groups': 'East Asian'}),
              P(3003, 2002, **{'Study Cases': 30, 'Study Controls': 40, 'P Value': '= 0.4', 'Ratio Stat Type': 'Unknown', 'Ratio Stat': '', 'Biogeographical Groups': 'Unknown'}),
              P(3004, 2101, **{'Study Cases': 10, 'Study Controls': 12, 'P Value': '= 1E-5', 'Ratio Stat Type': 'OR', 'Ratio Stat': 0.4, 'Biogeographical Groups': 'African American/Afro-Caribbean'})]
    zip_('variantAnnotations', '2026-10-05', {'var_drug_ann.tsv': _tsv(VAR_DRUG_HDR, drug), 'var_pheno_ann.tsv': _tsv(VAR_PHENO_HDR, pheno), 'var_fa_ann.tsv': _tsv(VAR_FA_HDR, fa),
                                               'study_parameters.tsv': _tsv(STUDY_HDR, params)})
    g1 = {'citations': [], 'guideline': {'objCls': 'GuidelineAnnotation', 'id': 'PA1000', 'name': 'Annotation of CPIC Guideline for testdrug and TESTG', 'source': 'CPIC', 'recommendation': True,
                                          'dosingInformation': True, 'hasTestingInfo': True, 'alternateDrugAvailable': False, 'otherPrescribingGuidance': False, 'cancerGenome': False, 'pediatric': False,
                                          'relatedGenes': [{'symbol': 'TESTG'}], 'relatedChemicals': [{'name': 'testdrug'}],
                                          'relatedAlleles': [{'objCls': 'Variant', 'symbol': 'rs705'}, {'objCls': 'Haplotype', 'symbol': 'TESTG*2'}], 'version': 3}}
    g2 = {'citations': [], 'guideline': {'objCls': 'GuidelineAnnotation', 'id': 'PA1001', 'name': 'Annotation of DPWG Guideline for other and TESTH', 'source': 'DPWG', 'recommendation': False,
                                          'relatedGenes': [{'symbol': 'TESTH'}, {'symbol': 'UNKNOWNX'}], 'relatedChemicals': [{'name': 'other'}], 'relatedAlleles': [], 'version': 1}}
    zip_('guidelineAnnotations.json', '2026-10-05', {'PA1000.json': json.dumps(g1), 'PA1001.json': json.dumps(g2)})
    zip_('drugLabels', '2026-10-05', {'drugLabels.tsv': _tsv(LABEL_HDR, [dict(zip(LABEL_HDR, ('PA2000', 'Annotation of FDA Label for testdrug and TESTG', 'FDA', 'Yes', 'Actionable PGx', 'Prescribing Info', '', '', '', '',
                                                                                             'Prescribing', 'testdrug', 'TESTG; TESTH', '', '2026-01-01')))])})
    return dict(k705=f'1:705:{b705}:{a705}', k710c=f'1:710:{b710}:C', k710t=f'1:710:{b710}:T', kdel='1:705:GCC:G')


def mini_mavedb(out):
    """Mini MaveDB dump metadata (main.json) and a mini Ensembl VEP MaveDB plugin file (same columns, plus two sparse replicate columns)."""
    import gzip
    os.makedirs(out, exist_ok=True)
    kw = [{'keyword': {'key': 'Phenotypic Assay Model System', 'label': 'Immortalized human cells'}}, {'keyword': {'key': 'Phenotypic Assay Mechanism', 'label': 'Loss of function'}},
          {'keyword': {'key': 'Variant Library Creation Method', 'label': 'Endogenous locus library method'}}]
    ext = [{'identifier': {'dbName': 'Ensembl', 'identifier': 'ENSG00000000001'}}]

    def sset(n, **kw2):
        d = {'urn': f'urn:mavedb:0000000{n}-a-1', 'title': f'Set {n}', 'numVariants': 10, 'supersedingScoreSet': None, 'mappingState': 'complete', 'publishedDate': '2024-01-01',
             'primaryPublicationIdentifiers': [], 'targetGenes': [{'name': 'TESTG', 'category': 'protein_coding', 'mappedHgncName': 'TESTG', 'externalIdentifiers': ext}], 'scoreCalibrations': []}
        d.update(kw2)
        return d
    cal = [{'urn': 'urn:mavedb:calibration-1', 'title': 'Test calibration', 'researchUseOnly': False, 'investigatorProvided': True, 'thresholdSources': [{'identifier': '99999'}],
            'functionalClassifications': [
                {'label': 'PS3 Moderate (2)', 'functionalClassification': 'abnormal', 'range': [None, 0.5], 'inclusiveLowerBound': False, 'inclusiveUpperBound': False,
                 'acmgClassification': {'criterion': 'PS3', 'evidenceStrength': 'MODERATE', 'points': 2}},
                {'label': 'BS3 Supporting (-1)', 'functionalClassification': 'normal', 'range': [1.5, None], 'inclusiveLowerBound': True, 'inclusiveUpperBound': False,
                 'acmgClassification': {'criterion': 'BS3', 'evidenceStrength': 'SUPPORTING', 'points': -1}}]},
           {'urn': 'urn:mavedb:calibration-2', 'title': 'Research calibration', 'researchUseOnly': True, 'investigatorProvided': False, 'thresholdSources': [],
            'functionalClassifications': [{'label': 'abnormal below 1', 'functionalClassification': 'abnormal', 'range': [None, 1.0], 'inclusiveLowerBound': False, 'inclusiveUpperBound': True,
                                           'acmgClassification': None}]}]
    sets = [sset(1, scoreCalibrations=cal), sset(2, supersedingScoreSet='urn:mavedb:00000006-a-1'), sset(3), sset(4), sset(5, scoreCalibrations=cal)]
    main = {'title': 'MaveDB public data', 'asOf': '2026-06-24T18:13:15+00:00', 'experimentSets': [{'urn': 'urn:mavedb:00000001', 'experiments': [
        {'urn': 'urn:mavedb:00000001-a', 'keywords': kw, 'primaryPublicationIdentifiers': [{'identifier': '12345', 'dbName': 'PubMed'}], 'scoreSets': sets}]}]}
    with open(os.path.join(out, 'main.json'), 'w') as f:
        json.dump(main, f)
    hdr = ['#chr', 'start', 'end', 'ref', 'alt', 'hgvs', 'urn', 'publish_date', 'refseq', 'pubmed', 'doi', 'url', 'accession', 'nt', 'splice', 'pro', 'score', 'sd', 'se', 'df', 'score_rep1', 'se_rep1']
    wrong = 'A' if b('1', 440) != 'A' else 'C'

    def row(urn_n, accn, chrom, start, ref, alt, hgvs, nt, pro, score, sd='', se='', rep=''):
        urn = f'urn:mavedb:0000000{urn_n}-a-1'
        return '\t'.join([chrom, str(start), str(start + len(ref) - 1), ref, alt, hgvs, urn, '2024-01-01', '', '12345', '', '', f'{urn}#{accn}', nt, '', pro, str(score), str(sd), str(se), '', rep, ''])
    r402 = b('1', 402), alt('1', 402)
    s500 = REF.fetch('1', 499, 501).upper()
    rows = [row(1, 1, '1', 402, r402[0], r402[1], 'NC_000001.11:g.402' + r402[0] + '>' + r402[1], 'c.1' + r402[0] + '>' + r402[1], 'p.Val1Ile', 0.2, 0.1, 0.05, rep='0.3'),
            row(1, 2, '1', 410, b('1', 410), alt('1', 410), 'NC_000001.11:g.410', 'c.9A>C', 'p.Ala3Pro', 1.8),
            row(1, 3, '1', 420, b('1', 420), alt('1', 420), 'NC_000001.11:g.420', 'c.19A>C', 'p.Ala7Pro', 1.0),
            row(1, 4, '1', 440, wrong, 'T' if wrong != 'T' else 'G', 'NC_000001.11:g.440', 'c.39A>C', 'p.Ala13Pro', 0.9),
            row(1, 1, '1', 402, r402[0], r402[1], 'NC_000001.11:g.402', 'c.1' + r402[0] + '>' + r402[1], 'p.Val1Ile', 0.2),          # the same allele and accession again
            row(3, 7, '1', 500, s500[:2], 'CC', 'NP_000001.1:p.Gly170Ala', '', 'p.Gly170Ala', -0.4),           # protein-level assay, candidate nucleotide allele (MNV)
            row(3, 7, '1', 500, s500[:2], 'TT', 'NP_000001.1:p.Gly170Ala', '', 'p.Gly170Ala', -0.4),
            row(2, 1, '1', 402, r402[0], r402[1], 'NC_000001.11:g.402', 'c.1' + r402[0] + '>' + r402[1], 'p.Val1Ile', 5.0),                       # superseded score set
            '\t'.join(['11', '102114330', '102114330', 'ga4gh:VA.GZCxxxx', 'A', 'NC_000011.10:g.102114330G>A', 'urn:mavedb:00000004-a-1', '2024-01-01', '', '', '', '', 'urn:mavedb:00000004-a-1#9',
                       'c.[1G>A;83C>G]', '', 'p.[Asp1Asn;Thr28Ser]', '-0.24', '', '', '', '', '', '']),
            row(1, 5, '1', 430, b('1', 430), alt('1', 430), 'NC_000001.11:g.430', 'c.29A>C', 'p.Ala10Pro', 'NA'),
            row(9, 1, '1', 430, b('1', 430), alt('1', 430), 'NC_000001.11:g.430', 'c.29A>C', 'p.Ala10Pro', 1.0)]
    with gzip.open(os.path.join(out, 'MaveDB_variants.tsv.gz'), 'wt') as f:
        f.write('\t'.join(hdr) + '\n' + '\n'.join(rows) + '\n')
    # MaveDB's own genomic mapping: annotations.csv rows with post-mapped genomic HGVS, and the matching scores.csv.gz
    os.makedirs(os.path.join(out, 'csv'), exist_ok=True)
    os.makedirs(os.path.join(out, 'scores_genomic_mapped'), exist_ok=True)
    ann_hdr = ['accession', 'hgvs_nt', 'hgvs_splice', 'hgvs_pro', 'mavedb.post_mapped_hgvs_g', 'mavedb.post_mapped_hgvs_p', 'mavedb.post_mapped_hgvs_c', 'mavedb.post_mapped_hgvs_at_assay_level',
               'mavedb.post_mapped_vrs_digest', 'vep.vep_functional_consequence', 'gnomad.gnomad_af', 'clingen.clingen_allele_id']
    u = 'urn:mavedb:00000005-a-1'
    g = lambda pos, r_, a_: f'NC_SYN001.1:g.{pos}{r_}>{a_}'
    ann = [[f'{u}#1', 'c.1G>A', 'NA', 'p.Val1Ile', g(402, r402[0], r402[1]), 'NA', 'NA', 'NA', 'digest1', 'missense_variant', 'NA', 'CA1'],
           [f'{u}#2', 'c.2G>A', 'NA', 'NA', g(440, wrong, 'T' if wrong != 'T' else 'G'), 'NA', 'NA', 'NA', 'digest2', 'NA', 'NA', 'CA2'],
           [f'{u}#3', 'c.[1G>A;9A>C]', 'NA', 'NA', f'NC_SYN001.1:g.[402{r402[0]}>{r402[1]};410{b("1", 410)}>{alt("1", 410)}]', 'NA', 'NA', 'NA', 'digest3', 'NA', 'NA', 'CA3'],
           [f'{u}#4', 'NA', 'NA', 'p.Ala3Pro', 'NA', 'NP_000001.1:p.Ala3Pro', 'NA', 'NA', 'digest4', 'NA', 'NA', 'CA4'],
           [f'{u}#5', 'c.5A>C', 'NA', 'NA', g(410, b('1', 410), alt('1', 410)), 'NA', 'NA', 'NA', 'digest5', 'NA', 'NA', 'CA5'],
           [f'{u}#6', 'c.7del', 'NA', 'NA', 'NC_SYN001.1:g.706_707del', 'NA', 'NA', 'NA', 'digest6', 'NA', 'NA', 'CA6'],
           [f'{u}#7', 'c.8=', 'NA', 'NA', 'NC_SYN001.1:g.705=', 'NA', 'NA', 'NA', 'digest7', 'NA', 'NA', 'CA7'],
           [f'{u}#8', 'c.9A>C', 'NA', 'NA', g(420, b('1', 420), alt('1', 420)), 'NA', 'NA', 'NA', 'digest8', 'NA', 'NA', 'CA8']]
    import csv as _csv
    with open(os.path.join(out, 'csv', 'urn-mavedb-00000005-a-1.annotations.csv'), 'w', newline='') as f:
        w = _csv.writer(f)
        w.writerow(ann_hdr)
        w.writerows(ann)
    sc_rows = [[f'{u}#1', 'c.1G>A', 'NA', 'p.Val1Ile', 0.3, 0.1, 0.02], [f'{u}#2', 'c.2G>A', 'NA', 'NA', 0.9, 0.1, 0.02], [f'{u}#3', 'c.[1G>A;9A>C]', 'NA', 'NA', 0.2, 0.1, 0.02],
               [f'{u}#4', 'NA', 'NA', 'p.Ala3Pro', 0.4, 0.1, 0.02], [f'{u}#6', 'c.7del', 'NA', 'NA', 1.7, 0.1, 0.02], [f'{u}#7', 'c.8=', 'NA', 'NA', 1.0, 0.1, 0.02],
           [f'{u}#8', 'c.9A>C', 'NA', 'NA', 'NA', 0.1, 0.02]]
    with gzip.open(os.path.join(out, 'scores_genomic_mapped', 'urn-mavedb-00000005-a-1.scores.csv.gz'), 'wt', newline='') as f:
        w = _csv.writer(f)
        w.writerow(['accession', 'hgvs_nt', 'hgvs_splice', 'hgvs_pro', 'scores.score', 'scores.sd', 'scores.se'])
        w.writerows(sc_rows)
    return dict(k402=f'1:402:{r402[0]}:{r402[1]}', k410=f'1:410:{b("1", 410)}:{alt("1", 410)}', kmnv1=f'1:500:{s500[:2]}:CC', kmnv2=f'1:500:{s500[:2]}:TT')


def mini_opentargets(out):
    """Mini Open Targets 26.09 genetics tables (biosample, study, l2g_prediction, credible_set) with the columns the adapter reads."""
    import pyarrow as pa
    import pyarrow.parquet as pq
    for d in ('biosample', 'study', 'l2g_prediction', 'credible_set'):
        os.makedirs(os.path.join(out, d), exist_ok=True)
    pq.write_table(pa.table({'biosampleId': ['EFO_0005292', 'UBERON_0000948'], 'biosampleName': ['lymphoblastoid cell line', 'heart']}), os.path.join(out, 'biosample', 'part-0.parquet'))
    st = [dict(studyId='GCST000001', geneId=None, projectId='GCST', studyType='gwas', traitFromSource='Test trait', traitFromSourceMappedIds=['EFO_0000001'], pubmedId='111', nSamples=1000, nCases=300,
               cohorts=['TestCohort'], discoverySamples=[{'sampleSize': 1000, 'ancestry': 'European'}], ldPopulationStructure=[{'ldPopulation': 'nfe', 'relativeSampleSize': 1.0}], hasSumstats=True,
               diseaseIds=['MONDO_0000001'], biosampleId=None),
          dict(studyId='FINNGEN_R12_TEST', geneId=None, projectId='FINNGEN_R12', studyType='gwas', traitFromSource='Finnish trait', traitFromSourceMappedIds=['EFO_0000002'], pubmedId=None, nSamples=500,
               nCases=None, cohorts=['FinnGen'], discoverySamples=[], ldPopulationStructure=[{'ldPopulation': 'fin', 'relativeSampleSize': 1.0}], hasSumstats=True, diseaseIds=[], biosampleId=None),
          dict(studyId='GTEX_EQTL_X', geneId='ENSG_TESTG', projectId='GTEx-v10', studyType='eqtl', traitFromSource='ENSG_TESTG', traitFromSourceMappedIds=[], pubmedId=None, nSamples=300, nCases=None,
               cohorts=['GTEx'], discoverySamples=[], ldPopulationStructure=[{'ldPopulation': 'nfe', 'relativeSampleSize': 1.0}], hasSumstats=False, diseaseIds=[], biosampleId='EFO_0005292')]
    pq.write_table(pa.Table.from_pylist(st), os.path.join(out, 'study', 'part-0.parquet'))
    wrong = 'A' if b('1', 440) != 'A' else 'C'
    vid = lambda pos, r=None, a=None: f'1_{pos}_{r or b("1", pos)}_{a or alt("1", pos)}'

    def member(variant, pip, beta=None, se=None, mant=None, exp=None, r2=None):
        return dict(is95CredibleSet=pip >= 0.1, is99CredibleSet=True, logBF=3.0, posteriorProbability=pip, variantId=variant, pValueMantissa=mant, pValueExponent=exp, beta=beta, standardError=se, r2Overall=r2)

    def cset(sl, study, lead, typ, members, method='SuSiE-inf', conf='SuSiE fine-mapped credible set with out-of-sample LD', idx=1):
        return dict(studyLocusId=sl, studyId=study, variantId=lead, chromosome='1', position=int(lead.split('_')[1]), beta=0.5, standardError=0.1, pValueMantissa=1.0, pValueExponent=-9,
                    effectAlleleFrequencyFromSource=0.25, finemappingMethod=method, credibleSetIndex=idx, credibleSetlog10BF=12.0, purityMeanR2=0.9, sampleSize=1000, confidence=conf, studyType=typ,
                    isTransQtl=False, qualityControls=[], locus=members)
    sets = [cset('SL1', 'GCST000001', vid(402), 'gwas', [member(vid(402), 0.8, 0.5, 0.1, 1.0, -9, 1.0), member(vid(410), 0.15, 0.4, 0.1, 2.0, -7, 0.8), member(vid(420), 0.05, 0.1, 0.1, 5.0, -3, 0.5),
                                                          member(vid(440, wrong, 'T' if wrong != 'T' else 'G'), 0.3, 0.2, 0.1, 3.0, -5, 0.7), member(vid(430), 0.2, -0.3, 0.1, 1.0, -400, 0.6)]),
            cset('SL2', 'GTEX_EQTL_X', vid(402), 'eqtl', [member(vid(402), 0.6, -0.8, 0.2, 4.0, -8, 1.0), member(vid(410), 0.4, -0.2, 0.2, 1.0, -4, 0.7)], method='SuSie', conf='SuSiE fine-mapped credible set with in-sample LD'),
            cset('SL3', 'FINNGEN_R12_TEST', vid(450), 'gwas', [member(vid(450), 0.04), member(vid(460), 0.03)], method='PICS', conf='PICS fine-mapped credible set based on reported top hit')]
    pq.write_table(pa.Table.from_pylist(sets), os.path.join(out, 'credible_set', 'part-0.parquet'))
    l2g = [dict(studyLocusId='SL1', geneId='ENSG_TESTG', score=0.9, features=[], shapBaseValue=0.1), dict(studyLocusId='SL1', geneId='ENSG_OTHER', score=0.3, features=[], shapBaseValue=0.1),
           dict(studyLocusId='SL3', geneId='ENSG_UNKNOWN', score=0.6, features=[], shapBaseValue=0.1)]
    pq.write_table(pa.Table.from_pylist(l2g), os.path.join(out, 'l2g_prediction', 'part-0.parquet'))
    return dict(k402=f'1:402:{b("1", 402)}:{alt("1", 402)}', k410=f'1:410:{b("1", 410)}:{alt("1", 410)}', k430=f'1:430:{b("1", 430)}:{alt("1", 430)}')


def mini_opentargets_gene(out):
    """Mini Open Targets 26.09 gene-level tables: target_essentiality, evidence_impc, evidence_crispr, evidence_crispr_screen, evidence_gene_burden."""
    import pyarrow as pa
    import pyarrow.parquet as pq
    for d in ('target_essentiality', 'evidence_impc', 'evidence_crispr', 'evidence_crispr_screen', 'evidence_gene_burden'):
        os.makedirs(os.path.join(out, d), exist_ok=True)
    scr = lambda *effects: [{'depmapId': f'ACH-{i}', 'cellLineName': f'L{i}', 'diseaseFromSource': 'x', 'diseaseCellLineId': None, 'mutation': None, 'geneEffect': e, 'expression': 1.0} for i, e in enumerate(effects)]
    ess = [dict(targetId='ENSG_TESTG', isEssential=True, depMapEssentiality=[{'tissueId': 'UBERON_1', 'tissueName': 'lung', 'screens': scr(-1.5, -1.2, -0.3, -1.1)},
                                                                           {'tissueId': 'UBERON_2', 'tissueName': 'skin of body', 'screens': scr(-1.3, -1.4)}]),
           dict(targetId='ENSG_TESTH', isEssential=False, depMapEssentiality=[{'tissueId': 'UBERON_1', 'tissueName': 'lung', 'screens': scr(0.1, 0.2, -0.2)}])]
    pq.write_table(pa.Table.from_pylist(ess), os.path.join(out, 'target_essentiality', 'part-0.parquet'))

    def impc(gene, comp, bg, disease, score, phenos, model='MGI:1'):
        return dict(targetId=gene, biologicalModelAllelicComposition=comp, biologicalModelGeneticBackground=bg, biologicalModelId=model, diseaseId=disease, diseaseFromSource=disease, resourceScore=score,
                    diseaseModelAssociatedModelPhenotypes=[{'id': f'MP:{i}', 'label': lab} for i, lab in enumerate(phenos)], literature=['111'], directionOnTrait='risk')
    rows = [impc('ENSG_TESTG', 'Testg<tm1A>/Testg<tm1A>', 'involves: C57BL/6J', 'MONDO_1', 40.0, ['short limbs', 'abnormal skull']),
            impc('ENSG_TESTG', 'Testg<tm1A>/Testg<tm1A>', 'involves: C57BL/6J', 'MONDO_2', 55.0, ['short limbs', 'embryonic lethality']),
            impc('ENSG_TESTG', 'Testg<tm1A>/Testg<+>', 'involves: C57BL/6J', 'MONDO_1', 20.0, ['small ears'], model='MGI:2'),
            impc('ENSG_UNKNOWN', 'H4c1<em1(IMPC)Ccpcz> hom early', 'C57BL/6N', 'MONDO_3', 30.0, ['abnormal gait'], model='IMPC:3')]
    pq.write_table(pa.Table.from_pylist(rows), os.path.join(out, 'evidence_impc', 'part-0.parquet'))
    cr = [dict(id='c1', targetId='ENSG_TESTG', diseaseId='EFO_1', diseaseFromSource='Test cancer', resourceScore=0.01, literature=['222'], datasourceId='crispr', publicationDate='2020-01-01',
               diseaseCellLines=[{'id': 'ACH-1', 'name': 'CELL1', 'tissue': 'lung', 'tissueId': 'UBERON_1'}], targetFromSourceId='TESTG', targetFromSource='TESTG', diseaseFromSourceMappedId='EFO_1', score=1.0)]
    pq.write_table(pa.Table.from_pylist(cr), os.path.join(out, 'evidence_crispr', 'part-0.parquet'))
    cs = [dict(id='s1', targetId='ENSG_TESTH', studyId='STUDY1', projectId='crispr_brain', crisprScreenLibrary='Genome-wide CRISPRi-v2', contrast='neuron', cellType='Glutamatergic Neuron', geneticBackground='iPSC',
               studyOverview='Description', diseaseId='MONDO_1', diseaseFromSource='ageing', literature=['333'], log2FoldChangeValue=4.5, resourceScore=0.005, statisticalTestTail='upper tail', datasourceId='crispr_screen')]
    pq.write_table(pa.Table.from_pylist(cs), os.path.join(out, 'evidence_crispr_screen', 'part-0.parquet'))
    gb = [dict(id='b1', targetId='ENSG_TESTG', diseaseId='EFO_2', diseaseFromSource='Test trait', oddsRatio=3.2, oddsRatioConfidenceIntervalLower=2.0, oddsRatioConfidenceIntervalUpper=5.1, beta=None,
               betaConfidenceIntervalLower=None, betaConfidenceIntervalUpper=None, pValueMantissa=2.0, pValueExponent=-9, ancestry='European', studySampleSize=400000, studyCases=1200, allelicRequirements=['monoallelic'],
               projectId='Genebass', cohortId='UKB', statisticalMethod='SKAT-O', releaseVersion='r1', studyCasesWithQualifyingVariants=30, sex=['male', 'female'], directionOnTrait='risk'),
          dict(id='b2', targetId='ENSG_TESTH', diseaseId='EFO_3', diseaseFromSource='Another trait', oddsRatio=None, oddsRatioConfidenceIntervalLower=None, oddsRatioConfidenceIntervalUpper=None, beta=-0.2,
               betaConfidenceIntervalLower=-0.4, betaConfidenceIntervalUpper=-0.1, pValueMantissa=0.5, pValueExponent=-1, ancestry='European', studySampleSize=400000, studyCases=None, allelicRequirements=[],
               projectId='Genebass', cohortId='UKB', statisticalMethod='burden', releaseVersion='r1', studyCasesWithQualifyingVariants=None, sex=[], directionOnTrait='protect')]
    pq.write_table(pa.Table.from_pylist(gb), os.path.join(out, 'evidence_gene_burden', 'part-0.parquet'))


def mini_events(out):
    """Mini versions of the gnomAD SV, gnomAD CNV, DGV and HGSVC3 event files (real layouts, GRCh38-style chr names on the synthetic contigs)."""
    import gzip
    os.makedirs(out, exist_ok=True)
    hdr8 = '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO'
    sv = [('chr1', 1000, 'gnomAD-SV_v3_DEL_chr1_aaa', '<DEL>', 'PASS', 'ALGORITHMS=manta,depth;CHR2=chr1;END=2000;EVIDENCE=PE,RD;SVLEN=1000;SVTYPE=DEL;AN=2000;AC=40;AF=0.02;N_BI_GENOS=1000;N_HOMALT=2;AF_afr=0.05;AF_nfe=0.01;PREDICTED_LOF=GENEA'),
          ('chr1', 1500, 'gnomAD-SV_v3_DUP_chr1_bbb', '<DUP>', 'HIGH_NCR', 'ALGORITHMS=depth;CHR2=chr1;END=2500;SVLEN=1000;SVTYPE=DUP;AN=2000;AC=4;AF=0.002;N_BI_GENOS=1000;LOW_CONFIDENCE_REPETITIVE_LARGE_DUP'),
          ('chr1', 100, 'gnomAD-SV_v3_BND_chr1_ccc', '<BND>', 'UNRESOLVED', 'ALGORITHMS=manta;CHR2=chr2;END=100;END2=900;SVLEN=-1;SVTYPE=BND;AN=2000;AC=100;AF=0.05;N_BI_GENOS=1000'),
          ('chr1', 2000, 'gnomAD-SV_v3_CNV_chr1_ddd', '<CNV>', 'PASS', 'ALGORITHMS=depth;CHR2=chr1;END=2400;MULTIALLELIC;SVLEN=400;SVTYPE=CNV;CN_NUMBER=900;CN_NONREF_FREQ=0.3;CN_FREQ=0.1,0.6,0.3'),
          ('chrUn_gl000220', 5, 'gnomAD-SV_v3_INS_chrUn_eee', '<INS>', 'PASS', 'ALGORITHMS=melt;END=5;SVLEN=300;SVTYPE=INS;AN=100;AC=1;AF=0.01;N_BI_GENOS=50'),
          ('chr2', 300, 'gnomAD-SV_v3_INS_chr2_fff', '<INS>', 'PASS', 'ALGORITHMS=melt;END=300;SVLEN=300;SVTYPE=INS;AN=100;AC=1;AF=0.01;N_BI_GENOS=50')]
    with gzip.open(os.path.join(out, 'gnomad.v4.1.sv.sites.vcf.gz'), 'wt') as f:
        f.write('##fileformat=VCFv4.2\n' + hdr8 + '\n' + ''.join('\t'.join([c, str(p_), i, 'N', a, '.', fl, inf]) + '\n' for c, p_, i, a, fl, inf in sv))
    cnv = [('chr1', 700, 'variant_is_1_1__DEL', '<DEL>', 'PASS', 'END=900;SVLEN=200;SVTYPE=DEL;POSMIN=690;POSMAX=700;ENDMIN=900;ENDMAX=910;Genes=GENEA;N_EXN_VAR=2;SC=3;SN=1000;SF=0.003;SF_nfe=0.004'),
           ('chr1', 700, 'variant_is_1_2__DUP', '<DUP>', 'PASS', 'END=950;SVLEN=250;SVTYPE=DUP;POSMIN=700;POSMAX=700;ENDMIN=950;ENDMAX=950;Genes=GENEA;SC=1;SN=1000;SF=0.001')]
    with gzip.open(os.path.join(out, 'gnomad.v4.1.cnv.all.vcf.gz'), 'wt') as f:
        f.write('##fileformat=VCFv4.2\n' + hdr8 + '\n' + ''.join('\t'.join([c, str(p_), i, 'N', a, '.', fl, inf]) + '\n' for c, p_, i, a, fl, inf in cnv))
    dh = ['variantaccession', 'chr', 'start', 'end', 'varianttype', 'variantsubtype', 'reference', 'pubmedid', 'method', 'platform', 'mergedvariants', 'supportingvariants', 'mergedorsample', 'frequency',
          'samplesize', 'observedgains', 'observedlosses', 'cohortdescription', 'genes', 'samples']
    dr = [['nsv1', '1', '1000', '2000', 'CNV', 'deletion', 'Study_A', '111', 'Sequencing', '', '', 'nssv1', 'M', '', '3202', '0', '5', '', '', ''],
          ['nsv2', '1', '1100', '2100', 'CNV', 'gain', 'Study_B', '222', 'Oligo aCGH', '', '', '', 'M', '0.5', '100', '2', '0', '', '', ''],
          ['nssv9', '1', '1000', '2000', 'CNV', 'deletion', 'Study_A', '111', 'Sequencing', '', '', '', 'S', '', '', '', '', '', '', 'HG0001'],
          ['nsv3', '2', '10', '50', 'OTHER', 'unusual thing', 'Study_C', '333', 'PCR', '', '', '', 'M', '', '10', '', '', '', '', ''],
          ['nsv4', '1', '500', '300', 'CNV', 'duplication', 'Study_D', '444', 'Sequencing', '', '', '', 'M', '', '10', '', '', '', '', '']]
    with open(os.path.join(out, 'GRCh38_hg38_variants_2025-12-01.txt'), 'w') as f:
        f.write('\t'.join(dh) + '\n' + ''.join('\t'.join(r) + '\n' for r in dr))
    hv = '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tH1\tH2\tH3'
    ins = [('chr1', 1200, 'chr1-1200-DEL-300', 'DEL', -300, 'ID=chr1-1200-DEL-300;SVTYPE=DEL;SVLEN=-300', ['0|1', '0|0', '1|.']),
           ('chr1', 1300, 'chr1-1300-INS-150', 'INS', 150, 'ID=chr1-1300-INS-150;SVTYPE=INS;SVLEN=150', ['.|.', '.|.', '.|.'])]
    with gzip.open(os.path.join(out, 'hgsvc3_insdel.vcf.gz'), 'wt') as f:
        f.write('##fileformat=VCFv4.2\n' + hv + '\n' + ''.join('\t'.join([c, str(p_), i, 'N', f'<{t}>', '.', '.', inf, 'GT'] + g) + '\n' for c, p_, i, t, l, inf, g in ins))
    with gzip.open(os.path.join(out, 'hgsvc3_inv.vcf.gz'), 'wt') as f:
        f.write('##fileformat=VCFv4.2\n' + hv + '\n' + '\t'.join(['chr1', '600', 'chr1-600-INV-200', 'N', '<INV>', '.', '.', 'ID=chr1-600-INV-200;SVTYPE=INV;SVLEN=200', 'GT', '1|1', '0|0', '0|1']) + '\n')


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
    write_all.clinpgx = mini_clinpgx(os.path.join(out, 'pharmgkb'))
    write_all.mavedb = mini_mavedb(os.path.join(out, 'mavedb'))
    write_all.opentargets = mini_opentargets(os.path.join(out, 'opentargets'))
    mini_opentargets_gene(os.path.join(out, 'opentargets'))
    mini_events(os.path.join(out, 'events'))
    return out


if __name__ == '__main__':
    import sys
    print('written', write_all(sys.argv[1] if len(sys.argv) > 1 else 'sources'))
