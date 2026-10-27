"""Table schemas and controlled vocabulary of the High-Value Human Variant Reference Index.

Design rules (they are tested):
  * Source records are immutable and stay in the acquired files. An evidence row points to its record (release_id, accession,
    locator) and carries the SHA-256 of the raw record text, so a reader can re-open the record and prove it is unchanged.
  * Identity: a small variant is its normalised GRCh38 key CHROM:POS:REF:ALT. rsIDs, ClinVar IDs, CA IDs and OT IDs are xrefs.
    Structural events are never merged by overlap; each source event is its own entity and equivalence is a separate, explained edge.
  * One evidence row is one assertion or measurement by one source about one subject at one level
    (variant | haplotype | gene | region | event | assay | study). Gene-level evidence stays on the gene; it reaches a variant only
    through an explicit relation row and is displayed as gene context, never as variant evidence.
  * Seven parallel axes; no composite score and no ordinal tier. A priority list is a view (a query) over this table.
  * Negative, null, uncertain and conflicting assertions are rows like any other (effect_assertion). Missing evidence is the absence
    of a row and means 'not assessed here', never 'no effect'.
  * Evidence rows that share a lineage_group descend from one underlying study or dataset and are not independent replication. A key is shared across databases where they
    describe the same thing (pmid:<PubMed id>, clinvar:variation:<VariationID>); see 'Lineage keys' in SCHEMA.md.
"""
import pyarrow as pa

from . import SCHEMA_VERSION

PRIMARY_CONTIGS = tuple([str(i) for i in range(1, 23)] + ['X', 'Y', 'MT'])
OTHER_PARTITION = 'NA'          # not '_': Arrow dataset discovery silently skips paths that start with '_' or '.'


def partition_of(chrom):
    """Directory under evidence/<level>/ that holds the rows of a contig: the contig itself on the primary assembly, 'NA' for every other contig (alternate loci,
    unlocalised and unplaced scaffolds, decoys). A row keeps its own chrom, so a reader that is asked for such a contig filters on subject_id or chrom inside 'NA'."""
    return chrom if chrom in PRIMARY_CONTIGS else OTHER_PARTITION


AXES = {
    'clinical_validity_actionability': 'Clinical classification, gene-disease validity, dosage sensitivity, clinical actionability and pharmacogenomic guidance, as asserted by curators or expert panels.',
    'empirical_functional_support': 'Measurements of the effect of the specific allele, or of a gene/region perturbation, in an assay: deep mutational scanning, MPRA, knockout and CRISPR screens, model-organism phenotypes.',
    'mechanistic_confidence': 'How well a particular molecular mechanism (loss of function, splicing, regulatory, dosage, protein-level) is supported for the allele; competing mechanisms stay as separate rows. Predictions are marked as such.',
    'effect_magnitude_direction': 'The size and sign of an effect with its units and uncertainty: beta or odds ratio for an effect allele, assay score and direction, allele or haplotype function.',
    'association_replication': 'Statistical association and its replication: p-value, sample metadata, independent studies. Rows sharing a lineage_group count once.',
    'fine_mapping_probability': 'Posterior inclusion probability, credible-set membership, locus-to-gene and colocalisation scores.',
    'applicability_context': 'Where the evidence applies: tissue or cell type, ancestry, phenotype, sex, age, germline or somatic origin, population frequency.',
}
LEVELS = ('variant', 'haplotype', 'gene', 'region', 'event', 'assay', 'study')
EFFECT_ASSERTIONS = {
    'effect': 'the source asserts or measures an effect (pathogenic, associated, functionally abnormal, loss of function)',
    'no_effect': 'the source asserts or measures no effect (benign, null association, functionally normal)',
    'uncertain': 'the source states uncertainty (VUS, inconclusive, low confidence)',
    'conflicting': 'the source itself reports conflicting interpretations',
    'not_applicable': 'context or descriptive row with no effect assertion (frequency, annotation)',
}
MAPPING_STATUS = {
    'mapped': 'resolved to an index entity of the stated level (for variants: a normalised key whose REF equals the reference)',
    'ambiguous': 'several candidate entities remain (for example a protein-level assay variant with several codon alternatives)',
    'unmapped': 'no entity could be derived from the record (no alleles, no position, unsupported identifier)',
    'rejected': 'the record contradicts the reference (REF mismatch) or its strand/allele orientation cannot be resolved',
}
MATCH_STATUS = {
    'EXACT_MATCH': 'the sample has a call of exactly this normalised allele that passes the quality policy',
    'FILTERED_CALL': 'the allele is called but fails FILTER/quality thresholds; presence is not asserted',
    'NO_CALL': 'the site is in the callset but the genotype is missing or incomplete',
    'VERIFIED_REFERENCE': 'reference-confidence evidence (gVCF block, callable region or reads) shows the allele is not carried',
    'OTHER_ALLELE_AT_SITE': 'a different allele is called at the same site and the indexed allele is not among the called alleles',
    'DELETED_BY_UPSTREAM': 'the haplotype is deleted at this position by a call elsewhere (spanning-deletion allele)',
    'ABSENT_FROM_CALLSET': 'no record and no reference-confidence evidence: the genotype is unknown, not reference',
    'UNRESOLVED': 'the question cannot be answered from the input (phase, structure, repeat length, difficult region); the reasons say exactly what is missing',
}
EVIDENCE_KINDS = ('direct_experimental', 'perturbation_experimental', 'population_observation', 'statistical_association', 'curated_expert_classification',
                  'computational_prediction', 'reference_annotation', 'experimental_structure', 'predicted_structure', 'aggregated_integration',
                  'experimental_annotation', 'case_observation')

S = pa.string()
LS = pa.list_(pa.string())

SOURCE_RELEASE = pa.schema([
    ('release_id', S), ('source_id', S), ('target_id', S), ('version', S), ('file', S), ('size_bytes', pa.int64()), ('sha256', S),
    ('url', S), ('retrieved_utc', S), ('licence', S), ('licence_class', S), ('manifest_line', pa.int32()),
])

ENTITY_VARIANT = pa.schema([
    ('variant_key', S), ('chrom', S), ('pos', pa.int32()), ('ref', S), ('alt', S), ('vtype', S),
])

XREF = pa.schema([
    ('namespace', S), ('identifier', S), ('entity_level', S), ('entity_id', S), ('release_id', S), ('relation', S),
])

ENTITY_EVENT = pa.schema([
    ('event_id', S), ('source_id', S), ('release_id', S), ('accession', S), ('event_type', S), ('chrom', S), ('start', pa.int32()), ('end', pa.int32()),
    ('start_ci_lo', pa.int32()), ('start_ci_hi', pa.int32()), ('end_ci_lo', pa.int32()), ('end_ci_hi', pa.int32()), ('svlen', pa.int32()),
    ('copy_number', pa.int32()), ('precise', pa.bool_()), ('population_af', pa.float64()), ('n_samples', pa.int32()), ('origin', S),
    ('raw_type', S), ('assembly', S), ('extra_json', S),
])

ENTITY_HAPLOTYPE = pa.schema([
    ('haplotype_id', S), ('system', S), ('gene_symbol', S), ('allele_name', S), ('is_reference', pa.bool_()), ('structural', pa.bool_()),
    ('release_id', S), ('accession', S), ('function', S), ('allele_level', S),
])

HAPLOTYPE_DEFINITION = pa.schema([
    ('haplotype_id', S), ('variant_key', S), ('hgvs', S), ('release_id', S), ('definition_scope', S),   # scope: defining | tolerated (the allele may carry this alternative or the reference base) | unresolved:<reason>
])

ENTITY_GENE = pa.schema([
    ('gene_id', S), ('symbol', S), ('hgnc_id', S), ('chrom', S), ('start', pa.int32()), ('end', pa.int32()), ('strand', S),
    ('mane_select', S), ('biotype', S), ('release_id', S),
])

ENTITY_REGION = pa.schema([
    ('region_id', S), ('kind', S), ('chrom', S), ('start', pa.int32()), ('end', pa.int32()), ('name', S), ('gene_ids', LS), ('release_id', S),
])

ENTITY_STUDY = pa.schema([
    ('study_id', S), ('source_id', S), ('lineage_group', S), ('study_type', S), ('pubmed', S), ('trait_ids', LS), ('trait_text', S),
    ('n_total', pa.int32()), ('n_cases', pa.int32()), ('ancestry', S), ('cohorts', S), ('biosample', S), ('has_sumstats', pa.bool_()), ('release_id', S),
])

ENTITY_ASSAY = pa.schema([
    ('assay_id', S), ('source_id', S), ('kind', S), ('target_gene', S), ('model_system', S), ('tissue_or_cell', S), ('perturbation', S),
    ('dose_context', S), ('readout', S), ('score_direction', S), ('lineage_group', S), ('release_id', S),
])

GENE_SET_MEMBER = pa.schema([
    ('set_id', S), ('source_id', S), ('kind', S), ('name', S), ('gene_id', S), ('release_id', S),      # one row per (set, gene): Reactome pathways today; annotation, never evidence
])

VARIANT_GENE = pa.schema([
    ('variant_key', S), ('gene_id', S), ('relation', S), ('consequence', S), ('release_id', S),
])

EVIDENCE = pa.schema([
    ('evidence_id', S), ('source_id', S), ('release_id', S), ('accession', S), ('locator', S), ('raw_sha256', S),
    ('level', S), ('subject_id', S), ('chrom', S), ('pos', pa.int32()), ('end', pa.int32()),
    ('axis', S), ('evidence_kind', S), ('claim_type', S), ('claim_value', S), ('claim_num', pa.float64()),
    ('effect_assertion', S), ('direction', S),
    ('condition_ids', LS), ('condition_text', S), ('gene_ids', LS),
    ('origin', S), ('inheritance', S), ('penetrance', S), ('tissue', S), ('cell_type', S), ('model_system', S), ('perturbation', S),
    ('dose', S), ('ancestry', S), ('sample_size', pa.int32()), ('n_cases', pa.int32()),
    ('effect_allele', S), ('other_allele', S), ('beta', pa.float64()), ('se', pa.float64()), ('ci_low', pa.float64()), ('ci_high', pa.float64()),
    ('odds_ratio', pa.float64()), ('p_value', pa.float64()), ('pip', pa.float64()), ('score', pa.float64()), ('score_unit', S), ('eaf', pa.float64()),
    ('review_status', S), ('review_stars', pa.int8()), ('n_submitters', pa.int32()), ('assertion_criteria', S), ('evaluated', S), ('curator', S),
    ('lineage_group', S), ('lineage_note', S), ('mapping_status', S), ('mapping_reason', S), ('unmapped_key', S), ('extra_json', S),
])

TABLES = {
    'source_release': SOURCE_RELEASE, 'entity_variant': ENTITY_VARIANT, 'xref': XREF, 'entity_event': ENTITY_EVENT,
    'entity_haplotype': ENTITY_HAPLOTYPE, 'haplotype_definition': HAPLOTYPE_DEFINITION, 'entity_gene': ENTITY_GENE,
    'entity_region': ENTITY_REGION, 'entity_study': ENTITY_STUDY, 'entity_assay': ENTITY_ASSAY, 'variant_gene': VARIANT_GENE, 'gene_set_member': GENE_SET_MEMBER, 'evidence': EVIDENCE,
}


def empty_table(name):
    return TABLES[name].empty_table()


def rows_to_table(name, rows):
    """List of dicts -> Table with the named schema (missing keys become null; unknown keys raise)."""
    sch = TABLES[name]
    unknown = {k for r in rows for k in r} - set(sch.names)
    if unknown:
        raise KeyError(f'{name}: unknown column(s) {sorted(unknown)}')
    return pa.Table.from_pylist([{n: r.get(n) for n in sch.names} for r in rows], schema=sch)


def validate_evidence_row(r):
    """Return a list of problems with one evidence row (empty when sound). Enforces the vocabulary and the level/axis rules."""
    p = []
    if r.get('level') not in LEVELS:
        p.append(f"level {r.get('level')!r}")
    if r.get('axis') not in AXES:
        p.append(f"axis {r.get('axis')!r}")
    if r.get('evidence_kind') not in EVIDENCE_KINDS:
        p.append(f"evidence_kind {r.get('evidence_kind')!r}")
    if r.get('effect_assertion') not in EFFECT_ASSERTIONS:
        p.append(f"effect_assertion {r.get('effect_assertion')!r}")
    if r.get('mapping_status') not in MAPPING_STATUS:
        p.append(f"mapping_status {r.get('mapping_status')!r}")
    for k in ('evidence_id', 'source_id', 'release_id', 'accession', 'locator', 'raw_sha256', 'claim_type'):
        if not r.get(k):
            p.append(f'missing {k}')
    if r.get('mapping_status') == 'mapped' and not r.get('subject_id'):
        p.append('mapped row without subject_id')
    if r.get('mapping_status') in ('unmapped', 'rejected') and r.get('subject_id'):
        p.append(f"{r.get('mapping_status')} row with a subject_id")
    if r.get('mapping_status') != 'mapped' and not r.get('mapping_reason'):
        p.append('unmapped/ambiguous/rejected row without mapping_reason')
    if r.get('level') == 'variant' and r.get('mapping_status') == 'mapped' and r.get('subject_id', '').count(':') != 3:
        p.append('variant subject_id is not a CHROM:POS:REF:ALT key')
    if r.get('p_value') is not None and not (0 <= r['p_value'] <= 1):
        p.append('p_value outside [0,1]')
    if r.get('pip') is not None and not (0 <= r['pip'] <= 1.0000001):
        p.append('pip outside [0,1]')
    if r.get('level') != 'variant' and r.get('axis') == 'fine_mapping_probability' and r.get('claim_type') == 'pip':
        p.append('pip on a non-variant subject')
    return p


def schema_markdown():
    """SCHEMA.md body, generated from the definitions so documentation cannot drift from the code."""
    L = [f'# Index schema (version {SCHEMA_VERSION})', '', 'Generated by `python3 -m vindex.schema > SCHEMA.md` from `vindex/schema.py`. Do not edit by hand.', '',
         '## Axes (parallel; there is no composite score and no ordinal tier)', '', '| Axis | Meaning |', '|---|---|']
    L += [f'| `{k}` | {v} |' for k, v in AXES.items()]
    L += ['', '## Levels', '', ', '.join(f'`{x}`' for x in LEVELS), '', '## effect_assertion', '', '| Value | Meaning |', '|---|---|']
    L += [f'| `{k}` | {v} |' for k, v in EFFECT_ASSERTIONS.items()]
    L += ['', '## mapping_status', '', '| Value | Meaning |', '|---|---|'] + [f'| `{k}` | {v} |' for k, v in MAPPING_STATUS.items()]
    L += ['', '## Match status (matcher output)', '', '| Value | Meaning |', '|---|---|'] + [f'| `{k}` | {v} |' for k, v in MATCH_STATUS.items()]
    L += ['', '## Physical layout', '',
          'One Parquet file per source in each entity table (`<table>/<source>.parquet`); `source_release.parquet`, `index_manifest.json` (build facts) and `table_digests.json` (SHA-256 of every file) at the top. '
          'Evidence is partitioned by level and contig: `evidence/<level>/<partition>/<source_id>.parquet`, and `entity_variant/<partition>.parquet` follows the same rule. The partition of a contig is the '
          f'contig itself for {", ".join(PRIMARY_CONTIGS)} and `{OTHER_PARTITION}` for every other contig (alternate loci, unlocalised and unplaced scaffolds, decoys) and for rows that carry no placement. '
          'A row keeps its own `chrom`, so a reader that wants one alternate locus filters on `subject_id` or `chrom` inside the `NA` partition. A mapped variant row is always filed under the partition '
          'of its key\'s contig (the builder refuses anything else and the validator checks it).']
    L += ['', '## Lineage keys', '',
          'Rows that share a `lineage_group` descend from one study or record and count once in the independence tally. A key is shared across databases where they describe the same thing: '
          '`pmid:<PubMed id>` (a paper, whichever of the GWAS Catalog, Open Targets, ClinPGx, MaveDB or ClinVar observations describes it; several GWAS Catalog studies of one paper count once), '
          '`clinvar:variation:<VariationID>` (the ClinVar aggregate, its submissions and the ClinGen expert-panel classification of that variation) and `gwas_study:<accession>` (a GWAS Catalog study with no '
          'PubMed id). Where no other database describes the same record the key is the source\'s own (`cpic_guideline:`, `clinpgx_guideline:`, `strchive:`, `gnomad_v4.1_joint`, ...). Two sources that '
          'describe one record under different keys are not linked (the CPIC guidelines and ClinPGx\'s annotations of them); `make_reports.py` lists the groups that more than one source writes.']
    L += ['', '## Tables']
    for name, sch in TABLES.items():
        L += ['', f'### `{name}`', '', '| Column | Type |', '|---|---|'] + [f'| `{f.name}` | {f.type} |' for f in sch]
    return '\n'.join(L) + '\n'


if __name__ == '__main__':
    print(schema_markdown(), end='')
