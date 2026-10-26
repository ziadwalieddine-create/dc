# Evidence profiles

*Generated 2026-10-26 15:25 UTC by `catalogue/build_catalogue.py` from `catalogue/sources/*.json`, `acquisition/targets.json`, `acquisition/manifest.jsonl`, `acquisition/validation_report.json` and `acquisition/compat_report.json`. Edit the sources, not this file.*

Every source record carries an `evidence_profile` that keeps apart what the resource is (a measurement, a perturbation, a population observation, an association, a curated classification, a prediction, a reference annotation), the resolution at which its records refer to alleles, the mechanisms it informs and on what basis, how strong the evidence is, how far a mechanistic reading can be trusted, what phenotype it implies and what absence means. These fields are independent: a curated classification can be strong and carry no mechanism; a strong association says nothing about the causal variant; a gene-level knockout or constraint metric is never allele-specific evidence; the absence of a study is never recorded as the absence of an effect. Hypotheses and competing explanations are listed per mechanism axis in `MECHANISM_MATRIX.md`.

Rules enforced by `build_catalogue.py --check`: a gene- or region-level resource cannot be direct evidence for a specific allele or small variant; a prediction basis needs the `computational_prediction` kind; every mechanism basis needs a matching evidence kind; a prediction resource cannot also carry direct or curated mechanism evidence unless it aggregates other sources; reference, catalogue and documentation resources provide context only; mechanism entries cannot name variant classes the record does not list; every overlap between sources is declared once, and identical, derived or conflicting records must say what is shared or contested.

## Evidence kinds

| Kind | Meaning | Sources |
|---|---|---|
| `direct_experimental` | A measurement of the effect of specific alleles in an assay (for example deep mutational scanning, MPRA, saturation genome editing) | `ensembl_vep`, `mavedb`, `mpravardb`, `proteingym` |
| `perturbation_experimental` | A gene- or region-level perturbation in cells or model organisms (knockout, knockdown, CRISPRi/a, overexpression); not allele-specific unless stated | `alliance_genome`, `biogrid_interactions`, `biogrid_orcs`, `depmap`, `impc`, `opentargets` |
| `population_observation` | Allele or genotype observations in cohorts (frequency, homozygosity, constraint derived from them); indirect evidence about function | `1000genomes`, `dbnsfp`, `dbvar`, `dgv`, `ensembl_vep`, `gnomad_constraint`, `gnomad_ld`, `gnomad_sv_cnv`, `gnomad_variants`, `hgsvc3`, `ldlink`, `trexplorer` |
| `statistical_association` | A statistical association between variants and a trait or molecular phenotype (GWAS, QTL, fine-mapping, colocalisation); association is not mechanism | `eqtl_catalogue`, `gtex`, `gwas_catalog`, `opentargets`, `pgs_catalog` |
| `curated_expert_classification` | A classification or assertion made by curators or expert panels from other evidence (clinical significance, gene-disease validity, dosage sensitivity, star-allele function) | `alliance_genome`, `biogrid_interactions`, `clingen`, `clinvar`, `cpic`, `dbvar`, `intact`, `ipd_imgt_hla`, `ipd_kir`, `mitomap`, `opentargets`, `pharmcat`, `pharmgkb`, `pharmvar`, `reactome`, `strchive`, `string_db`, `uniprot` |
| `computational_prediction` | A model or algorithm output; a hypothesis about effect, never a measurement | `alphamissense`, `cadd`, `dbnsfp`, `ensembl_vep`, `gnomad_sv_cnv`, `gnomad_variants`, `opentargets`, `string_db` |
| `reference_annotation` | Reference sequence, gene models, transcript sets, region annotations and identifier maps | `ensembl_genomes`, `gencode_mane`, `giab_stratifications`, `interpro`, `ipd_imgt_hla`, `ipd_kir`, `pdbe_sifts`, `refseq`, `trexplorer`, `uniprot` |
| `experimental_structure` | Atomic structures or residue coverage determined experimentally | `pdbe_sifts` |
| `predicted_structure` | Structure models predicted computationally | `alphafold_db` |
| `aggregated_integration` | A resource that aggregates and harmonises records from other sources; evidence kinds of the parts persist | `alliance_genome`, `dbnsfp`, `dbvar`, `encode`, `ensembl_vep`, `hpa`, `interpro`, `mpravardb`, `opentargets`, `pharmgkb`, `proteingym`, `trexplorer` |
| `catalogue_metadata` | Descriptions of resources (catalogues of databases and tools), not biological evidence | `biotools`, `fairsharing`, `pgs_catalog` |
| `documentation` | Educational or background material with no variant-level data | `nhgri_education` |
| `experimental_annotation` | Assays of the biochemical activity of the genome (chromatin, binding, expression) that annotate regions or genes without measuring the effect of a particular allele | `encode`, `gtex`, `hpa` |
| `case_observation` | Observations of variants in affected individuals or families (case series, cohorts ascertained on phenotype, segregation); ascertained on phenotype, so frequency in cases is not a test of association | `decipher` |

## Allele resolution

The resolution at which a source's records refer to alleles. Gene-, region- and protein-level records describe a gene, interval or protein as a whole.

| Resolution | Meaning | Sources |
|---|---|---|
| `allele_specific` | Records refer to specific alleles (chromosome, position, reference, alternate or an exact protein change) | `1000genomes`, `alphamissense`, `cadd`, `clinvar`, `dbnsfp`, `ensembl_vep`, `eqtl_catalogue`, `gnomad_ld`, `gnomad_variants`, `gtex`, `gwas_catalog`, `ldlink`, `mavedb`, `mitomap`, `mpravardb`, `pgs_catalog`, `pharmgkb`, `proteingym` |
| `haplotype_level` | Records refer to haplotypes, star alleles or named alleles composed of several variants | `cpic`, `ipd_imgt_hla`, `ipd_kir`, `pharmcat`, `pharmvar` |
| `variant_class_level` | Records refer to classes of change (for example any loss-of-function change in a gene) | — |
| `gene_level` | Records refer to genes or proteins as a whole (constraint, knockout phenotypes, essentiality); they do not describe any particular allele | `alliance_genome`, `biogrid_interactions`, `biogrid_orcs`, `clingen`, `depmap`, `gnomad_constraint`, `hpa`, `impc`, `intact`, `opentargets`, `reactome`, `string_db` |
| `region_level` | Records refer to genomic intervals or loci (regulatory elements, CNV and repeat regions, difficult regions) | `dbvar`, `decipher`, `dgv`, `encode`, `ensembl_genomes`, `gencode_mane`, `giab_stratifications`, `gnomad_sv_cnv`, `hgsvc3`, `refseq`, `strchive`, `trexplorer` |
| `protein_residue_level` | Records refer to residues, domains or structures of a protein | `alphafold_db`, `interpro`, `pdbe_sifts`, `uniprot` |
| `none` | No variant-level or gene-level content | `biotools`, `fairsharing`, `nhgri_education` |

## Per-source profiles

### `1000genomes` · 1000 Genomes Project / IGSR (International Genome Sample Resource)

| | |
|---|---|
| Evidence kinds | population_observation |
| Allele resolution | allele_specific |
| Mechanism evidence | `population_constraint` · population · snv, indel · population_frequency<br>`structural_rearrangement` · population · sv · structural_variation_cnv |
| Biological context | Openly consented samples from people who declared themselves healthy, grouped by population label. The callsets describe germline variation in reference cohorts; no tissue, expression or phenotype data accompany the variants. |
| Evidence strength | Observational genotype calls with sample and population metadata. No effect size or per-variant uncertainty is attached, and call quality differs between SNV, indel and SV callsets. |
| Mechanistic confidence | None for mechanism. Presence of an allele in a healthy-declared cohort is indirect evidence of tolerance only and cannot separate a benign allele from one with recessive, late-onset or reduced-penetrance effects. |
| Phenotype implication | Observation in healthy-declared individuals argues against a high-penetrance, early-onset dominant effect of that allele; it does not exclude recessive, late-onset or reduced-penetrance effects. |
| Absence means | An allele missing from the callset was not observed in this cohort. Coverage, callset design and population representation decide what could be seen; absence is not evidence of rarity or harm. |

### `alliance_genome` · Alliance of Genome Resources (model-organism and human gene knowledge)

| | |
|---|---|
| Evidence kinds | perturbation_experimental, curated_expert_classification, aggregated_integration |
| Allele resolution | gene_level |
| Mechanism evidence | `gene_dosage` · direct · gene_level · perturbation_evidence (model-organism knockout and mutant phenotypes, orthology-mapped)<br>`outcome_association` · curated · gene_level · perturbation_evidence (curated disease-model annotations) |
| Biological context | Phenotypes and disease models of engineered or natural mutants in other species, mapped to human genes by orthology; the genotype is an organism-specific allele. |
| Evidence strength | Curated from model-organism databases; the evidence is as strong as the original experiment and weakened by orthology mapping. |
| Mechanistic confidence | Low to moderate for the gene's role in the model; none for a specific human allele. |
| Phenotype implication | Model-organism phenotype terms; supportive of a human phenotype only when concordant with human genetics. |
| Absence means | A gene without a model phenotype was not studied or has no ortholog in the model; it is not shown to be dispensable. |
| Overlaps (declared here) | independent_related · `impc` (mouse knockout phenotypes appear in both through MGI; count a mouse line once) |

### `alphafold_db` · AlphaFold Protein Structure Database

| | |
|---|---|
| Evidence kinds | predicted_structure |
| Allele resolution | protein_residue_level |
| Mechanism evidence | `protein_function` · context · protein_level · protein_structure_function |
| Biological context | Single-chain structure models of UniProt reference sequences. Ligands, binding partners, membranes, modifications and alternative conformations are not modelled. |
| Evidence strength | Computational model. Per-residue confidence (pLDDT) is the only uncertainty provided and measures confidence of the model, not experimental support. |
| Mechanistic confidence | Low on its own. A model can place a residue in a buried, exposed or interface-like environment; it cannot show that a substitution changes stability, binding or activity. |
| Phenotype implication | None. Supports reasoning about missense and in-frame changes only in combination with experimental structures, domain annotation and variant-level assays. |
| Absence means | A protein or isoform without a model, or a low-confidence region, has no usable structural context here; this is not evidence of disorder or of the lack of a structure. |
| Overlaps (declared here) | independent_related · `pdbe_sifts` (predicted models versus experimentally determined structures: compare, never merge) |

### `alphamissense` · AlphaMissense

| | |
|---|---|
| Evidence kinds | computational_prediction |
| Allele resolution | allele_specific |
| Mechanism evidence | `protein_function` · prediction · snv, protein_level · protein_structure_function |
| Biological context | Predicts a pathogenicity class for single amino-acid substitutions per transcript or isoform. No tissue, cell-type, dosage, inheritance or mechanism context is carried. |
| Evidence strength | Computational prediction. Class thresholds were chosen to reach 90% precision on ClinVar labels although the model was not trained on ClinVar variants, so agreement with ClinVar is calibration and not independent confirmation. No per-variant uncertainty is provided. |
| Mechanistic confidence | Low. The score does not say whether a change causes loss or gain of function, destabilisation or a dominant-negative effect, and it cannot see effects that act through splicing or regulation. |
| Phenotype implication | A pathogenicity likelihood for the substitution; no phenotype, penetrance or mode of inheritance. |
| Absence means | No score exists for synonymous, nonsense, splice, indel, structural or non-coding variants. Absence of a prediction is not evidence of a benign effect. |

### `biogrid_interactions` · BioGRID protein and genetic interactions (release 5.0.262, human file)

| | |
|---|---|
| Evidence kinds | curated_expert_classification, perturbation_experimental |
| Allele resolution | gene_level |
| Mechanism evidence | `genetic_interaction` · direct · gene_level · pathways_interactions_epistasis (curated genetic interactions from double-perturbation experiments; whole-gene, not allele-specific)<br>`genetic_interaction` · context · gene_level, protein_level · pathways_interactions_epistasis (physical interaction network)<br>`protein_function` · curated · protein_level · pathways_interactions_epistasis (protein-protein interaction partners) |
| Components | physical interactions (affinity capture, two-hybrid and related systems) — curated_expert_classification · gene_level · from native<br>genetic interactions (synthetic lethality, dosage rescue, suppression and related systems) — curated_expert_classification, perturbation_experimental · gene_level · from native |
| Biological context | Interactions from published experiments in many systems (cell lines, yeast and others); each row records the experimental system, throughput and publication. The partner proteins are whole proteins, so an interaction is not tied to a particular isoform, domain or variant. |
| Evidence strength | Curated from the literature; the supporting publications and experimental systems are the evidence, and high-throughput studies are flagged. Counts of supporting publications are an indication of reproducibility, not a probability. |
| Mechanistic confidence | Low to moderate. A genetic interaction shows that perturbing two genes together has an effect that differs from the single perturbations, which is the basis for epistasis and compensation hypotheses; it does not identify the allele, the pathway step or the direction of compensation. A physical interaction does not show that a variant disrupts it. |
| Phenotype implication | Cell-level or organism-level phenotype of the double perturbation, as described in the supporting study; no human disease phenotype. |
| Absence means | A missing interaction means it was not reported in curated publications, which reflects what has been tested; it does not show that two genes do not interact, nor that perturbing both has no combined effect. |
| Overlaps (declared here) | independent_related · `opentargets` (Open Targets interactions come from IntAct, Signor, Reactome and STRING; BioGRID is a separate curation, so the same pair can appear in both and should be counted once per publication) |

### `biogrid_orcs` · BioGRID ORCS (Open Repository of CRISPR Screens)

| | |
|---|---|
| Evidence kinds | perturbation_experimental |
| Allele resolution | gene_level |
| Mechanism evidence | `gene_dosage` · direct · gene_level · perturbation_evidence (whole-gene perturbation in the screened system; not allele-specific) |
| Biological context | Curated results of CRISPR screens in cell lines and cell types of several model organisms. Every hit is conditional on the screen's cell line, library, perturbation type and phenotype readout. |
| Evidence strength | Curated from publications. Hit calls and scores follow each study's own method and threshold and are not uniformly reprocessed. |
| Mechanistic confidence | Moderate for the effect of perturbing the whole gene in the tested context. No particular allele is described, and compensation by paralogs or pathways is not resolved. |
| Phenotype implication | Gene-level cellular phenotype in vitro (for example viability or drug response); not an organismal or human phenotype. |
| Absence means | A gene that is not a hit was not scored as a hit in the curated screens; power, context and threshold limit this. It is not shown to be dispensable elsewhere, and its alleles are not shown to be tolerated. |
| Overlaps (declared here) | independent_related · `depmap` (both report cancer cell-line perturbation screens; keep both and compare, do not count a gene twice when the same screen appears in both) |

### `biotools` · bio.tools and OpenEBench tool registry

| | |
|---|---|
| Evidence kinds | catalogue_metadata |
| Allele resolution | none |
| Mechanism evidence | none (no mechanism content) |
| Biological context | Not biological evidence: the records describe software and services and their EDAM annotations. |
| Evidence strength | Not applicable. Annotation depth and accuracy differ between entries. |
| Mechanistic confidence | Not applicable; the registry makes no biological assertion. |
| Phenotype implication | None: registry entries carry no phenotype information. |
| Absence means | A tool missing from the registry is unregistered, not unavailable; absent annotation terms mean undeclared, not unsupported. |
| Overlaps (declared here) | independent_related · `fairsharing` (both describe resources; use only for discovery) |

### `cadd` · CADD: combined annotation dependent depletion

| | |
|---|---|
| Evidence kinds | computational_prediction |
| Allele resolution | allele_specific |
| Mechanism evidence | `protein_function` · prediction · snv, indel · annotation_integration (combined deleteriousness score) |
| Biological context | A genome-wide score that combines conservation, regulatory and protein-level annotations; it carries no tissue or cell-type context. |
| Evidence strength | Computational prediction with no per-variant uncertainty; the score is a rank among possible substitutions, not a probability. |
| Mechanistic confidence | Low. It does not distinguish protein, splicing and regulatory mechanisms and is not calibrated to a particular phenotype. |
| Phenotype implication | A deleteriousness rank; no phenotype. |
| Absence means | A variant without a score was not precomputed (for example an indel) or lies outside the scored set; absence is not a benign prediction. |
| Overlaps (declared here) | independent_related · `gnomad_variants` (gnomAD INFO carries CADD scores for observed variants, so the score and the frequency come from one file for those variants)<br>independent_related · `dbnsfp` (dbNSFP carries CADD among its predictors) |

### `clingen` · ClinGen (dosage sensitivity, gene-disease validity, region curation)

| | |
|---|---|
| Evidence kinds | curated_expert_classification |
| Allele resolution | gene_level |
| Mechanism evidence | `gene_dosage` · curated · gene_level, region_level, cnv · natural_knockouts_dosage (haploinsufficiency and triplosensitivity scores for whole genes and regions)<br>`structural_rearrangement` · curated · region_level, cnv · structural_variation_cnv (recurrent CNV regions)<br>`outcome_association` · curated · gene_level · clinical_interpretation (gene-disease validity) |
| Biological context | Expert-panel curation of whole genes and genomic regions: dosage sensitivity (loss or gain of an entire gene or region), recurrent CNV regions and gene-disease validity. Gene- and region-level assertions; recurrent regions are intervals. |
| Evidence strength | Graded by curators from published cases. Dosage scores encode evidence levels; 30 stands for a gene associated with an autosomal recessive phenotype and 40 for dosage sensitivity unlikely. A curated judgement, not a measurement, and not allele-specific. |
| Mechanistic confidence | Moderate for the dosage mechanism of whole-gene or whole-region loss or gain. It says nothing about the mechanism of a particular small variant, and contiguous-gene or position effects of multi-gene CNVs are not separated. |
| Phenotype implication | Curated association of dosage change (and, in the validity records, of gene disruption) with a named phenotype; penetrance and expressivity vary. |
| Absence means | A gene or region without a record has not been curated, which differs from a score meaning no evidence (0) or dosage sensitivity unlikely (40). Absence is not evidence of dosage tolerance. |

### `clinvar` · ClinVar

| | |
|---|---|
| Evidence kinds | curated_expert_classification |
| Allele resolution | allele_specific |
| Mechanism evidence | `outcome_association` · curated · snv, indel, sv, cnv · clinical_interpretation |
| Biological context | Submitter-asserted interpretations of variants for named conditions. Each submission carries its own review status, classification, condition and evidence summary; germline and other classifications are separate fields. |
| Evidence strength | Aggregate classification with review status and conflicting-interpretation flags. These grade the assertion, not the underlying evidence; submissions often rest on the same upstream evidence, so agreement between submitters is not independent confirmation. |
| Mechanistic confidence | None as structured data. A classification asserts clinical significance and does not establish mechanism; any functional evidence is in free text. |
| Phenotype implication | Condition-specific clinical significance (pathogenic to benign, or uncertain significance) as named by the submitter. |
| Absence means | A variant missing from ClinVar has not been submitted, which reflects ascertainment (tested genes, populations, laboratories) and not benignity. Uncertain significance means unresolved, not benign. |
| Overlaps (declared here) | may_conflict · `uniprot` (humsavar categories are curated independently from the literature; disagreement with ClinVar is information) |

### `cpic` · CPIC (Clinical Pharmacogenetics Implementation Consortium)

| | |
|---|---|
| Evidence kinds | curated_expert_classification |
| Allele resolution | haplotype_level |
| Mechanism evidence | `pharmacogenomic_haplotype` · curated · snv, indel, cnv, gene_level · functional_haplotypes_pgx (star-allele function and diplotype-phenotype translation)<br>`mitochondrial` · curated · mtdna · functional_haplotypes_pgx (MT-RNR1 alleles and their drug-response guideline) |
| Biological context | Gene-drug pairs with guidelines. Function is assigned to named alleles (star alleles; for X-linked, mitochondrial and HLA genes single alleles) and diplotypes are translated to phenotype categories; categories are specific to the gene and to the drug in the guideline. Allele function records carry the strength of evidence, citations and population frequencies. |
| Evidence strength | Expert-curated allele function and guideline recommendations with CPIC's own evidence levels; assertions derived from published studies, not measurements of a given sample's alleles. |
| Mechanistic confidence | Moderate at the haplotype level for well-studied alleles: function belongs to the haplotype and is not shown separately for each variant. Low for variants outside the allele definition tables, structural alleles defined without breakpoints, and unresolved phase. |
| Phenotype implication | Predicted metaboliser or transporter phenotype and a prescribing recommendation for a diplotype; drug-specific, not general. |
| Absence means | An allele missing from the definition tables has unassigned function, which is different from normal function. A sample with no defining variant detected defaults to the reference allele by convention, an assumption about untested positions. |
| Overlaps (declared here) | independent_related · `pharmvar` (both describe star-allele definitions of the same genes; compare definitions, never sum counts) |

### `dbnsfp` · dbNSFP

| | |
|---|---|
| Evidence kinds | computational_prediction, population_observation, aggregated_integration |
| Allele resolution | allele_specific |
| Mechanism evidence | `protein_function` · prediction · snv · protein_structure_function<br>`splicing` · prediction · snv · annotation_integration<br>`population_constraint` · population · snv · population_frequency |
| Components | deleteriousness predictors (37 algorithms) — computational_prediction · allele_specific · from external: third-party algorithms<br>allele frequencies from 1000 Genomes, gnomAD and others — population_observation · allele_specific · from gnomad_variants |
| Biological context | Per-variant annotations for non-synonymous and splice-site SNVs. Predictors ignore tissue and cell type; they may share training data, so agreement among them is not independent confirmation. |
| Evidence strength | Mixed: predictions (hypotheses), conservation and population frequencies (indirect). Score scales and thresholds differ between algorithms and are not calibrated to each other. |
| Mechanistic confidence | Low. Predictors report a damaging or tolerated call and not a mechanism, and cover only missense and splice-site SNVs. |
| Phenotype implication | Likelihood of a deleterious effect on the protein or on splicing; no phenotype. |
| Absence means | Indels, synonymous and non-coding SNVs, copy-number and structural variants are outside its scope. A missing score for a covered SNV means that algorithm did not score it, not that it is benign. |
| Overlaps (declared here) | aggregates · `alphamissense`<br>aggregates · `gnomad_variants`<br>aggregates · `1000genomes` |

### `dbvar` · dbVar (NCBI database of genomic structural variation)

| | |
|---|---|
| Evidence kinds | population_observation, curated_expert_classification, aggregated_integration |
| Allele resolution | region_level |
| Mechanism evidence | `structural_rearrangement` · population · sv, cnv · structural_variation_cnv<br>`outcome_association` · curated · sv, cnv · structural_variation_cnv (clinical assertions where submitted) |
| Biological context | Archive of SV and CNV calls submitted by studies of differing design, each with the submitting study's method, samples and assertion. Calls are intervals with breakpoint uncertainty. |
| Evidence strength | An archive, not a curated truth set: every submission keeps its own method and quality, calls on the same samples can appear in more than one study, and clinical assertions are submitter judgements. |
| Mechanistic confidence | Low. Records give interval, type and, where submitted, an assertion; gene content, breakpoints and regulatory effects must be annotated separately. |
| Phenotype implication | Clinical assertions, where submitted, link a call to a condition; population studies link none. |
| Absence means | Absence of an overlapping call reflects what submitted studies could detect and report, not rarity or safety. |
| Overlaps (declared here) | identical_records · `dgv` (studies shown in DGV are archived either in dbVar (accessions nsv, nssv) or in DGVa (esv, essv); count each study once) |

### `decipher` · DECIPHER

| | |
|---|---|
| Evidence kinds | case_observation |
| Allele resolution | region_level |
| Mechanism evidence | `outcome_association` · case · cnv, snv, indel · clinical_interpretation (phenotype-linked case data; not assessed beyond the resource description) |
| Biological context | Phenotype-linked CNV and variant data from clinical genomics projects, with public and consented tiers of access. Not assessed beyond the resource description because the site pages could not be retrieved (access denied). |
| Evidence strength | Not assessed. Case-level observations from clinical cohorts are ascertained by clinical referral and are not a population sample. |
| Mechanistic confidence | Not assessed; case data do not establish mechanism. |
| Phenotype implication | Links individual observations to patient phenotypes; interpretation needs phenotype ontology terms, access and consent conditions that have not been verified here. |
| Absence means | Absence reflects consent, referral and data-sharing choices and is not evidence that a CNV is benign or rare. |

### `depmap` · DepMap (Cancer Dependency Map)

| | |
|---|---|
| Evidence kinds | perturbation_experimental |
| Allele resolution | gene_level |
| Mechanism evidence | `gene_dosage` · direct · gene_level · perturbation_evidence (whole-gene knockout or knockdown in cancer cell lines; not allele-specific) |
| Biological context | Dependency screens in cancer cell lines (CRISPR knockout, as documented by Open Targets). Dependency is conditional on lineage, genotype and culture conditions. Other screen types and the context data are unverified because the portal could not be read. |
| Evidence strength | Measured effect of perturbing a gene on cell fitness in cultured cancer lines. The scoring model is unverified (portal unreadable). |
| Mechanistic confidence | Moderate for whole-gene loss of function in the tested lines. It says nothing about a particular allele, other tissues or germline dosage, and does not separate a direct dependency from buffering or compensation. |
| Phenotype implication | Cell-fitness effect in vitro; not an organismal phenotype. Essentiality in cancer cells does not imply a human germline disease phenotype. |
| Absence means | Non-essential in the screened lines means no fitness effect was detected there. It does not show the gene is dispensable in other cell types, conditions or in development. |

### `dgv` · DGV (Database of Genomic Variants)

| | |
|---|---|
| Evidence kinds | population_observation |
| Allele resolution | region_level |
| Mechanism evidence | `structural_rearrangement` · population · cnv, sv · structural_variation_cnv |
| Biological context | Variants reported by studies of healthy-control cohorts. Each study used its own platform and calling method, so breakpoints and frequencies are study-specific. |
| Evidence strength | Observational calls with per-study sample sizes. Frequencies are not harmonised across studies, and overlapping calls from different studies can describe one event. |
| Mechanistic confidence | None. DGV records intervals and types seen in controls and does not describe mechanisms. |
| Phenotype implication | Observation in control cohorts is evidence that a CNV or SV is compatible with absence of overt disease in those individuals, weakened by incomplete phenotyping, reduced penetrance and late onset. |
| Absence means | Not catalogued means not reported in these studies; it does not make a CNV rare, de novo or pathogenic. |
| Overlaps (declared here) | independent_related · `gnomad_sv_cnv` (independent population SV catalogues of different design; compare frequencies, do not pool) |

### `encode` · ENCODE and SCREEN candidate cis-regulatory elements

| | |
|---|---|
| Evidence kinds | experimental_annotation, aggregated_integration |
| Allele resolution | region_level |
| Mechanism evidence | `transcription_regulation` · context · region_level · regulation_expression_splicing (biochemical activity of regions; no allele-specific effect) |
| Biological context | Biochemical assays (chromatin accessibility, transcription-factor binding, histone modification, RNA expression) in biosamples, and candidate cis-regulatory elements derived from them. Activity is specific to the biosample. |
| Evidence strength | Direct measurements of biochemical activity of regions in a biosample. Candidate element classes are computed from those assays and are annotations, not demonstrations of regulatory function. |
| Mechanistic confidence | Low to moderate for the statement that a region is biochemically active in a biosample. No allele-specific effect, target gene or direction of effect is given, and activity does not prove that the element is required for regulation. |
| Phenotype implication | None directly; places variants in regions that are candidates for regulatory effects in particular cell types. |
| Absence means | A variant outside every annotated element may still act in untested cell types, states or at long range. Absence of a peak reflects the assayed biosamples and thresholds. |

### `ensembl_genomes` · Ensembl genome registry, assemblies and gene annotation (human)

| | |
|---|---|
| Evidence kinds | reference_annotation |
| Allele resolution | region_level |
| Mechanism evidence | `splicing` · context · gene_level, region_level · reference_annotation (exon and intron structure of transcript models)<br>`truncation_nmd` · context · gene_level, region_level · reference_annotation (reading frame and exon junctions needed to evaluate decay)<br>`paralog_complex_locus` · context · complex_locus · reference_annotation (alternate, patch and haplotype sequences) |
| Biological context | Reference assemblies (primary, alternate and patch sequences, T2T-CHM13 and haplotype assemblies) and gene annotation. The GRCh38 reference is a single consensus sequence and its annotation is a set of transcript models, not a map of which isoform is expressed where. |
| Evidence strength | Reference data. The properties to verify are coordinates and sequence (done here by published checksums and reference-allele concordance), not an effect or association strength. |
| Mechanistic confidence | Not applicable. It supplies the sequence and transcript frame needed to define consequences (exon boundaries, reading frame, splice sites) and asserts no effect. |
| Phenotype implication | None: reference sequence and gene models carry no phenotype information. |
| Absence means | A position outside annotated genes is not shown to be non-functional. Sequence that is carried only on alternate or haplotype assemblies may be unmappable on the primary reference. |
| Overlaps (declared here) | derived_from · `gencode_mane` (the human gene annotation distributed by Ensembl is merged with the GENCODE gene set) |

### `ensembl_vep` · Ensembl VEP and release-116 plugin data

| | |
|---|---|
| Evidence kinds | computational_prediction, direct_experimental, population_observation, aggregated_integration |
| Allele resolution | allele_specific |
| Mechanism evidence | `protein_function` · prediction · snv, indel · annotation_integration (consequence terms and predictor plugins)<br>`truncation_nmd` · prediction · snv, indel · annotation_integration (rule-based consequence terms (for example frameshift, stop gained) against transcript models)<br>`splicing` · prediction · snv, indel · regulation_expression_splicing (SpliceAI and SpliceVault plugins)<br>`transcription_regulation` · prediction · snv · regulation_expression_splicing (Enformer plugin and promoter annotation)<br>`post_transcriptional` · context · snv, indel · regulation_expression_splicing (RiboSeqORFs plugin)<br>`protein_function` · direct · snv · direct_experimental_variant_effects (MaveDB plugin carries measured scores)<br>`population_constraint` · population · gene_level · natural_knockouts_dosage (LOEUF plugin)<br>`structural_rearrangement` · prediction · sv · annotation_integration (consequence terms for structural variants) |
| Components | consequence terms against transcript models — computational_prediction · allele_specific · from native<br>predictor plugins (for example SpliceAI, SpliceVault, Enformer, AlphaMissense) — computational_prediction · allele_specific · from external: published predictors<br>MaveDB plugin — direct_experimental · allele_specific · from mavedb<br>LOEUF plugin — population_observation · gene_level · from gnomad_constraint |
| Biological context | A tool that annotates variants against transcript models and attaches plugin data. Every annotation is computed or looked up per transcript, so one variant can carry different consequences across transcripts and tools. |
| Evidence strength | Mixed, and inherited from the plugin: consequence terms are deterministic rules against annotation, predictor plugins are hypotheses, MaveDB values are measurements and LOEUF is a gene-level population estimate. The evidence kind of each part persists. |
| Mechanistic confidence | Low to moderate depending on the plugin. A consequence term names a molecular consequence class and does not show loss versus gain of function or tissue context; agreement between predictor plugins is not independent confirmation because they can share training data. |
| Phenotype implication | None directly. |
| Absence means | A missing plugin value means that plugin has no data for the variant (type, transcript or build); for example precomputed SpliceAI files cover SNVs in MANE Select transcripts and no indels. Absence is not a negative prediction. |
| Overlaps (declared here) | aggregates · `alphamissense` (AlphaMissense plugin)<br>aggregates · `mavedb` (MaveDB plugin)<br>aggregates · `gnomad_constraint` (LOEUF plugin) |

### `eqtl_catalogue` · eQTL Catalogue

| | |
|---|---|
| Evidence kinds | statistical_association |
| Allele resolution | allele_specific |
| Mechanism evidence | `transcription_regulation` · association · snv, indel · regulation_expression_splicing (gene-expression QTL)<br>`splicing` · association · snv, indel · regulation_expression_splicing (splicing QTL)<br>`transcription_regulation` · association · snv, indel · association_finemapping_coloc (SuSiE fine-mapped credible sets) |
| Biological context | Public human studies, each with its own tissue or cell type, condition, quantification method and sample size, processed uniformly. Effects are changes in an expression or splicing phenotype per alternate allele in that context. |
| Evidence strength | Statistical association with summary statistics and SuSiE fine-mapped credible sets. Credible sets narrow the candidate variants but do not establish the causal variant. |
| Mechanistic confidence | Low to moderate for the regulatory axis. A QTL says a variant, or a linked variant, is associated with expression or splicing of a gene in that context; colocalisation with a trait and functional assays are needed to connect it to a phenotype. |
| Phenotype implication | A molecular phenotype (gene expression or splice-event usage) in the studied tissue or cell type; not a disease phenotype. |
| Absence means | No QTL means no association was detected in the sampled contexts at the available power. Tissues, cell states, developmental stages and rare alleles that were not sampled are unaddressed. |

### `fairsharing` · FAIRsharing

| | |
|---|---|
| Evidence kinds | catalogue_metadata |
| Allele resolution | none |
| Mechanism evidence | none (no mechanism content) |
| Biological context | Not biological evidence: the records describe standards, databases and policies. |
| Evidence strength | Not applicable. |
| Mechanistic confidence | Not applicable; the registry makes no biological assertion. |
| Phenotype implication | None: registry entries carry no phenotype information. |
| Absence means | A resource missing from the registry is unregistered, not unavailable. |

### `gencode_mane` · GENCODE and MANE transcript sets

| | |
|---|---|
| Evidence kinds | reference_annotation |
| Allele resolution | region_level |
| Mechanism evidence | `truncation_nmd` · context · gene_level, region_level · reference_annotation (transcript structure on which decay rules are evaluated)<br>`splicing` · context · gene_level, region_level · reference_annotation (exon and splice-junction coordinates)<br>`protein_function` · context · gene_level · reference_annotation (the transcript that defines protein-level variant nomenclature) |
| Biological context | Transcript models. The GENCODE comprehensive set lists supported transcripts; MANE Select designates one matched Ensembl and RefSeq transcript per gene as a clinical reference, and MANE Plus Clinical adds transcripts required to report variants. |
| Evidence strength | Curated and automated annotation of transcript structure. Choosing a MANE transcript is a convention and not evidence of the dominant isoform in any tissue. |
| Mechanistic confidence | Not applicable. It defines the frame (exons, coding sequence, start and stop) in which consequences are computed, and the consequence depends on the transcript chosen. |
| Phenotype implication | None directly. A consequence in a non-MANE isoform can be misread as benign or damaging, so consequences are reported per transcript. |
| Absence means | A gene without a MANE Select transcript, or a variant outside MANE transcripts, has no clinical-reference consequence here; this reflects the transcript set and not benignity. |

### `giab_stratifications` · GIAB genome stratifications v3.6 (GRCh38): difficult and complex regions

| | |
|---|---|
| Evidence kinds | reference_annotation |
| Allele resolution | region_level |
| Mechanism evidence | `paralog_complex_locus` · context · complex_locus, region_level · complex_loci_callability (segmental duplications, falsely duplicated and collapsed regions, MHC, KIR and VDJ loci)<br>`repeat_expansion` · context · region_level · complex_loci_callability (tandem-repeat, satellite and homopolymer regions where small-variant calls are unreliable) |
| Biological context | Regions defined from reference-sequence features and read mappability, not from the biology of any sample: segmental duplications and paralogous regions, highly polymorphic immune loci (MHC, KIR, VDJ), falsely duplicated or collapsed reference regions, tandem repeats and reference gaps. |
| Evidence strength | Region annotations from the GIAB and GA4GH benchmarking effort. Membership of a region is a property of the reference and of short-read mappability, not a measurement of any variant call. |
| Mechanistic confidence | Not applicable to mechanism. A region flag says a call there is less reliable or that a paralog may be involved; it cannot separate a true variant from a mis-mapped paralogous sequence variant. |
| Phenotype implication | None: region membership carries no phenotype information. |
| Absence means | A position outside every listed region is not guaranteed to be reliably callable: the sets stratify benchmark performance and do not certify individual calls. A call inside a region must not be discarded or trusted on that basis alone. |

### `gnomad_constraint` · gnomAD constraint metrics (LoF/missense) and pext

| | |
|---|---|
| Evidence kinds | population_observation |
| Allele resolution | gene_level |
| Mechanism evidence | `population_constraint` · population · gene_level · natural_knockouts_dosage<br>`truncation_nmd` · population · gene_level · natural_knockouts_dosage (depletion of predicted loss-of-function variants)<br>`protein_function` · population · gene_level · natural_knockouts_dosage (depletion of missense variants) |
| Biological context | Gene- and transcript-level observed/expected counts of predicted loss-of-function, missense and synonymous variants in the gnomAD cohorts. Reflects selection against heterozygous variants, not against recessive-only effects. |
| Evidence strength | Statistical depletion estimates with confidence bounds. They depend on the expectation model (mutation and methylation rates), coverage and sample size, and are weak for small genes. |
| Mechanistic confidence | Low for any individual variant. Constraint is a gene property and does not say that a specific allele is damaging; strong constraint supports haploinsufficiency or severe missense sensitivity as a hypothesis. |
| Phenotype implication | Intolerance to loss of function suggests, but does not show, a severe early-onset dominant phenotype. Tolerance does not exclude recessive disease, gain of function or missense mechanisms. |
| Absence means | A gene without metrics (too few expected variants or poor coverage) is unscored, not unconstrained. A tolerant score is not evidence that a particular allele is benign. |
| Overlaps (declared here) | derived_from · `gnomad_variants` (constraint is computed from the gnomAD callsets; the counts are the same observations) |

### `gnomad_ld` · gnomAD v2.1.1 linkage disequilibrium matrices

| | |
|---|---|
| Evidence kinds | population_observation |
| Allele resolution | allele_specific |
| Mechanism evidence | `outcome_association` · population · snv · association_finemapping_coloc (LD separates a causal variant from linked variants) |
| Biological context | Correlation of common alleles across haplotypes in the gnomAD genome callset, by genetic ancestry group. |
| Evidence strength | Computed from genotype calls of a large sample; precise for common variants. |
| Mechanistic confidence | None for mechanism; LD identifies linked variants. |
| Phenotype implication | None: linkage disequilibrium carries no phenotype information. |
| Absence means | A variant outside the common set or the callset has no LD entry. |
| Overlaps (declared here) | derived_from · `gnomad_variants` (LD is computed from the gnomAD callset; the same individuals) |

### `gnomad_sv_cnv` · gnomAD structural variants and copy-number variants

| | |
|---|---|
| Evidence kinds | population_observation, computational_prediction |
| Allele resolution | region_level |
| Mechanism evidence | `structural_rearrangement` · population · sv, cnv · structural_variation_cnv<br>`gene_dosage` · prediction · sv, cnv · structural_variation_cnv (PREDICTED_LOF, PREDICTED_COPY_GAIN and related gene-overlap annotations) |
| Biological context | Germline SV calls (genome callset) and CNV calls (exome callset) in the gnomAD cohorts. Algorithm and evidence fields record how each call was made, and PREDICTED_* INFO fields annotate the gene-level consequence of each call. |
| Evidence strength | Population allele counts and frequencies for calls made by combining several algorithms. Sensitivity differs by size and type, and PREDICTED_* consequences are annotation predictions and not assays. |
| Mechanistic confidence | Low to moderate. Dosage-type consequences (whole-gene deletion, duplication, exon overlap) follow from gene overlap; effects of breakpoints on regulation, gene fusion and position effects are not assessed. |
| Phenotype implication | Frequency in population controls bounds pathogenic claims for dosage-sensitive genes; it does not give a phenotype. |
| Absence means | Absence from the callset reflects detectability (size, type, repeat context, coverage) as well as rarity; poorly mappable regions are under-called. |

### `gnomad_variants` · gnomAD v4.1 (exomes and genomes) allele frequencies

| | |
|---|---|
| Evidence kinds | population_observation, computational_prediction |
| Allele resolution | allele_specific |
| Mechanism evidence | `population_constraint` · population · snv, indel · population_frequency<br>`protein_function` · prediction · snv · annotation_integration (in-silico scores carried in INFO)<br>`splicing` · prediction · snv · annotation_integration (SpliceAI and Pangolin scores carried in INFO) |
| Components | allele counts and frequencies per ancestry group — population_observation · allele_specific · from native<br>in-silico scores carried in INFO (CADD, REVEL, SpliceAI, Pangolin, PhyloP, SIFT, PolyPhen) — computational_prediction · allele_specific · from external: published predictors |
| Biological context | Aggregated germline variant calls from exome and genome sequencing cohorts, with allele counts per genetic ancestry group. No phenotype data accompany the variants. |
| Evidence strength | Observed allele counts after quality filtering. Frequencies are reliable only where coverage and callability are adequate (FILTER and coverage fields) and differ by ancestry group; the in-silico scores are predictions and not measurements. |
| Mechanistic confidence | None for mechanism. Frequency is indirect evidence: a high frequency argues against high-penetrance early-onset dominant pathogenicity, while a low frequency is consistent with pathogenicity without supporting it. |
| Phenotype implication | Population frequency bounds the plausible penetrance and prevalence of an allele; no phenotype is attached. |
| Absence means | A variant that is not found was not observed among the sequenced individuals, whose coverage and ancestry representation are uneven. It is not evidence of rarity or harm. |

### `gtex` · GTEx (Genotype-Tissue Expression)

| | |
|---|---|
| Evidence kinds | statistical_association, experimental_annotation |
| Allele resolution | allele_specific |
| Mechanism evidence | `transcription_regulation` · association · snv, indel · association_finemapping_coloc (eQTL)<br>`splicing` · association · snv, indel · regulation_expression_splicing (sQTL)<br>`transcription_regulation` · context · gene_level · regulation_expression_splicing (tissue-specific gene expression) |
| Biological context | Bulk RNA-seq and whole-genome-sequencing genotypes from donors, across non-diseased tissues. Expression and QTL effects are averages over the cell types present in each tissue, and the number of donors differs between tissues. |
| Evidence strength | Association statistics per tissue; release v11 lists 946 subjects (943 with eQTL data) in 54 tissues, 50 of them with eQTL results. Fine-mapping and colocalisation are needed before causal claims. |
| Mechanistic confidence | Low to moderate. An eQTL or sQTL associates a variant, or a linked variant, with expression or splicing in a tissue; allele-specific expression collections add allele-resolved support but remain observational. |
| Phenotype implication | A molecular phenotype in non-diseased tissue; disease relevance needs colocalisation with a trait and the relevant cell type or state. |
| Absence means | No eQTL in a tissue means none was detected at the available power in bulk tissue. Cell-type-specific, context-specific and rare-allele effects are not covered. |

### `gwas_catalog` · NHGRI-EBI GWAS Catalog

| | |
|---|---|
| Evidence kinds | statistical_association |
| Allele resolution | allele_specific |
| Mechanism evidence | `outcome_association` · association · snv, indel · association_finemapping_coloc |
| Biological context | Published genome-wide association results with ontology-annotated traits, study and ancestry tables. Each association belongs to a study with its own cohort, ancestry, genotyping platform and imputation. |
| Evidence strength | Reported association statistics (p-value, effect size and risk-allele frequency where reported). Significance thresholds apply per study, and lead variants tag a locus. |
| Mechanistic confidence | None. An associated locus does not identify the causal variant, gene or mechanism; mapped genes are positional. |
| Phenotype implication | Statistical association with a trait or disease in the studied ancestry. Effects of common variants are generally small, and transfer across ancestries needs evidence. |
| Absence means | A variant with no association in the catalogue has not reached the reporting threshold in catalogued studies, whose trait coverage, ancestry and power vary. It is not evidence of no effect. |

### `hgsvc3` · HGSVC3 release 1.0: assembly-based structural variants from 65 phased genomes (GRCh38)

| | |
|---|---|
| Evidence kinds | population_observation |
| Allele resolution | region_level |
| Mechanism evidence | `structural_rearrangement` · population · sv · structural_variation_cnv (insertions, deletions and inversions in 65 assembled genomes)<br>`paralog_complex_locus` · context · sv, complex_locus · complex_loci_callability (REF_SD and REF_TRF flag calls in segmental duplications and tandem repeats) |
| Biological context | Germline structural variants of 65 individuals of diverse ancestry, called by comparing phased long-read-based assemblies with the reference; the sample set is a research panel and carries no phenotype data. |
| Evidence strength | Observational allele calls with per-haplotype genotypes from assembly alignment. Breakpoints and sizes are sequence-based, which makes them more exact than short-read estimates; sensitivity depends on how well the assemblies align. |
| Mechanistic confidence | None for mechanism. It supplies which structural alleles exist and where; effects on genes, regulation or dosage must come from annotation and from the curated dosage and clinical sources. |
| Phenotype implication | Presence in a research panel without phenotype data; no disease implication. |
| Absence means | A structural variant not listed was not seen in these 65 genomes or lay where assemblies did not align; this is not evidence that the variant is rare in populations or harmless. |
| Overlaps (declared here) | independent_related · `gnomad_sv_cnv` (independent callsets of structural variants (long-read assemblies versus short-read population sequencing); compare by type and reciprocal overlap, never pool counts)<br>independent_related · `1000genomes` (the release also genotypes its variants in 1000 Genomes samples (Genotyping_1kGP directory), so the same individuals can appear in 1000 Genomes-derived callsets; count an individual once) |

### `hpa` · Human Protein Atlas v25.1: tissue and cell-type RNA expression and protein subcellular location

| | |
|---|---|
| Evidence kinds | experimental_annotation, aggregated_integration |
| Allele resolution | gene_level |
| Mechanism evidence | `transcription_regulation` · context · gene_level · regulation_expression_splicing (tissue and cell-type resolved RNA expression of the gene)<br>`protein_function` · context · gene_level · protein_structure_function (subcellular location of the protein) |
| Biological context | RNA sequencing of tissues and single-cell types, normalised per gene, and protein localisation by immunofluorescence imaging of cultured cell lines (the download page lists the image-based location data). Expression is a property of the sampled tissue or cell type and its donors. |
| Evidence strength | Measured expression and localisation summarised per gene; reliability categories grade the localisation evidence. The consensus RNA table combines datasets (see unverified fields), so it is not independent of GTEx. |
| Mechanistic confidence | Low. It shows where a gene is expressed or its protein found, which bears on whether a tissue-specific regulatory or loss-of-function mechanism is plausible; it says nothing about the effect of any allele. |
| Phenotype implication | None directly; tissue specificity informs which tissues could manifest an effect. |
| Absence means | A gene without a value in a tissue or cell type is not shown to be silent there: detection limits, dilution in bulk tissue and tissues or states that were not sampled apply. A gene absent from the localisation table was not localised (13,603 genes), not shown to be non-localised. |
| Overlaps (declared here) | aggregates · `gtex` (the consensus tissue RNA table is described on the download page beside separate HPA and GTEx RNA tables; do not count GTEx expression twice (the combination method was not read)) |

### `impc` · IMPC (International Mouse Phenotyping Consortium)

| | |
|---|---|
| Evidence kinds | perturbation_experimental |
| Allele resolution | gene_level |
| Mechanism evidence | `gene_dosage` · direct · gene_level · perturbation_evidence (mouse knockout phenotype with allele, zygosity and sex; orthology to the human gene is an inference)<br>`gene_dosage` · direct · gene_level · natural_knockouts_dosage (viability of the knockout (lethal, subviable, viable) as the mouse counterpart of loss-of-function intolerance) |
| Biological context | Knockout mouse lines phenotyped in a standardised pipeline. The genotype is an engineered allele (named by allele accession, zygosity and sex) in one organism and strain background; phenotype calls are made per procedure and parameter. For human use the result is gene-level, through orthology. |
| Evidence strength | Direct measurement of phenotypes in knockout mice, with statistical phenotype calls and release-specific thresholds. |
| Mechanistic confidence | Moderate for the gene's role in the mouse. It is not allele-specific, differences between mouse and human physiology limit transfer, and compensation by paralogs is not resolved. |
| Phenotype implication | Mouse phenotype terms. They support a human loss-of-function phenotype only when concordant with human genetics. |
| Absence means | No assertion means none was significant in the procedures performed (screen breadth, sex and genotype coverage vary); a line that was not phenotyped, or was lethal before adult procedures, has no adult phenotype. Absence is not shown to be dispensable. |

### `intact` · IntAct molecular interaction database

| | |
|---|---|
| Evidence kinds | curated_expert_classification |
| Allele resolution | gene_level |
| Mechanism evidence | `genetic_interaction` · context · gene_level, protein_level · pathways_interactions_epistasis (physical interaction network)<br>`protein_function` · curated · protein_level · pathways_interactions_epistasis (protein interaction partners) |
| Biological context | Interactions curated from publications with the detection method of each; they are properties of whole proteins in the studied system. |
| Evidence strength | Literature curation with an interaction confidence score; independent detections raise the score. |
| Mechanistic confidence | Low for variant effects: an interaction says two proteins associate, not that a particular change disrupts it. |
| Phenotype implication | None directly. |
| Absence means | A pair that is not listed was not reported in curated publications; it is not shown not to interact. |
| Overlaps (declared here) | identical_records · `opentargets` (Open Targets carries IntAct interactions as its intact source; use one route and count each pair once)<br>independent_related · `biogrid_interactions` (separate curation of physical interactions; two thirds of the curated Open Targets pairs (IntAct, Signor, Reactome) also appear in BioGRID) |

### `interpro` · InterPro

| | |
|---|---|
| Evidence kinds | reference_annotation, aggregated_integration |
| Allele resolution | protein_residue_level |
| Mechanism evidence | `protein_function` · context · protein_level · protein_structure_function (family, domain and site membership) |
| Biological context | Protein family, domain and site annotations derived from sequence signatures of member databases (Pfam, SMART, PROSITE, CDD, PANTHER, CATH-Gene3D and others). |
| Evidence strength | Signature-based computational annotation integrated from member databases. Membership is a sequence-similarity statement and not a measurement of function. |
| Mechanistic confidence | Low. It places a residue in a family, domain or site; whether a substitution there alters function needs assays, structure or conservation. |
| Phenotype implication | None: family, domain and site membership carries no phenotype information. |
| Absence means | A residue outside every signature is unannotated, not functionally unimportant; signature coverage is incomplete for disordered and lineage-specific regions. |
| Overlaps (declared here) | independent_related · `uniprot` (UniProt entries cross-reference InterPro; domain boundaries from curated features and from signatures can differ, so compare) |

### `ipd_imgt_hla` · IPD-IMGT/HLA Database (HLA allele sequences and nomenclature)

| | |
|---|---|
| Evidence kinds | reference_annotation, curated_expert_classification |
| Allele resolution | haplotype_level |
| Mechanism evidence | `paralog_complex_locus` · context · complex_locus · complex_loci_callability (named allele sequences needed to interpret calls in the HLA region) |
| Biological context | Curated sequences of named HLA alleles from submitted samples of many populations; each allele is a full haplotype of the gene, not a single variant, and the nomenclature encodes protein, synonymous and non-coding differences. |
| Evidence strength | Curated sequence records with nomenclature committee review; sequence completeness differs per allele. |
| Mechanistic confidence | Not applicable. It names and sequences alleles; it asserts no functional or clinical effect. |
| Phenotype implication | None: sequences carry no phenotype information. |
| Absence means | An allele that is not listed has not been submitted and named; a sequence that matches no listed allele may be novel or mis-assembled. Absence says nothing about function. |
| Overlaps (declared here) | independent_related · `ipd_kir` (sister database for the KIR genes under the same terms) |

### `ipd_kir` · IPD-KIR Database (killer-cell immunoglobulin-like receptor allele sequences)

| | |
|---|---|
| Evidence kinds | reference_annotation, curated_expert_classification |
| Allele resolution | haplotype_level |
| Mechanism evidence | `paralog_complex_locus` · context · complex_locus · complex_loci_callability (named allele sequences for the KIR genes) |
| Biological context | Curated sequences of named KIR alleles from submitted samples; the genes are highly similar to each other and differ in presence between individuals. |
| Evidence strength | Curated sequence records; completeness differs per allele. |
| Mechanistic confidence | Not applicable. It names and sequences alleles and asserts no functional or clinical effect. |
| Phenotype implication | None: sequences carry no phenotype information. |
| Absence means | An allele that is not listed has not been submitted and named. Absence of a KIR gene in a person is a real polymorphism and not a gap in the database. |

### `ldlink` · LDlink (NCI): linkage disequilibrium tools for 1000 Genomes populations

| | |
|---|---|
| Evidence kinds | population_observation |
| Allele resolution | allele_specific |
| Mechanism evidence | `outcome_association` · population · snv, indel · association_finemapping_coloc (LD separates a causal variant from linked variants) |
| Biological context | Correlation of alleles across haplotypes in reference population samples; the structure differs between ancestry groups. |
| Evidence strength | Computed from phased reference genotypes; reliable for common variants and small for rare ones. |
| Mechanistic confidence | None for mechanism. LD says which variants travel together and so which are candidates for a shared association signal. |
| Phenotype implication | None: linkage disequilibrium carries no phenotype information. |
| Absence means | A variant not in the reference panel has no LD estimate; low LD in one population does not imply low LD in another. |
| Overlaps (declared here) | derived_from · `1000genomes` (LD is computed from 1000 Genomes genotypes; the individuals are the same) |

### `mavedb` · MaveDB (Multiplexed Assays of Variant Effect database)

| | |
|---|---|
| Evidence kinds | direct_experimental |
| Allele resolution | allele_specific |
| Mechanism evidence | `protein_function` · direct · snv, indel, protein_level · direct_experimental_variant_effects (deep mutational scanning and related protein assays)<br>`transcription_regulation` · direct · snv, indel · direct_experimental_variant_effects (MPRA and related reporter assays) |
| Biological context | Each score set measures variants in a defined assay: a target sequence (gene, exon or regulatory element), an expression system, a selection or readout and a library. Scores are meaningful relative to that assay's controls and direction conventions. |
| Evidence strength | Direct measurements with assay-specific replicates and, where deposited, standard errors. Scales are not harmonised across score sets, and a score set may include variants measured with few observations. |
| Mechanistic confidence | Moderate to high for the assayed molecular function in the assayed system (for example growth, binding, abundance or reporter output). Mechanisms that give the same readout are not distinguished, and splicing or regulatory effects are missed unless the assay includes them. |
| Phenotype implication | A molecular or cellular phenotype of the assay. Mapping to a disease phenotype needs calibration against clinical labels, which the repository holds only where authors supplied it. |
| Absence means | A variant without a score was not measured, was filtered for low counts, or lies outside the assayed region; a position with no data is untested and not tolerant. A neutral score means no effect on the assayed function only. |

### `mitomap` · MITOMAP: human mitochondrial genome database

| | |
|---|---|
| Evidence kinds | curated_expert_classification |
| Allele resolution | allele_specific |
| Mechanism evidence | `mitochondrial` · curated · mtdna · clinical_interpretation (curated disease-associated mitochondrial DNA variants (not assessed: the site was not readable)) |
| Biological context | Not assessed: the site was not readable. Mitochondrial variants act through heteroplasmy, tissue distribution and the mitochondrial genetic code. |
| Evidence strength | Not assessed. |
| Mechanistic confidence | Not assessed. |
| Phenotype implication | Not assessed. |
| Absence means | Not assessed; a mitochondrial variant missing from every catalogued source is simply uncovered. |

### `mpravardb` · MPRAVarDB: massively parallel reporter assay variant database

| | |
|---|---|
| Evidence kinds | direct_experimental, aggregated_integration |
| Allele resolution | allele_specific |
| Mechanism evidence | `transcription_regulation` · direct · snv, indel · direct_experimental_variant_effects (MPRA allelic regulatory activity)<br>`transcription_regulation` · direct · snv · regulation_expression_splicing (variants at GWAS loci and eQTLs) |
| Components | published MPRA measurements, collected — direct_experimental · allele_specific · from external: the original MPRA studies |
| Biological context | Reporter-assay measurements of variant sequences in cultured cell lines; the regulatory effect is that of the sequence in an episomal reporter, outside its native chromatin context. |
| Evidence strength | Direct allele-specific measurements compiled from studies; the quality control differs between studies. Details of the stored measure were not read. |
| Mechanistic confidence | Moderate for regulatory activity of the tested sequence in the tested cell line; it does not show the effect at the native locus, in other cell types or on a target gene. |
| Phenotype implication | A molecular phenotype (reporter activity); the disease label names the GWAS trait of the tested locus and is not a measured phenotype. |
| Absence means | A variant that was not tested has no data; an untested cell type or state is unknown, and a variant with no effect in one cell line may act in another. |
| Overlaps (declared here) | independent_related · `mavedb` (MaveDB holds MPRA score sets; the same experiment can appear in both and should be counted once) |

### `nhgri_education` · NHGRI educational pages (Talking Glossary, 1000 Genomes Project page)

| | |
|---|---|
| Evidence kinds | documentation |
| Allele resolution | none |
| Mechanism evidence | none (no mechanism content) |
| Biological context | Not biological evidence: educational definitions and narration, and an archived project description. |
| Evidence strength | Not applicable; educational material with no variant-level data. |
| Mechanistic confidence | Not applicable. |
| Phenotype implication | None: educational pages carry no phenotype data. |
| Absence means | A concept not covered is not defined here; the pages carry no data from which absence could be inferred. |

### `opentargets` · Open Targets Platform

| | |
|---|---|
| Evidence kinds | aggregated_integration, statistical_association, curated_expert_classification, perturbation_experimental, computational_prediction |
| Allele resolution | gene_level |
| Mechanism evidence | `outcome_association` · association · snv, indel · association_finemapping_coloc (credible sets, studies and colocalisation)<br>`outcome_association` · prediction · snv, indel · association_finemapping_coloc (locus-to-gene predictions)<br>`outcome_association` · curated · snv, indel, gene_level · clinical_interpretation (evidence from ClinVar, ClinGen and UniProt variants)<br>`gene_dosage` · direct · gene_level · perturbation_evidence (essentiality, CRISPR screens and knockout mouse phenotypes)<br>`transcription_regulation` · prediction · region_level · regulation_expression_splicing (enhancer-to-gene links)<br>`pharmacogenomic_haplotype` · curated · snv, gene_level · functional_haplotypes_pgx (pharmacogenomics dataset)<br>`genetic_interaction` · context · gene_level, protein_level · pathways_interactions_epistasis (interaction and pathway context)<br>`transcription_regulation` · association · snv, indel · regulation_expression_splicing (molecular-QTL credible sets: eQTL, single-cell eQTL and transcript-usage QTL)<br>`splicing` · association · snv, indel · regulation_expression_splicing (sQTL credible sets)<br>`protein_function` · association · snv, indel · association_finemapping_coloc (pQTL credible sets (protein abundance)) |
| Components | credible_set, study and colocalisation (GWAS and molecular QTL) — statistical_association · allele_specific · from external: published GWAS and QTL studies<br>l2g_prediction and enhancer_to_gene — computational_prediction · allele_specific · from native<br>evidence_eva — curated_expert_classification · allele_specific · from clinvar<br>evidence_clingen — curated_expert_classification · gene_level · from clingen<br>evidence_uniprot_variants — curated_expert_classification · allele_specific · from uniprot<br>target_essentiality — perturbation_experimental · gene_level · from depmap<br>evidence_impc — perturbation_experimental · gene_level · from impc |
| Biological context | Integrated target-disease evidence. The target-disease association is gene-level; variant-level datasets (credible sets, colocalisation, evidence from ClinVar and UniProt) are allele-specific but remain those of the originating source, whose profile governs. |
| Evidence strength | Aggregated, with the evidence kind of each part retained: association statistics, curated assertions, perturbation results and predictions. Overall association scores combine datasets with fixed weights and are ranking aids, not evidence strength. |
| Mechanistic confidence | Inherited from the originating dataset and generally low for allele-level mechanism: gene-level scores say nothing about a specific allele, and a locus-to-gene prediction is a hypothesis about the target gene. |
| Phenotype implication | Disease and phenotype associations of genes and loci, ontology-mapped; each part carries the phenotype semantics of its source. |
| Absence means | A target without evidence for a disease has no integrated evidence in this release, which reflects what the integrated sources cover and not the absence of a relationship. |
| Overlaps (declared here) | aggregates · `alphamissense`<br>aggregates · `clingen`<br>aggregates · `clinvar` (ClinVar via EVA)<br>aggregates · `eqtl_catalogue`<br>aggregates · `gtex`<br>aggregates · `gwas_catalog`<br>aggregates · `gnomad_variants` (population allele frequencies)<br>aggregates · `reactome`<br>aggregates · `uniprot`<br>aggregates · `pharmgkb` (pharmacogenetics data from ClinPGx)<br>aggregates · `depmap` (target_essentiality)<br>aggregates · `impc` (evidence_impc) |

### `pdbe_sifts` · PDBe SIFTS: PDB chain to UniProt residue mappings

| | |
|---|---|
| Evidence kinds | reference_annotation, experimental_structure |
| Allele resolution | protein_residue_level |
| Mechanism evidence | `protein_function` · context · protein_level · protein_structure_function (which UniProt residues are covered and observed in an experimental structure) |
| Biological context | Residue-level mapping between PDB chains and UniProtKB sequences. The structures are experimentally determined constructs, often fragments or domains with ligands or engineered changes; the mapping says which UniProt residues each chain covers and which are observed. |
| Evidence strength | Experimental structures with mapping tables. Resolution and model quality belong to each structure entry and are not in the mapping. |
| Mechanistic confidence | Moderate for the structural context of observed residues (position, burial, contacts) when the construct is relevant. It does not show the effect of a substitution, and unobserved residues have no coordinates. |
| Phenotype implication | None: residue mappings carry no phenotype information. |
| Absence means | A residue not covered by any PDB chain, or not observed in the covering chains, has no experimental coordinates; this is not evidence of disorder. |

### `pgs_catalog` · Polygenic Score (PGS) Catalog

| | |
|---|---|
| Evidence kinds | statistical_association, catalogue_metadata |
| Allele resolution | allele_specific |
| Mechanism evidence | `outcome_association` · association · snv, indel · association_finemapping_coloc |
| Biological context | Each score is a weighted sum of effect-allele dosages from a published model, developed for a stated trait in stated development samples and evaluated in other samples; ancestry and cohort come with the score. |
| Evidence strength | Reported performance of the score as a whole in its evaluation samples. A single weight is an association estimate for one allele in a model and carries no individual evidence strength. |
| Mechanistic confidence | None. Weights are statistical effect estimates in a prediction model; they identify neither causal variants nor mechanisms. |
| Phenotype implication | Predisposition to a trait in the ancestry the score was developed and evaluated in. A raw score is meaningful only against a reference distribution of the same score and ancestry. |
| Absence means | A variant that is not a weight of a score was not selected or estimated by that model. Absence from a score is not evidence of no effect, and a missing genotype is not a zero dosage. |
| Overlaps (declared here) | independent_related · `gwas_catalog` (scores are developed from genome-wide association results, which may or may not be among the GWAS Catalog studies; compare, do not pool) |

### `pharmcat` · PharmCAT (Pharmacogenomics Clinical Annotation Tool)

| | |
|---|---|
| Evidence kinds | curated_expert_classification |
| Allele resolution | haplotype_level |
| Mechanism evidence | `pharmacogenomic_haplotype` · curated · snv, indel, cnv, gene_level · functional_haplotypes_pgx (applies curated allele definitions and guideline tables; a method, not a new assertion) |
| Biological context | A method that applies published allele definitions and guidelines to a genotype; the clinical meaning is that of the guidelines, specific to the gene and the drug. |
| Evidence strength | Inherited from the guidelines it applies; a calling pipeline adds the risk of wrong calls at positions it cannot see. |
| Mechanistic confidence | As for the CPIC and PharmVar tables it uses: moderate at the haplotype level, low for variants outside the definitions. |
| Phenotype implication | Predicted metaboliser or risk phenotype and a recommendation for a diplotype, drug-specific. |
| Absence means | A gene, drug or allele outside the covered set receives no report; an uncalled position is treated according to the tool's rules, which can default to the reference allele. |
| Overlaps (declared here) | derived_from · `cpic` (applies the CPIC allele definitions and guideline tables; the evidence is the CPIC evidence, so count it once) |

### `pharmgkb` · ClinPGx / PharmGKB

| | |
|---|---|
| Evidence kinds | curated_expert_classification, aggregated_integration |
| Allele resolution | allele_specific |
| Mechanism evidence | `pharmacogenomic_haplotype` · curated · snv, indel, gene_level · functional_haplotypes_pgx<br>`outcome_association` · curated · snv, indel · clinical_interpretation (variant-drug annotations with levels of evidence) |
| Components | variant and clinical annotations from the literature — curated_expert_classification · allele_specific · from native<br>clinical guideline annotations — curated_expert_classification · haplotype_level · from cpic |
| Biological context | Annotations linking variants, haplotypes, genes and drugs from published studies. Each annotation records the studied variant or haplotype, the drug, a phenotype category, the population and a level of evidence. |
| Evidence strength | Literature-derived annotations with levels of evidence for clinical annotations. Individual studies are small and often conflict, and the annotation level is not an effect size. |
| Mechanistic confidence | Low to moderate. Relationships state an association between genotype and drug response; molecular mechanism is annotated only where the study reports it. |
| Phenotype implication | Drug-specific response phenotypes (efficacy, toxicity, dosage) in studied populations; not transferable across drugs or ancestries without evidence. |
| Absence means | No annotation means no curated study for that variant-drug pair and not that genotype is irrelevant; rare variants and alleles outside the studied panels have no annotation. |
| Overlaps (declared here) | aggregates · `cpic` (ClinPGx comprises PharmGKB, CPIC and PharmCAT; clinical guideline annotations include CPIC guidelines, so count a guideline once) |

### `pharmvar` · PharmVar (Pharmacogene Variation Consortium)

| | |
|---|---|
| Evidence kinds | curated_expert_classification |
| Allele resolution | haplotype_level |
| Mechanism evidence | `pharmacogenomic_haplotype` · context · snv, indel, gene_level · functional_haplotypes_pgx (defining variants of star alleles; function is assigned elsewhere) |
| Biological context | Nomenclature and defining variants for star alleles of 15 pharmacogenes, on RefSeqGene, GRCh37 and GRCh38 coordinates. Structural variation and complex variant types are described separately from sequence haplotypes. |
| Evidence strength | Consortium-curated nomenclature. An allele definition lists the variants that define a haplotype; a name is a definition and not evidence of function. |
| Mechanistic confidence | Not applicable to function. A star allele defines which variants occur together; function and drug response are assigned by other resources (CPIC, PharmGKB) and the mechanism of each variant is not given. |
| Phenotype implication | None directly. |
| Absence means | A variant not in any allele definition is not a named allele, so its function is unassigned and not normal. Genes outside the 15 are not covered. |

### `proteingym` · ProteinGym v1.3: deep mutational scanning benchmark assays (substitutions and indels)

| | |
|---|---|
| Evidence kinds | direct_experimental, aggregated_integration |
| Allele resolution | allele_specific |
| Mechanism evidence | `protein_function` · direct · protein_level, indel · direct_experimental_variant_effects (assay-specific readout such as growth, binding, abundance or proteolytic stability)<br>`genetic_interaction` · direct · protein_level · direct_experimental_variant_effects (assays with multiple mutants allow non-additivity (epistasis) to be estimated) |
| Components | published deep mutational scanning measurements, re-formatted — direct_experimental · allele_specific · from external: the original studies |
| Biological context | Each assay measures a variant library of one protein or domain in a specific system; the organism, selection assay and selection type are given per assay in the reference file. The score reflects the assayed function under that selection. |
| Evidence strength | Direct measurements curated and processed by the benchmark maintainers: scores are as published in assay-specific units and the binarised label uses a per-assay cutoff. Replicate-level uncertainty is not retained. |
| Mechanistic confidence | Moderate for the assayed molecular function of the substitution in the assayed system. Stability, activity and binding effects are not separated unless the selection type says so. |
| Phenotype implication | A molecular or cellular phenotype of the assay. Correspondence to a human disease phenotype needs calibration, and most assays are of non-human proteins or of protein stability. |
| Absence means | A mutation missing from an assay table was not measured or was filtered by the maintainers; an untested position is unknown and not tolerant. A protein with no assay has no curated measurement here. |
| Overlaps (declared here) | independent_related · `mavedb` (some assays are the same experiments as MaveDB score sets; match by publication and target before combining, and count each experiment once) |

### `reactome` · Reactome

| | |
|---|---|
| Evidence kinds | curated_expert_classification |
| Allele resolution | gene_level |
| Mechanism evidence | `genetic_interaction` · context · gene_level, protein_level · pathways_interactions_epistasis (pathway membership and reactions; no allele-specific information) |
| Biological context | Curated human pathways and reactions with gene, protein and small-molecule participants and interaction tables. Pathways are idealised and not tissue-specific. |
| Evidence strength | Expert curation from the literature. Membership states that a participant takes part in a reaction, not how much a change in it matters. |
| Mechanistic confidence | Low for variants. Pathway context helps frame hypotheses on compensation and buffering but does not show the effect of any allele. |
| Phenotype implication | None directly; pathway enrichment is hypothesis-generating. |
| Absence means | A gene absent from the pathways is uncurated, not unrelated; coverage is biased towards well-studied biology. |

### `refseq` · NCBI RefSeq: GRCh38.p14 accession map, alternate-locus placement and gene annotation (release RS_2025_08)

| | |
|---|---|
| Evidence kinds | reference_annotation |
| Allele resolution | region_level |
| Mechanism evidence | `truncation_nmd` · context · gene_level, region_level · reference_annotation (RefSeq transcript structures on which ClinVar HGVS is expressed)<br>`splicing` · context · gene_level, region_level · reference_annotation (exon structure of RefSeq transcripts)<br>`paralog_complex_locus` · context · complex_locus · complex_loci_callability (alternate-locus and patch scaffolds and the regions they cover) |
| Biological context | Reference accessions and transcript models of one haploid reference with alternate loci and patches. The annotation is a set of transcript models and not a map of which isoform is expressed where. |
| Evidence strength | Reference data from an automated and curated annotation pipeline; correctness is checked here by published checksums, structural consistency of the GFF and joins to the other reference files. |
| Mechanistic confidence | Not applicable. It supplies identifiers, roles and transcript structure needed to define consequences and to recognise regions with competing representations; it asserts no effect. |
| Phenotype implication | None: accession maps and gene models carry no phenotype information. |
| Absence means | A sequence not in the report is not part of GRCh38.p14; a position outside annotated features is not shown to be non-functional; a region without alternate loci or patches is not shown to be free of duplication. |
| Overlaps (declared here) | independent_related · `gencode_mane` (MANE Select pairs one RefSeq with one Ensembl transcript per gene; the RefSeq half of MANE should reappear in this annotation, and the two gene sets are otherwise independent annotations of the same genome)<br>independent_related · `ensembl_genomes` (independent annotation of the same assembly; compare transcript sets and boundaries, never merge counts) |

### `strchive` · STRchive and related tandem-repeat catalogues

| | |
|---|---|
| Evidence kinds | curated_expert_classification |
| Allele resolution | region_level |
| Mechanism evidence | `repeat_expansion` · curated · str · repeat_expansions_str<br>`outcome_association` · curated · str · repeat_expansions_str (locus-disease association with inheritance mode) |
| Biological context | Curated tandem-repeat loci associated with Mendelian disease: locus coordinates (hg38, hg19, T2T-CHM13), repeat motif and locus structure, inheritance mode, and normal and pathogenic length thresholds. |
| Evidence strength | Expert curation from the literature. Thresholds differ between loci and are summarised as single numbers; intermediate ranges, repeat interruptions and somatic instability are not captured by one threshold. |
| Mechanistic confidence | Moderate for repeat length and motif as the disease-relevant variable at curated loci. The molecular mechanism (RNA toxicity, polyglutamine proteotoxicity, silencing) varies by locus and is not resolved by length alone. |
| Phenotype implication | Length above the pathogenic threshold associates with the named disease; onset, penetrance and severity depend on length and modifiers. |
| Absence means | A locus missing from the catalogue is not curated as disease-associated. It is not shown to be benign, and the catalogue covers only loci with published disease associations. |

### `string_db` · STRING protein-protein association networks

| | |
|---|---|
| Evidence kinds | curated_expert_classification, computational_prediction |
| Allele resolution | gene_level |
| Mechanism evidence | `genetic_interaction` · context · gene_level, protein_level · pathways_interactions_epistasis (functional association network) |
| Biological context | Associations from experiments, curated pathway databases, text mining, co-expression and genomic context, combined into one score; not tissue-specific. |
| Evidence strength | A mixture of curated, experimental and predicted evidence combined into a confidence; the channels differ greatly in reliability. |
| Mechanistic confidence | Low. It places proteins in a network; it does not show that a variant affects a link. |
| Phenotype implication | None directly. |
| Absence means | A pair without a link has no combined evidence above the stored threshold; it is not shown to be unrelated. |
| Overlaps (declared here) | identical_records · `opentargets` (Open Targets carries STRING associations as its string source; use one route and count each pair once) |

### `trexplorer` · TRExplorer: genome-wide tandem-repeat catalogue

| | |
|---|---|
| Evidence kinds | reference_annotation, population_observation, aggregated_integration |
| Allele resolution | region_level |
| Mechanism evidence | `repeat_expansion` · population · str · repeat_expansions_str (polymorphism of repeat loci in long-read genomes)<br>`repeat_expansion` · context · str, region_level · repeat_expansions_str (locus and motif definitions) |
| Biological context | Repeat loci defined from the reference sequence with motifs, annotated with polymorphism in assembled genomes of diverse ancestry (HPRC). |
| Evidence strength | A curated and computed catalogue; polymorphism summaries come from long-read genomes of a research panel. |
| Mechanistic confidence | None for mechanism. A locus being polymorphic or near a gene does not make an allele pathogenic; thresholds for disease come from curated sources such as STRchive. |
| Phenotype implication | None directly. |
| Absence means | A repeat not in the catalogue was not defined as a locus (for example a very long or impure repeat); this is not evidence that it is stable. |
| Overlaps (declared here) | independent_related · `strchive` (STRchive curates disease loci with thresholds; TRExplorer extends the locus list genome-wide) |

### `uniprot` · UniProtKB (Swiss-Prot and TrEMBL)

| | |
|---|---|
| Evidence kinds | curated_expert_classification, reference_annotation |
| Allele resolution | protein_residue_level |
| Mechanism evidence | `protein_function` · curated · protein_level, gene_level · protein_structure_function (function, domain, site and modification annotations)<br>`outcome_association` · curated · protein_level · clinical_interpretation (humsavar categories for missense variants) |
| Components | sequence, function, domain, site and modification annotation — curated_expert_classification, reference_annotation · protein_residue_level · from native<br>humsavar variant index (ACMG/AMP-style categories) — curated_expert_classification · allele_specific · from native |
| Biological context | Reviewed (Swiss-Prot) human protein entries with literature-curated annotation. Annotations refer to the canonical sequence, and isoforms are described by their differences. |
| Evidence strength | Curated with evidence codes and literature references per annotation. The humsavar categories are curated from the literature with ACMG/AMP-style classes (likely pathogenic or pathogenic, likely benign or benign, uncertain significance). |
| Mechanistic confidence | Moderate for annotated functional sites and domains; variant categories are classifications and not mechanisms. |
| Phenotype implication | Disease associations are listed per variant where curated; a variant with several diseases appears once per disease. |
| Absence means | A residue without annotation is uncurated, not functionally unimportant. A variant missing from humsavar is not curated as disease-related or benign. |

## Overlaps and deduplication

Policy: identical records are combined once; independent experiments are never merged only because they concern the same gene or variant; assertions that can disagree are kept side by side, and disagreement is information. For an aggregate, reuse the aggregate and consult the original for fields the aggregate drops; count the underlying evidence once. A pair is declared once, by the source that depends on or restates the other.

| Source | Relation | With | What is shared or contested |
|---|---|---|---|
| `dbnsfp` | aggregates | `1000genomes` |  |
| `dbnsfp` | aggregates | `alphamissense` |  |
| `dbnsfp` | aggregates | `gnomad_variants` |  |
| `ensembl_vep` | aggregates | `alphamissense` | AlphaMissense plugin |
| `ensembl_vep` | aggregates | `gnomad_constraint` | LOEUF plugin |
| `ensembl_vep` | aggregates | `mavedb` | MaveDB plugin |
| `hpa` | aggregates | `gtex` | the consensus tissue RNA table is described on the download page beside separate HPA and GTEx RNA tables; do not count GTEx expression twice (the combination method was not read) |
| `opentargets` | aggregates | `alphamissense` |  |
| `opentargets` | aggregates | `clingen` |  |
| `opentargets` | aggregates | `clinvar` | ClinVar via EVA |
| `opentargets` | aggregates | `depmap` | target_essentiality |
| `opentargets` | aggregates | `eqtl_catalogue` |  |
| `opentargets` | aggregates | `gnomad_variants` | population allele frequencies |
| `opentargets` | aggregates | `gtex` |  |
| `opentargets` | aggregates | `gwas_catalog` |  |
| `opentargets` | aggregates | `impc` | evidence_impc |
| `opentargets` | aggregates | `pharmgkb` | pharmacogenetics data from ClinPGx |
| `opentargets` | aggregates | `reactome` |  |
| `opentargets` | aggregates | `uniprot` |  |
| `pharmgkb` | aggregates | `cpic` | ClinPGx comprises PharmGKB, CPIC and PharmCAT; clinical guideline annotations include CPIC guidelines, so count a guideline once |
| `ensembl_genomes` | derived_from | `gencode_mane` | the human gene annotation distributed by Ensembl is merged with the GENCODE gene set |
| `gnomad_constraint` | derived_from | `gnomad_variants` | constraint is computed from the gnomAD callsets; the counts are the same observations |
| `gnomad_ld` | derived_from | `gnomad_variants` | LD is computed from the gnomAD callset; the same individuals |
| `ldlink` | derived_from | `1000genomes` | LD is computed from 1000 Genomes genotypes; the individuals are the same |
| `pharmcat` | derived_from | `cpic` | applies the CPIC allele definitions and guideline tables; the evidence is the CPIC evidence, so count it once |
| `dbvar` | identical_records | `dgv` | studies shown in DGV are archived either in dbVar (accessions nsv, nssv) or in DGVa (esv, essv); count each study once |
| `intact` | identical_records | `opentargets` | Open Targets carries IntAct interactions as its intact source; use one route and count each pair once |
| `string_db` | identical_records | `opentargets` | Open Targets carries STRING associations as its string source; use one route and count each pair once |
| `alliance_genome` | independent_related | `impc` | mouse knockout phenotypes appear in both through MGI; count a mouse line once |
| `alphafold_db` | independent_related | `pdbe_sifts` | predicted models versus experimentally determined structures: compare, never merge |
| `biogrid_interactions` | independent_related | `opentargets` | Open Targets interactions come from IntAct, Signor, Reactome and STRING; BioGRID is a separate curation, so the same pair can appear in both and should be counted once per publication |
| `biogrid_orcs` | independent_related | `depmap` | both report cancer cell-line perturbation screens; keep both and compare, do not count a gene twice when the same screen appears in both |
| `biotools` | independent_related | `fairsharing` | both describe resources; use only for discovery |
| `cadd` | independent_related | `dbnsfp` | dbNSFP carries CADD among its predictors |
| `cadd` | independent_related | `gnomad_variants` | gnomAD INFO carries CADD scores for observed variants, so the score and the frequency come from one file for those variants |
| `cpic` | independent_related | `pharmvar` | both describe star-allele definitions of the same genes; compare definitions, never sum counts |
| `dgv` | independent_related | `gnomad_sv_cnv` | independent population SV catalogues of different design; compare frequencies, do not pool |
| `hgsvc3` | independent_related | `1000genomes` | the release also genotypes its variants in 1000 Genomes samples (Genotyping_1kGP directory), so the same individuals can appear in 1000 Genomes-derived callsets; count an individual once |
| `hgsvc3` | independent_related | `gnomad_sv_cnv` | independent callsets of structural variants (long-read assemblies versus short-read population sequencing); compare by type and reciprocal overlap, never pool counts |
| `intact` | independent_related | `biogrid_interactions` | separate curation of physical interactions; two thirds of the curated Open Targets pairs (IntAct, Signor, Reactome) also appear in BioGRID |
| `interpro` | independent_related | `uniprot` | UniProt entries cross-reference InterPro; domain boundaries from curated features and from signatures can differ, so compare |
| `ipd_imgt_hla` | independent_related | `ipd_kir` | sister database for the KIR genes under the same terms |
| `mpravardb` | independent_related | `mavedb` | MaveDB holds MPRA score sets; the same experiment can appear in both and should be counted once |
| `pgs_catalog` | independent_related | `gwas_catalog` | scores are developed from genome-wide association results, which may or may not be among the GWAS Catalog studies; compare, do not pool |
| `proteingym` | independent_related | `mavedb` | some assays are the same experiments as MaveDB score sets; match by publication and target before combining, and count each experiment once |
| `refseq` | independent_related | `ensembl_genomes` | independent annotation of the same assembly; compare transcript sets and boundaries, never merge counts |
| `refseq` | independent_related | `gencode_mane` | MANE Select pairs one RefSeq with one Ensembl transcript per gene; the RefSeq half of MANE should reappear in this annotation, and the two gene sets are otherwise independent annotations of the same genome |
| `trexplorer` | independent_related | `strchive` | STRchive curates disease loci with thresholds; TRExplorer extends the locus list genome-wide |
| `clinvar` | may_conflict | `uniprot` | humsavar categories are curated independently from the literature; disagreement with ClinVar is information |

Relations: **identical_records** — Both sources carry the same underlying records; combine them once, not twice; **derived_from** — This source is built from, or re-expresses, the other; treat as one lineage; **aggregates** — This source aggregates the other among others; prefer the aggregate for reuse and the original for fields the aggregate drops; **independent_related** — Independent measurements or curations of related phenomena; keep both and compare; **may_conflict** — Independent assertions that can disagree; never merged, disagreement is information

## Components of mixed sources

Where one resource carries parts of different evidence kinds, each part keeps its own kind and resolution. The origin tells whether the part is native, taken from another catalogued source, or from an external method.

| Source | Component | Evidence kind | Allele resolution | Origin |
|---|---|---|---|---|
| `biogrid_interactions` | physical interactions (affinity capture, two-hybrid and related systems) | curated_expert_classification | gene_level | native |
| `biogrid_interactions` | genetic interactions (synthetic lethality, dosage rescue, suppression and related systems) | curated_expert_classification, perturbation_experimental | gene_level | native |
| `dbnsfp` | deleteriousness predictors (37 algorithms) | computational_prediction | allele_specific | external: third-party algorithms |
| `dbnsfp` | allele frequencies from 1000 Genomes, gnomAD and others | population_observation | allele_specific | gnomad_variants |
| `ensembl_vep` | consequence terms against transcript models | computational_prediction | allele_specific | native |
| `ensembl_vep` | predictor plugins (for example SpliceAI, SpliceVault, Enformer, AlphaMissense) | computational_prediction | allele_specific | external: published predictors |
| `ensembl_vep` | MaveDB plugin | direct_experimental | allele_specific | mavedb |
| `ensembl_vep` | LOEUF plugin | population_observation | gene_level | gnomad_constraint |
| `gnomad_variants` | allele counts and frequencies per ancestry group | population_observation | allele_specific | native |
| `gnomad_variants` | in-silico scores carried in INFO (CADD, REVEL, SpliceAI, Pangolin, PhyloP, SIFT, PolyPhen) | computational_prediction | allele_specific | external: published predictors |
| `mpravardb` | published MPRA measurements, collected | direct_experimental | allele_specific | external: the original MPRA studies |
| `opentargets` | credible_set, study and colocalisation (GWAS and molecular QTL) | statistical_association | allele_specific | external: published GWAS and QTL studies |
| `opentargets` | l2g_prediction and enhancer_to_gene | computational_prediction | allele_specific | native |
| `opentargets` | evidence_eva | curated_expert_classification | allele_specific | clinvar |
| `opentargets` | evidence_clingen | curated_expert_classification | gene_level | clingen |
| `opentargets` | evidence_uniprot_variants | curated_expert_classification | allele_specific | uniprot |
| `opentargets` | target_essentiality | perturbation_experimental | gene_level | depmap |
| `opentargets` | evidence_impc | perturbation_experimental | gene_level | impc |
| `pharmgkb` | variant and clinical annotations from the literature | curated_expert_classification | allele_specific | native |
| `pharmgkb` | clinical guideline annotations | curated_expert_classification | haplotype_level | cpic |
| `proteingym` | published deep mutational scanning measurements, re-formatted | direct_experimental | allele_specific | external: the original studies |
| `uniprot` | sequence, function, domain, site and modification annotation | curated_expert_classification, reference_annotation | protein_residue_level | native |
| `uniprot` | humsavar variant index (ACMG/AMP-style categories) | curated_expert_classification | allele_specific | native |
