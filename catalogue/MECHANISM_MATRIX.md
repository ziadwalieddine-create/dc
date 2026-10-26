# Mechanism matrix: variant types by mechanism axis

*Generated 2026-10-26 02:55 UTC by `catalogue/build_catalogue.py` from `catalogue/sources/*.json`, `acquisition/targets.json`, `acquisition/manifest.jsonl`, `acquisition/validation_report.json` and `acquisition/compat_report.json`. Edit the sources, not this file.*

Which mechanisms the catalogued evidence can inform for each type of called variant, on what basis, and what is usable now. The axes and their competing explanations are defined in `catalogue/vocab.json`; a hypothesis is never recorded as a finding, and a gene-level result is never recorded as an allele-specific one.

Cell = number of **usable** sources by basis: **D** direct measurement, **C** curated assertion, **A** statistical association, **K** case-level observation, **P** prediction, **F** population frequency or constraint, **X** context only; `(+n)` = n further sources declare the entry but are not yet usable (not acquired, validated or mapped, or without a verified licence); `·` = nothing; `n/a` = the axis does not apply to that class. Gene-, region- and protein-level entries are counted in the last column only and never fill an allele-level cell.

| Mechanism axis | snv | indel | cnv | sv | str | complex_locus | mtdna | gene / region / protein level |
|---|---|---|---|---|---|---|---|---|
| Clinical or trait outcome (mechanism-agnostic) (`outcome_association`) | A2 C2 P1 (+4) | A2 C2 P1 (+3) | C1 (+2) | C1 (+1) | C1 | · | · | C2 (+1) |
| Protein function, stability and interactions (`protein_function`) | A1 D2 P1 (+4) | A1 D2 (+2) | n/a | n/a | n/a | n/a | n/a | C2 D2 F1 X5 (+2) |
| Truncation, frameshift and nonsense-mediated decay (`truncation_nmd`) | · (+1) | · (+1) | n/a | n/a | n/a | n/a | n/a | F1 X3 |
| Splicing (`splicing`) | A2 P1 (+3) | A2 (+2) | n/a | n/a | n/a | n/a | n/a | X3 |
| Transcriptional regulation (`transcription_regulation`) | A1 D1 (+6) | A1 D1 (+4) | n/a | · | n/a | n/a | n/a | P1 X3 |
| UTR, RNA stability and translation (`post_transcriptional`) | · (+1) | · (+1) | n/a | n/a | n/a | n/a | n/a | · |
| Gene dosage (`gene_dosage`) | n/a | · | C1 P1 | P1 | n/a | n/a | n/a | C1 D4 (+2) |
| Structural rearrangement (`structural_rearrangement`) | n/a | n/a | F2 (+2) | F3 (+3) | n/a | · | n/a | · (+1) |
| Repeat expansion and instability (`repeat_expansion`) | n/a | n/a | n/a | n/a | C1 (+2) | n/a | n/a | X1 (+1) |
| Pharmacogenomic haplotypes (`pharmacogenomic_haplotype`) | C3 X1 (+1) | C2 X1 (+1) | C1 (+1) | · | n/a | · | n/a | C3 X1 (+1) |
| Paralogs and complex loci (`paralog_complex_locus`) | · | · | · | X1 | n/a | X4 (+2) | n/a | X1 |
| Population frequency and constraint (`population_constraint`) | F1 (+2) | F1 (+1) | · | · | · | n/a | · | F1 (+1) |
| Genetic interaction, epistasis and compensation (`genetic_interaction`) | · | · | · | n/a | n/a | n/a | n/a | D2 X3 (+2) |
| Mitochondrial DNA variation (`mitochondrial`) | n/a | n/a | n/a | n/a | n/a | n/a | C1 (+1) | n/a |

## Computed gaps

An applicable (axis, class) cell is a gap when no usable source offers direct, curated, associational or case-level evidence for it. Context, population and prediction entries can inform but do not fill it.

| Mechanism axis | Class | Usable but non-empirical | Declared, not yet usable (first missing stage) |
|---|---|---|---|
| Clinical or trait outcome (mechanism-agnostic) | complex_locus | — | — |
| Clinical or trait outcome (mechanism-agnostic) | mtdna | — | — |
| Truncation, frameshift and nonsense-mediated decay | snv | — | `ensembl_vep` (P; acquired) |
| Truncation, frameshift and nonsense-mediated decay | indel | — | `ensembl_vep` (P; acquired) |
| Transcriptional regulation | sv | — | — |
| UTR, RNA stability and translation | snv | — | `ensembl_vep` (X; acquired) |
| UTR, RNA stability and translation | indel | — | `ensembl_vep` (X; acquired) |
| Gene dosage | sv | `gnomad_sv_cnv` (P) | — |
| Gene dosage | indel | — | — |
| Structural rearrangement | sv | `dgv` (F), `gnomad_sv_cnv` (F), `hgsvc3` (F) | `1000genomes` (F; acquired), `dbvar` (F; acquired), `ensembl_vep` (P; acquired) |
| Structural rearrangement | cnv | `dgv` (F), `gnomad_sv_cnv` (F) | `clingen` (C; mapped), `dbvar` (F; acquired) |
| Structural rearrangement | complex_locus | — | — |
| Pharmacogenomic haplotypes | sv | — | — |
| Pharmacogenomic haplotypes | complex_locus | — | — |
| Paralogs and complex loci | snv | — | — |
| Paralogs and complex loci | indel | — | — |
| Paralogs and complex loci | cnv | — | — |
| Paralogs and complex loci | sv | `hgsvc3` (X) | — |
| Paralogs and complex loci | complex_locus | `ensembl_genomes` (X), `giab_stratifications` (X), `hgsvc3` (X), `refseq` (X) | `ipd_imgt_hla` (X; acquired), `ipd_kir` (X; acquired) |
| Population frequency and constraint | snv | `gnomad_variants` (F) | `1000genomes` (F; acquired), `dbnsfp` (F; acquired) |
| Population frequency and constraint | indel | `gnomad_variants` (F) | `1000genomes` (F; acquired) |
| Population frequency and constraint | cnv | — | — |
| Population frequency and constraint | sv | — | — |
| Population frequency and constraint | str | — | — |
| Population frequency and constraint | mtdna | — | — |
| Genetic interaction, epistasis and compensation | snv | — | — |
| Genetic interaction, epistasis and compensation | indel | — | — |
| Genetic interaction, epistasis and compensation | cnv | — | — |

28 applicable cells have no usable empirical evidence.

## Clinical or trait outcome (mechanism-agnostic) (`outcome_association`)

Assertions or associations that link variants or genes to a disease, trait or molecular phenotype without establishing how: clinical significance, gene-disease validity, GWAS and QTL associations.

- **Competing explanations:** causal variant versus a linked variant; ascertainment and testing bias; population-specific effects; reverse causation or confounding in association; classification circularity when predictors were calibrated on the same labels
- **What discriminates them:** fine-mapping and colocalisation, independent replication, functional assays, segregation, review status and conflicting-interpretation flags
- **Applies to:** snv, indel, cnv, sv, str, complex_locus, mtdna; levels: gene_level, region_level

| Source | Basis | Classes | Dimension | Usable now | Note |
|---|---|---|---|---|---|
| `gwas_catalog` | association | snv, indel | association_finemapping_coloc | ✓ |  |
| `opentargets` | association | snv, indel | association_finemapping_coloc | ✓ | credible sets, studies and colocalisation |
| `decipher` | case | cnv, snv, indel | clinical_interpretation | no (acquired) | phenotype-linked case data; not assessed beyond the resource description |
| `alliance_genome` | curated | gene_level | perturbation_evidence | no (acquired) | curated disease-model annotations |
| `clingen` | curated | gene_level | clinical_interpretation | ✓ | gene-disease validity |
| `clinvar` | curated | snv, indel, sv, cnv | clinical_interpretation | ✓ |  |
| `dbvar` | curated | sv, cnv | structural_variation_cnv | no (acquired) | clinical assertions where submitted |
| `opentargets` | curated | snv, indel, gene_level | clinical_interpretation | ✓ | evidence from ClinVar, ClinGen and UniProt variants |
| `pharmgkb` | curated | snv, indel | clinical_interpretation | no (acquired) | variant-drug annotations with levels of evidence |
| `strchive` | curated | str | repeat_expansions_str | ✓ | locus-disease association with inheritance mode |
| `uniprot` | curated | protein_level | clinical_interpretation | ✓ | humsavar categories for missense variants |
| `gnomad_ld` | population | snv | association_finemapping_coloc | no (acquired) | LD separates a causal variant from linked variants |
| `ldlink` | population | snv, indel | association_finemapping_coloc | no (acquired) | LD separates a causal variant from linked variants |
| `opentargets` | prediction | snv, indel | association_finemapping_coloc | ✓ | locus-to-gene predictions |

## Protein function, stability and interactions (`protein_function`)

Altered activity, stability, folding, localisation or interactions from amino-acid substitutions, in-frame indels and start or stop changes.

- **Competing explanations:** loss versus gain of function; destabilisation versus active-site or interface disruption versus mislocalisation; dominant-negative versus haploinsufficient; a coding change acting instead through splicing or translation
- **What discriminates them:** allele-specific functional assays, structure and interface context, domain annotation, expression and dosage data, segregation and clinical classification
- **Applies to:** snv, indel; levels: protein_level, gene_level

| Source | Basis | Classes | Dimension | Usable now | Note |
|---|---|---|---|---|---|
| `opentargets` | association | snv, indel | association_finemapping_coloc | ✓ | pQTL credible sets (protein abundance) |
| `alphafold_db` | context | protein_level | protein_structure_function | ✓ |  |
| `gencode_mane` | context | gene_level | reference_annotation | ✓ | the transcript that defines protein-level variant nomenclature |
| `hpa` | context | gene_level | protein_structure_function | ✓ | subcellular location of the protein |
| `interpro` | context | protein_level | protein_structure_function | ✓ | family, domain and site membership |
| `pdbe_sifts` | context | protein_level | protein_structure_function | ✓ | which UniProt residues are covered and observed in an experimental structure |
| `biogrid_interactions` | curated | protein_level | pathways_interactions_epistasis | ✓ | protein-protein interaction partners |
| `intact` | curated | protein_level | pathways_interactions_epistasis | no (acquired) | protein interaction partners |
| `uniprot` | curated | protein_level, gene_level | protein_structure_function | ✓ | function, domain, site and modification annotations |
| `ensembl_vep` | direct | snv | direct_experimental_variant_effects | ✓ | MaveDB plugin carries measured scores |
| `mavedb` | direct | snv, indel, protein_level | direct_experimental_variant_effects | ✓ | deep mutational scanning and related protein assays |
| `proteingym` | direct | protein_level, indel | direct_experimental_variant_effects | ✓ | assay-specific readout such as growth, binding, abundance or proteolytic stability |
| `gnomad_constraint` | population | gene_level | natural_knockouts_dosage | ✓ | depletion of missense variants |
| `alphamissense` | prediction | snv, protein_level | protein_structure_function | no (acquired) |  |
| `cadd` | prediction | snv, indel | annotation_integration | no (acquired) | combined deleteriousness score |
| `dbnsfp` | prediction | snv | protein_structure_function | no (acquired) |  |
| `ensembl_vep` | prediction | snv, indel | annotation_integration | no (acquired) | consequence terms and predictor plugins |
| `gnomad_variants` | prediction | snv | annotation_integration | ✓ | in-silico scores carried in INFO |

## Truncation, frameshift and nonsense-mediated decay (`truncation_nmd`)

Premature termination, frameshifts, start loss and stop loss: decay of the transcript (haploinsufficiency) versus stable truncated proteins.

- **Competing explanations:** decay versus escape (last exon, 50-nucleotide rule); alternative start codons; exon skipping that restores frame; transcript-specific effects; genes that tolerate loss of function
- **What discriminates them:** transcript structure and decay rules, expression across transcripts (pext), loss-of-function curation, constraint, population loss-of-function counts, clinical classification
- **Applies to:** snv, indel; levels: gene_level

| Source | Basis | Classes | Dimension | Usable now | Note |
|---|---|---|---|---|---|
| `ensembl_genomes` | context | gene_level, region_level | reference_annotation | ✓ | reading frame and exon junctions needed to evaluate decay |
| `gencode_mane` | context | gene_level, region_level | reference_annotation | ✓ | transcript structure on which decay rules are evaluated |
| `refseq` | context | gene_level, region_level | reference_annotation | ✓ | RefSeq transcript structures on which ClinVar HGVS is expressed |
| `gnomad_constraint` | population | gene_level | natural_knockouts_dosage | ✓ | depletion of predicted loss-of-function variants |
| `ensembl_vep` | prediction | snv, indel | annotation_integration | no (acquired) | rule-based consequence terms (for example frameshift, stop gained) against transcript models |

## Splicing (`splicing`)

Disrupted or created splice sites, exonic and intronic regulators, branch points and deep-intronic pseudoexons; exon skipping, inclusion and cryptic sites.

- **Competing explanations:** a splicing effect versus a protein effect of the same exonic change; partial versus complete; tissue-specific; rescue by cryptic sites
- **What discriminates them:** splice predictors, tissue RNA-seq splice events and sQTLs, minigene or MAVE splicing assays, transcript structure
- **Applies to:** snv, indel; levels: gene_level, region_level

| Source | Basis | Classes | Dimension | Usable now | Note |
|---|---|---|---|---|---|
| `eqtl_catalogue` | association | snv, indel | regulation_expression_splicing | no (acquired) | splicing QTL |
| `gtex` | association | snv, indel | regulation_expression_splicing | ✓ | sQTL |
| `opentargets` | association | snv, indel | regulation_expression_splicing | ✓ | sQTL credible sets |
| `ensembl_genomes` | context | gene_level, region_level | reference_annotation | ✓ | exon and intron structure of transcript models |
| `gencode_mane` | context | gene_level, region_level | reference_annotation | ✓ | exon and splice-junction coordinates |
| `refseq` | context | gene_level, region_level | reference_annotation | ✓ | exon structure of RefSeq transcripts |
| `dbnsfp` | prediction | snv | annotation_integration | no (acquired) |  |
| `ensembl_vep` | prediction | snv, indel | regulation_expression_splicing | no (acquired) | SpliceAI and SpliceVault plugins |
| `gnomad_variants` | prediction | snv | annotation_integration | ✓ | SpliceAI and Pangolin scores carried in INFO |

## Transcriptional regulation (`transcription_regulation`)

Promoters, enhancers, silencers, insulators, transcription-factor binding, chromatin accessibility, 3D contacts and expression QTLs.

- **Competing explanations:** causal variant versus linkage with the causal variant; tissue- or state-specific effects; target-gene ambiguity; enhancer redundancy and buffering
- **What discriminates them:** fine-mapping and colocalisation, MPRA and CRISPRi/a, chromatin and cCRE annotation, allelic imbalance, tissue specificity
- **Applies to:** snv, indel, sv; levels: region_level, gene_level

| Source | Basis | Classes | Dimension | Usable now | Note |
|---|---|---|---|---|---|
| `eqtl_catalogue` | association | snv, indel | regulation_expression_splicing | no (acquired) | gene-expression QTL |
| `eqtl_catalogue` | association | snv, indel | association_finemapping_coloc | no (acquired) | SuSiE fine-mapped credible sets |
| `gtex` | association | snv, indel | association_finemapping_coloc | no (acquired) | eQTL |
| `opentargets` | association | snv, indel | regulation_expression_splicing | ✓ | molecular-QTL credible sets: eQTL, single-cell eQTL and transcript-usage QTL |
| `encode` | context | region_level | regulation_expression_splicing | ✓ | biochemical activity of regions; no allele-specific effect |
| `gtex` | context | gene_level | regulation_expression_splicing | ✓ | tissue-specific gene expression |
| `hpa` | context | gene_level | regulation_expression_splicing | ✓ | tissue and cell-type resolved RNA expression of the gene |
| `mavedb` | direct | snv, indel | direct_experimental_variant_effects | ✓ | MPRA and related reporter assays |
| `mpravardb` | direct | snv, indel | direct_experimental_variant_effects | no (acquired) | MPRA allelic regulatory activity |
| `mpravardb` | direct | snv | regulation_expression_splicing | no (acquired) | variants at GWAS loci and eQTLs |
| `ensembl_vep` | prediction | snv | regulation_expression_splicing | no (acquired) | Enformer plugin and promoter annotation |
| `opentargets` | prediction | region_level | regulation_expression_splicing | ✓ | enhancer-to-gene links |

## UTR, RNA stability and translation (`post_transcriptional`)

5 prime and 3 prime UTR variants, miRNA and RNA-binding sites, polyadenylation, RNA structure, upstream open reading frames and translation initiation.

- **Competing explanations:** transcript abundance versus RNA stability versus translation efficiency
- **What discriminates them:** ribosome profiling, UTR reporter assays, conservation, binding-site annotation
- **Applies to:** snv, indel; levels: region_level

| Source | Basis | Classes | Dimension | Usable now | Note |
|---|---|---|---|---|---|
| `ensembl_vep` | context | snv, indel | regulation_expression_splicing | no (acquired) | RiboSeqORFs plugin |

## Gene dosage (`gene_dosage`)

Copy-number loss or gain, or loss of function, causing haploinsufficiency or triplosensitivity; imprinting; compound heterozygosity.

- **Competing explanations:** dosage versus disruption of regulatory elements (position effect); contiguous-gene effects of multi-gene CNVs; mosaicism; incomplete penetrance
- **What discriminates them:** dosage-sensitivity curation, constraint, population CNV frequencies, gene content of the interval
- **Applies to:** cnv, sv, indel; levels: gene_level, region_level

| Source | Basis | Classes | Dimension | Usable now | Note |
|---|---|---|---|---|---|
| `clingen` | curated | gene_level, region_level, cnv | natural_knockouts_dosage | ✓ | haploinsufficiency and triplosensitivity scores for whole genes and regions |
| `alliance_genome` | direct | gene_level | perturbation_evidence | no (acquired) | model-organism knockout and mutant phenotypes, orthology-mapped |
| `biogrid_orcs` | direct | gene_level | perturbation_evidence | ✓ | whole-gene perturbation in the screened system; not allele-specific |
| `depmap` | direct | gene_level | perturbation_evidence | no (acquired) | whole-gene knockout or knockdown in cancer cell lines; not allele-specific |
| `impc` | direct | gene_level | perturbation_evidence | ✓ | mouse knockout phenotype with allele, zygosity and sex; orthology to the human gene is an inference |
| `impc` | direct | gene_level | natural_knockouts_dosage | ✓ | viability of the knockout (lethal, subviable, viable) as the mouse counterpart of loss-of-function intolerance |
| `opentargets` | direct | gene_level | perturbation_evidence | ✓ | essentiality, CRISPR screens and knockout mouse phenotypes |
| `gnomad_sv_cnv` | prediction | sv, cnv | structural_variation_cnv | ✓ | PREDICTED_LOF, PREDICTED_COPY_GAIN and related gene-overlap annotations |

## Structural rearrangement (`structural_rearrangement`)

Inversions, translocations, insertions and complex events: gene disruption, fusion transcripts, disruption of regulatory domains.

- **Competing explanations:** breakpoint gene disruption versus enhancer hijacking versus dosage change versus mobile-element effects
- **What discriminates them:** breakpoint mapping, domain (TAD) annotation, expression, long-read resolution
- **Applies to:** sv, cnv, complex_locus; levels: region_level

| Source | Basis | Classes | Dimension | Usable now | Note |
|---|---|---|---|---|---|
| `clingen` | curated | region_level, cnv | structural_variation_cnv | no (mapped) | recurrent CNV regions |
| `1000genomes` | population | sv | structural_variation_cnv | no (acquired) |  |
| `dbvar` | population | sv, cnv | structural_variation_cnv | no (acquired) |  |
| `dgv` | population | cnv, sv | structural_variation_cnv | ✓ |  |
| `gnomad_sv_cnv` | population | sv, cnv | structural_variation_cnv | ✓ |  |
| `hgsvc3` | population | sv | structural_variation_cnv | ✓ | insertions, deletions and inversions in 65 assembled genomes |
| `ensembl_vep` | prediction | sv | annotation_integration | no (acquired) | consequence terms for structural variants |

## Repeat expansion and instability (`repeat_expansion`)

Short tandem repeat and VNTR length change: toxic RNA, polyglutamine proteotoxicity, silencing, somatic instability.

- **Competing explanations:** length threshold versus repeat interruptions and purity versus somatic mosaicism versus modifiers; motif identity
- **What discriminates them:** locus catalogues with thresholds and motifs, population allele-size distributions, long-read repeat sequence
- **Applies to:** str; levels: region_level

| Source | Basis | Classes | Dimension | Usable now | Note |
|---|---|---|---|---|---|
| `giab_stratifications` | context | region_level | complex_loci_callability | ✓ | tandem-repeat, satellite and homopolymer regions where small-variant calls are unreliable |
| `trexplorer` | context | str, region_level | repeat_expansions_str | no (acquired) | locus and motif definitions |
| `strchive` | curated | str | repeat_expansions_str | ✓ |  |
| `trexplorer` | population | str | repeat_expansions_str | no (acquired) | polymorphism of repeat loci in long-read genomes |

## Pharmacogenomic haplotypes (`pharmacogenomic_haplotype`)

Star alleles, including copy-number and hybrid alleles, that set enzyme or transporter activity and drug response.

- **Competing explanations:** function assigned by haplotype versus by the individual variant; phase ambiguity; rare variants outside the haplotype definitions; structural alleles; drug-specific effects
- **What discriminates them:** allele-definition tables, function and phenotype assignments, guideline evidence, long-read phasing
- **Applies to:** snv, indel, cnv, sv, complex_locus; levels: gene_level

| Source | Basis | Classes | Dimension | Usable now | Note |
|---|---|---|---|---|---|
| `pharmvar` | context | snv, indel, gene_level | functional_haplotypes_pgx | ✓ | defining variants of star alleles; function is assigned elsewhere |
| `cpic` | curated | snv, indel, cnv, gene_level | functional_haplotypes_pgx | ✓ | star-allele function and diplotype-phenotype translation |
| `opentargets` | curated | snv, gene_level | functional_haplotypes_pgx | ✓ | pharmacogenomics dataset |
| `pharmcat` | curated | snv, indel, cnv, gene_level | functional_haplotypes_pgx | no (acquired) | applies curated allele definitions and guideline tables; a method, not a new assertion |
| `pharmgkb` | curated | snv, indel, gene_level | functional_haplotypes_pgx | ✓ |  |

## Paralogs and complex loci (`paralog_complex_locus`)

Variants in segmental duplications and paralog or pseudogene-rich or highly polymorphic loci (HLA, KIR, SMN1 and SMN2, CYP2D6 and CYP2D7, GBA and GBAP1): mis-mapping, gene conversion, copy-number polymorphism.

- **Competing explanations:** a true variant versus a paralogous sequence variant that was mis-mapped; gene conversion versus point mutation; copy-number-dependent genotype; a reference-sequence artefact (collapsed or falsely duplicated region)
- **What discriminates them:** region annotation of duplications and mappability, assembly-aware calling, long reads, haplotype-specific references
- **Applies to:** snv, indel, cnv, sv, complex_locus; levels: region_level

| Source | Basis | Classes | Dimension | Usable now | Note |
|---|---|---|---|---|---|
| `ensembl_genomes` | context | complex_locus | reference_annotation | ✓ | alternate, patch and haplotype sequences |
| `giab_stratifications` | context | complex_locus, region_level | complex_loci_callability | ✓ | segmental duplications, falsely duplicated and collapsed regions, MHC, KIR and VDJ loci |
| `hgsvc3` | context | sv, complex_locus | complex_loci_callability | ✓ | REF_SD and REF_TRF flag calls in segmental duplications and tandem repeats |
| `ipd_imgt_hla` | context | complex_locus | complex_loci_callability | no (acquired) | named allele sequences needed to interpret calls in the HLA region |
| `ipd_kir` | context | complex_locus | complex_loci_callability | no (acquired) | named allele sequences for the KIR genes |
| `refseq` | context | complex_locus | complex_loci_callability | ✓ | alternate-locus and patch scaffolds and the regions they cover |

## Population frequency and constraint (`population_constraint`)

Indirect evidence from allele frequency, homozygosity, constraint and selection.

- **Competing explanations:** tolerance versus pathogenicity inferred from frequency; founder effects; ascertainment bias; absence from a callset is not evidence of harm
- **What discriminates them:** ancestry-specific frequencies, coverage, constraint with confidence intervals, curated loss-of-function variants
- **Applies to:** snv, indel, cnv, sv, str, mtdna; levels: gene_level

| Source | Basis | Classes | Dimension | Usable now | Note |
|---|---|---|---|---|---|
| `1000genomes` | population | snv, indel | population_frequency | no (acquired) |  |
| `dbnsfp` | population | snv | population_frequency | no (acquired) |  |
| `ensembl_vep` | population | gene_level | natural_knockouts_dosage | no (acquired) | LOEUF plugin |
| `gnomad_constraint` | population | gene_level | natural_knockouts_dosage | ✓ |  |
| `gnomad_variants` | population | snv, indel | population_frequency | ✓ |  |

## Genetic interaction, epistasis and compensation (`genetic_interaction`)

Epistasis, pathway buffering, paralog compensation, modifier loci and oligogenic or digenic effects.

- **Competing explanations:** paralog compensation versus pathway redundancy versus modifier alleles; a gene-level knockout result does not describe a particular allele
- **What discriminates them:** combinatorial perturbation screens, interaction and pathway context, modifier association evidence
- **Applies to:** snv, indel, cnv; levels: gene_level, protein_level

| Source | Basis | Classes | Dimension | Usable now | Note |
|---|---|---|---|---|---|
| `biogrid_interactions` | context | gene_level, protein_level | pathways_interactions_epistasis | ✓ | physical interaction network |
| `intact` | context | gene_level, protein_level | pathways_interactions_epistasis | no (acquired) | physical interaction network |
| `opentargets` | context | gene_level, protein_level | pathways_interactions_epistasis | ✓ | interaction and pathway context |
| `reactome` | context | gene_level, protein_level | pathways_interactions_epistasis | ✓ | pathway membership and reactions; no allele-specific information |
| `string_db` | context | gene_level, protein_level | pathways_interactions_epistasis | no (acquired) | functional association network |
| `biogrid_interactions` | direct | gene_level | pathways_interactions_epistasis | ✓ | curated genetic interactions from double-perturbation experiments; whole-gene, not allele-specific |
| `proteingym` | direct | protein_level | direct_experimental_variant_effects | ✓ | assays with multiple mutants allow non-additivity (epistasis) to be estimated |

## Mitochondrial DNA variation (`mitochondrial`)

mtDNA variants: heteroplasmy, tRNA, rRNA and protein-coding effects under the mitochondrial genetic code, maternal inheritance.

- **Competing explanations:** heteroplasmy level and tissue distribution versus nuclear modifiers; threshold effects
- **What discriminates them:** heteroplasmy-aware callsets, mtDNA curation resources
- **Applies to:** mtdna

| Source | Basis | Classes | Dimension | Usable now | Note |
|---|---|---|---|---|---|
| `cpic` | curated | mtdna | functional_haplotypes_pgx | ✓ | MT-RNR1 alleles and their drug-response guideline |
| `mitomap` | curated | mtdna | clinical_interpretation | no (acquired) | curated disease-associated mitochondrial DNA variants (not assessed: the site was not readable) |

## Sources with no mechanism content

`biotools`, `fairsharing`, `nhgri_education`
