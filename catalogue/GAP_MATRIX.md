# Gap matrix

*Generated 2026-10-25 11:02 UTC by `catalogue/build_catalogue.py` from `catalogue/sources/*.json`, `acquisition/targets.json`, `acquisition/manifest.jsonl`, `acquisition/validation_report.json` and `acquisition/compat_report.json`. Edit the sources, not this file.*

For each dimension: the best stage reached by any source, what is actually usable, what is missing and the next action. The first three columns are computed; the gap and next-action text is curated in `catalogue/gaps.json`.

| Dimension | Priority | Best stage reached | Usable evidence now | Main gap | Next action |
|---|---|---|---|---|---|
| D01 Direct experimental variant-effect measurements | highest | usable | `ensembl_vep`, `mavedb` | Genomic placement is incomplete. Ensembl expands protein-level variants into candidate codon alleles (990 of 1,402 sets); MaveDB's own mapper assigns genomic HGVS to only about 317k of 3.07M annotated variants. 197 of our sets are absent from the CC0 dump (non-CC0 or newer). Superseded versions (51 human-target sets) exist only as a labelled extension. Scores are not comparable across sets and calibrations exist for few sets. | Build one variant-level table keyed by VRS digest/ClinGen allele ID that links assay variant, mapped alleles and annotations; treat superseded versions explicitly; do not pool raw scores across sets. |
| D02 Protein structure, stability, activity and interactions | high | usable | `uniprot` | UniProt reviewed human table is acquired and joins to 98% of MaveDB proteins; AlphaFold models, InterPro matches and variant-level stability predictors (Ensembl ProtVar 1.3 GB, mutfunc 2.0 GB) are not fetched. | Fetch AlphaFold models and InterPro matches per MaveDB target protein through their APIs; evaluate ProtVar/mutfunc for variant-level structural effects. |
| D03 Transcription, splicing, regulatory elements, tissue-specific expression | high | usable | `gtex` | Only GTEx median expression and SCREEN cCRE intervals are acquired. No eQTL/sQTL tables, no splicing predictions (SpliceAI 27 GB, SpliceVault 843 MB) and no sequence-to-function predictions (Enformer 2.8 GB). cCREs are predicted elements, not variant effects. | Pull Open Targets credible_set/colocalisation/enhancer_to_gene subsets; evaluate SpliceVault; use MPRA score sets already in MaveDB for regulatory variants. |
| D04 Knockout, knockdown, CRISPRi/a, overexpression, knock-in and rescue experiments | high | usable | `biogrid_orcs` | BioGRID ORCS human screens are acquired and validated (licence unconfirmed). DepMap is behind a human-verification wall and was not fetched. IMPC (mouse) is not fetched. No variant-level perturbation (base-editing) source is catalogued. | User downloads DepMap Public from its downloads page (or use Open Targets target_essentiality); fetch Open Targets evidence_crispr/evidence_impc; assess base-editor screens. |
| D05 Human natural knockouts, constraint and dosage sensitivity | highest | usable | `clingen`, `gnomad_constraint` | gnomAD v4.1.2 constraint and ClinGen HI/TS gene and region sets are acquired. Not covered: gnomAD LoF curation, homozygous LoF counts, Open Targets gene-burden evidence. | Fetch gnomAD lof_curation and Open Targets evidence_gene_burden. |
| D06 Genetic association, fine-mapping, LD and colocalisation | high | usable | `gwas_catalog` | Only the flat GWAS Catalog associations table is acquired. Fine-mapping, colocalisation and locus-to-gene evidence (Open Targets) and eQTL Catalogue credible sets are catalogued but not fetched. No LD reference panel is acquired (1000 Genomes genotypes not fetched). | Pull selected Open Targets Parquet datasets (study, credible_set, colocalisation, l2g_prediction) for the MaveDB gene set; decide the LD source. |
| D07 Functional haplotypes and pharmacogenomics | medium | usable | `cpic`, `pharmgkb` | ClinPGx and CPIC tables are acquired. PharmVar licence and bulk route are unverified, so haplotype definitions at sequence level are missing. | Read PharmVar terms manually; use its API if permitted. |
| D08 Structural variation and copy-number variation | high | usable | `dgv`, `gnomad_sv_cnv` | gnomAD exome CNV and ClinGen recurrent CNV/dosage regions are acquired. Not fetched: gnomAD SV (1.74 GB), DGV (364 MB), 1000 Genomes/HGSVC SV callsets. dbVar currency is unknown and DECIPHER denied access. | Fetch gnomAD SV and DGV GRCh38; inspect dbVar contents and dates. |
| D09 Repeat expansions and short tandem repeats (STR/VNTR) | medium | usable | `strchive` | STRchive disease loci (82, hg38) are acquired. Population allele-size distributions (gnomAD TR) and genome-wide reference catalogues (PlatinumTRs, 7.8M loci) are not. | Fetch gnomAD tandem-repeat and PlatinumTRs data when STR calling is in scope. |
| D10 Pathways, molecular interactions, epistasis and compensation | medium | usable | `reactome` | Reactome mappings are acquired. No protein/genetic interaction database is acquired (BioGRID interactions, IntAct, STRING not assessed) and nothing in the library quantifies epistasis or compensation. | Assess interaction resources and digenic/oligogenic disease databases; Open Targets interaction dataset is a candidate. |
| D11 Clinical variant interpretation and disease association | highest | usable | `clingen`, `clinvar` | ClinVar (2026-10-19) and ClinGen gene-disease validity are acquired and join to MaveDB variants (10.1% of mapped SNV alleles have an exact ClinVar record). OMIM is licensed, HPO/Monarch are not assessed, DECIPHER is access-controlled. | Assess HPO/Monarch/MONDO; keep ClinVar release dates in every result. |
| D12 Population allele frequencies and haplotype structure | highest | accessible | — | No frequency table is acquired directly. gnomAD v4.1 VCFs are hundreds of GiB per chromosome and are deliberately deferred; the dump annotations carry gnomAD allele frequency for only 5,942 mapped variants. 1000 Genomes contributes a 3,202-sample population table only. | Use regional tabix/API access or Ensembl/dbNSFP frequency fields for the variants of interest; avoid mirroring gnomAD. |
| D13 Integrated variant annotation and predictors | high | usable | `alphamissense` | The best integrators are not yet usable: dbNSFP needs institutional registration, Ensembl VEP caches/plugin data are not fetched (MaveDB plugin file is), Open Targets Parquet is not fetched. AlphaMissense gene-level scores are acquired; the 643 MB variant file is tier B. | Ask the user to register for dbNSFP; fetch the AlphaMissense variant file and the VEP GRCh38 cache when annotation begins. |
| D14 Reference genomes, transcripts and gene models | highest | usable | `gencode_mane` | MANE v1.5 summary and GTF are acquired and bridge Ensembl, RefSeq, HGNC and NCBI IDs for all joins. GENCODE v50 and the genome FASTA are not fetched. | Fetch GENCODE/FASTA when sequence extraction is needed. |
| D15 Resource catalogues and classifications | low | accessible | — | OpenEBench tool counts are derived (bio.tools API is bot-protected). FAIRsharing records are not harvested (API needs an account). | Optional: user supplies a FAIRsharing token. |
| D16 Background documentation and terminology | low | validated | — | Background pages only; no variant-level evidence. | None. |

## Resources that could not be accessed or reused as-is

| Source | Reason | What is needed |
|---|---|---|
| `biotools` | OpenEBench reachable; bio.tools API challenge page returned HTTP 200 with no data. | bio.tools API is behind a bot challenge; OpenEBench is used instead. |
| `dbnsfp` | Site and licence pages readable; bulk download gated by registration. | User registers with an institutional email at dbnsfp.org and requests download links (about 50 GB); licence CC BY-NC-ND 4.0 forbids redistributing derivatives. |
| `decipher` | Access Denied. | User applies for DECIPHER data access. |
| `depmap` | Automated access met a human-verification page; not bypassed. | User downloads DepMap Public files manually from the DepMap downloads page (portal asks automated clients not to scrape), or use Open Targets target_essentiality. |
| `eqtl_catalogue` | FTP listing verified; API retired. | Use FTP (API retired); metadata tables are hosted in a GitHub repository outside this session's scope. |
| `fairsharing` | Web pages rendered. | User creates a FAIRsharing account/token if programmatic harvesting is wanted. |
| `pharmvar` | Home and download pages readable; terms page not found. | Read the PharmVar terms of use manually, then use its API or downloads. |

## Not yet assessed

Candidates for later evaluation (no verified facts recorded yet):

- HPO / Monarch Initiative / MONDO (phenotype and disease ontologies and annotations)
- OMIM (licensed; restricted agreement)
- gnomAD tandem repeats, LoF curation and local-ancestry datasets (listed on the downloads page, contents not opened)
- Haplotype and LD reference panels (HGDP+1kGP, 1000 Genomes high-coverage genotypes)
- Interaction and network resources (IntAct, STRING, BioGRID interactions)
- Epistasis, digenic and oligogenic disease databases
- Protein stability databases and benchmarks (ThermoMutDB, ProteinGym) and PTM databases
- Standalone predictors outside dbNSFP (CADD, REVEL, SpliceAI as distributed by Illumina, ESM1b/popEVE scores)
- ClinGen Allele Registry and ClinGen Variant Curation Expert Panel classifications beyond the files acquired
- Regulatory variant-effect datasets outside MaveDB (MPRA/STARR-seq repositories)
- Ensembl VEP plugin list and per-plugin licences (documentation page returned HTTP 504)
