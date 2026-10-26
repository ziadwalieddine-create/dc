"""High-Value Human Variant Reference Index (vindex).

Source records -> normalised entities -> evidence rows (parallel axes, no composite score) -> matcher for called variants.
Modules:
  refgenome   reference FASTA access, contig naming, assembly validation
  normalise   allele normalisation (left-aligned, parsimonious), atomisation, identity keys
  vcfio       sample VCF/gVCF reading: multiallelic splitting that keeps genotype, phase and quality
  callability reference-confidence providers (gVCF blocks, callable-region BED, read pileup)
  schema      table schemas and controlled vocabulary of the index
"""
SCHEMA_VERSION = '0.1.0'
