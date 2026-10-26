# Release 1.0 variant callset
Prepared by Peter Audano and Carolyn Paisie

## Description
Variants from 65 phased genomes were called by PAV 2.4.0.1 (development version
before the projected 3.0 release). Merging was performed with SV-Pop. Support
intersects against a variety of other callers was done and used for QC.
Unreliable variant regions including centromeres and telomeres were discarded.

## References
Variant calls were generated on two references using identical data, tools,
and versions:

* GRCh38: GRCh38-NoALT. Primary assembly only (chromosome scaffolds, unplaced contigs, and unlocalized contigs)
* T2T-CHM13: The T2T CHM13 reference version 2.0 (T2T-CHM13v2.0, or hs1 on UCSC).

## Data
There is a subdirectory for each reference.

The VCF files in the reference directory (i.e. not in a subdirectory) contain
2 phased genotypes per sample for each variant call where the first allele is
"haplotype 1" (h1) and the second is "haplotype 2" (h2) in the phased assemblies
produced for this project. For example, genotype "0|1" means the variant was
found in the h2 assembly, but not the h1 assembly. A "." in a genotype allele
means the locus was not callable because assemblies were not aligned there. For
example, genotype ".|1" means the variant was not callable in h1.

The variant annotation table Excel file contains descriptions for the fields
in tables "annotation_table" subdirectory for each reference.

## Subdirectories

* annotation_table: The same variants in the VCFs in tabular format.
* vcf_with_unphased: Same VCF files as the reference directory with three
  haplotypes derived from the phased assemblies: h1, h2, and unphased. The
  genotype order is the same with a third genotype for unphased.
* complex_sv: Merged Complex SVs called with a development version of PAV.
