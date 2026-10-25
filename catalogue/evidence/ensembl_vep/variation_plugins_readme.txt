Date: 10/2025

About
-----

SpliceAI is a deep learning model developed by Illumina that predicts the impact of genetic variants on splicing (https://github.com/Illumina/SpliceAI).
This README explains how Ensembl generates pre-computed SpliceAI scores for use with the VEP SpliceAI plugin (https://github.com/Ensembl/VEP_plugins/blob/main/SpliceAI.pm)

Ensembl provides SpliceAI scores for all possible SNVs (single nucleotide variants) in GRCh38 for all MANE Select transcripts.
The precomputed Ensembl SNV scores allow for rapid annotation of variants without requiring real-time computation.

NOTE: Ensembl does NOT generate SpliceAI scores for insertions and deletions (indels). These are only available from Illumina, see https://github.com/Ensembl/VEP_plugins/blob/main/SpliceAI.pm for Illumina data download instructions.


How did we generate the scores?
------------------------------
- Assembly: GRCh38
- Gene annotation file: includes MANE select transcript for each protein coding gene
- Maximum distance between the variant and gained/lost splice site: 50 (default)

We generated the scores for all possible SNVs overlapping MANE Select transcripts using the SpliceAI tool.
We run the tool to calculate both raw and masked scores. The type of scores are indicated in the file name, example: spliceai_scores.raw.snv.ensembl_mane_v1.4.grch38.vcf.gz includes only raw scores
In the file name is also indicated the MANE version for which the scores where calculated (mane_v1.4).
