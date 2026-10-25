# File descriptions


### AlphaMissense_hg19.tsv.gz, AlphaMissense_hg38.tsv.gz


Predictions for all possible single nucleotide missense variants (71M) from 19k human protein-coding genes (canonical transcripts) for both hg19 and hg38 coordinates. These files are sorted by genomic coordinates.




### AlphaMissense_gene_hg19.tsv.gz, AlphaMissense_gene_hg38.tsv.gz


Gene-level average predictions, which were computed by taking the mean alphamissense_pathogenicity over all possible missense variants in a transcript (canonical transcript). 




### AlphaMissense_aa_substitutions.tsv.gz


Predictions for all possible single amino acid substitutions within 20k UniProt canonical isoforms (216M protein variants). These are a superset of the amino acid substitutions induced by single nucleotide missense variants. This file uses UniProt accession numbers for proteins and does not have genomic coordinates. 




### AlphaMissense_isoforms_hg38.tsv.gz


Predictions for all possible missense variants for 60k non-canonical transcript isoforms (hg38, GENCODE V32). This file has transcript_id but no UniProt accession numbers. Predictions for non-canonical isoforms were not thoroughly evaluated and should be used with caution. This file is sorted by genomic coordinates.




### AlphaMissense_isoforms_aa_substitutions.tsv.gz


Predictions for all possible single amino acid substitutions for 60k non-canonical transcript isoforms (GENCODE V32). These are a superset of the amino acid substitutions induced by single nucleotide missense variants.This file has transcript_id but no UniProt accession numbers. 


All transcript annotations are based on GENCODE V27 (hg19) or V32 (hg38). Canonical transcripts are defined as described in the [publication](https://www.science.org/doi/10.1126/science.adg7492).


All files are compressed with bgzip.
## Column descriptions 


Note that not all columns are present in every file.




**CHROM**


The chromosome as a string: chr<N>, where N is one of [1-22, X, Y, M].


**POS**


Genome position (1-based).


**REF**


The reference nucleotide (GRCh38.p13 for hg38, GRCh37.p13 for hg19).


**ALT**


The alternative nucleotide.


**genome**


The genome build, hg38 or hg19.


**uniprot_id**


UniProtKB accession number of the protein in which the variant induces a single amino-acid substitution ([UniProt release 2021_02](https://www.uniprot.org/release-notes/2021-04-07-release)).


**transcript_id**


Ensembl transcript ID from GENCODE V27 (hg19) or V32 (hg38).


**protein_variant**


Amino acid change induced by the alternative allele, in the format <Reference amino acid><POS_aa><Alternative amino acid> (e.g. V2L). POS_aa is the 1-based position of the residue within the protein amino acid sequence.


**am_pathogenicity**


Calibrated AlphaMissense pathogenicity scores (ranging between 0 and 1), which can be interpreted as the predicted probability of a variant being clinically pathogenic.


**am_class**


Classification of the protein_variant into one of three discrete categories: 'likely_benign', 'likely_pathogenic', or 'ambiguous'. These are derived using the following thresholds: 'likely_benign' if alphamissense_pathogenicity < 0.34; 'likely_pathogenic' if alphamissense_pathogenicity > 0.564; and 'ambiguous' otherwise.


**mean_am_pathogenicity**


The average alphamissense_pathogenicity of all missense variants per transcript.




## Citation/license and disclaimer


AlphaMissense Database Copyright (2023) DeepMind Technologies Limited. All predictions are provided for non-commercial research use only under [CC BY-NC-SA](https://creativecommons.org/licenses/by-nc-sa/4.0/) license. Researchers interested in predictions not yet provided, and for non-commercial use, can send an expression of interest to alphamissense@google.com.


**Disclaimer**: The AlphaMissense Database and other information provided on this site is for theoretical modelling only, caution should be exercised in use. It is provided “as-is” without any warranty of any kind, whether express or implied. For clarity, no warranty is given that use of the information shall not infringe the rights of any third party. The information provided is not intended to be a substitute for professional medical advice, diagnosis, or treatment, and does not constitute medical or other professional advice.


The predictions in the AlphaMissense Database are predictions only, with varying levels of confidence and should be interpreted carefully. 


If you use this resource for your research please cite the [following publication](https://www.science.org/doi/10.1126/science.adg7492):


“Accurate proteome-wide missense variant effect prediction with AlphaMissense”


Jun Cheng, Guido Novati, Joshua Pan, Clare Bycroft, Akvilė Žemgulytė, Taylor Applebaum, Alexander Pritzel, Lai Hong Wong, Michal Zielinski, Tobias Sargeant, Rosalia G. Schneider, Andrew W. Senior, John Jumper, Demis Hassabis, Pushmeet Kohli, Žiga Avsec 




## Data format samples


`AlphaMissense_hg19.tsv.gz, AlphaMissense_hg38.tsv.gz`
```
#CHROM  POS        REF    ALT       genome        uniprot_id  transcript_id      protein_variant  am_pathogenicity  am_class
chr1    69094      G      T         hg38          Q8NH21      ENST00000335137.4  V2L              0.2937            likely_benign
chr1    69094      G      C         hg38          Q8NH21      ENST00000335137.4  V2L              0.2937            likely_benign
chr1    69094      G      A         hg38          Q8NH21      ENST00000335137.4  V2M              0.3296            likely_benign
chr1    69095      T      C         hg38          Q8NH21      ENST00000335137.4  V2A              0.2609            likely_benign
```


`AlphaMissense_aa_substitutions.tsv.gz`
```
uniprot_id        protein_variant        am_pathogenicity        am_class
A0A024R1R8        M1A        0.4673        ambiguous
A0A024R1R8        M1C        0.3828        ambiguous
```


`AlphaMissense_gene_hg19.tsv.gz, AlphaMissense_gene_hg38.tsv.gz`
```
transcript_id      mean_am_pathogenicity
ENST00000000233.5  0.7422697635438503
ENST00000000412.3  0.37834258163288265
ENST00000001008.4  0.4222901115567318
ENST00000001146.2  0.4666058543393151
```


`AlphaMissense_isoforms_hg38.tsv.gz`
```
#CHROM  POS        REF    ALT       genome        transcript_id      protein_variant  am_pathogenicity  am_class
chr1    65568      A      C         hg38          ENST00000641515.2  K2Q              0.0938            likely_benign
chr1    65568      A      G         hg38          ENST00000641515.2  K2E              0.0766            likely_benign
chr1    65569      A      G         hg38          ENST00000641515.2  K2R              0.0756            likely_benign
chr1    65569      A      T         hg38          ENST00000641515.2  K2M              0.1732            likely_benign
```


`AlphaMissense_isoforms_aa_substitutions.tsv.gz`
```
transcript_id       protein_variant  am_pathogenicity  am_class
ENST00000000442.11  M1A              0.2808            likely_benign
ENST00000000442.11  M1C              0.1724            likely_benign
ENST00000000442.11  M1D              0.7278            likely_pathogenic
ENST00000000442.11  M1E              0.5328            ambiguous
```