The Polygenic Score (PGS) Catalog

The PGS FTP is composed of 2 main directories: "metadata/" and "scores/"

The file "pgs_scores_list.txt" at the root of the PGS FTP contains the list of the Polygenic Score IDs.

#----------#
# metadata #
#----------#

Directory containing all the metadata from the PGS Catalog (the "bulk" metadata).
The metadata of the current release is available in different file formats: Excel (several sheets), CSV and an archived version of all (tar.gz).
Metadata of previous releases are available as archives ("previous_releases/").
For some of the larger studies, specific metadata files are available in the "publications/" subdirectory.

Here is a representation of the structure and files:

metadata/
  |
  |-- pgs_all_metadata.xlsx
  |
  |-- pgs_all_metadata_[sheet_name].csv (7 files)
  |
  |-- pgs_all_metadata.tar.gz (xlsx + csv files)
  |
  |-- publications/
  |   |
  |   |-- PGP000244/
  |   .   |
  |   .   |-- PGP000244_metadata.xlsx
  |   .   |
  |   .   |-- PGP000244_metadata_[sheet_name].csv (7 files)
  |   .   |
  |   .   |-- PGP000244_metadata.tar.gz (xlsx + csv files)
  |   .   |
  |   .   |-- archived_versions/
  |   .
  |   .
  |   |-- PGP00XXXX/
  |       |-- ...
  |       |-- ...
  |
  |-- previous_releases/


#--------#
# scores #
#--------#

This directory lists the individual Polygenic Score entries.
It is indexed by Polygenic Score (PGS) ID to allow programmatic access to the Scoring files and metadata for each PGS, e.g. "PGS000001/, PGS000002/, etc"

The PGS Score metadata files (e.g. in "PGS000001/Metadata/") are in the same format at the ones whitin the "metadata" directory with the exception that they only contain data related to one PGS Score entry instead of all the PGS Scores.
Archived versions of the metadata are available if it has changed since the first release of the PGS Score.

Scoring files have been reformatted from the author/paper-provided tables to have consistent column names and formats and places in the `ScoringFiles` sub-directory.
Files are named according to their PGS Catalog Score ID, e.g. "PGS000001.txt.gz".
Please refer to https://www.pgscatalog.org/downloads/#dl_ftp_scoring for more information on the PGS formatted scoring files.

The variant-level information in PGS is often heterogeneously described and may lack chromosome/position information, contain a mix of positions and/or rsIDs, or be mapped to a genome build different from your sample genotypes. To make PGS easier to apply we have created a set of files that contain harmonized variant information (chromosome name and base pair position) and variant identifiers (updated rsID), in commonly used genome builds (GRCh37/hg19 and GRCh38/hg38) to make variant matching and PGS calculation easier. These files are available under the "Harmonized" sub-directory.
Please refer to https://www.pgscatalog.org/downloads/#dl_ftp_scoring_hm_pos for more information on the PGS harmonized scoring files.


Here is a representation of the structure and files:

scores/
  |
  |-- PGS000001/
  |   |
  |   |-- Metadata/
  |   |   |
  |   |   |-- PGS000001_metadata.xlsx
  |   |   |
  |   |   |-- PGS000001_metadata_[sheet_name].csv (7 files)
  |   |   |
  |   |   |-- PGS000001_metadata.tar.gz (xlsx + csv files)
  |   |   |
  |   |   |-- archived_versions/
  |   |
  |   |-- ScoringFiles/
  |       |
  |       |-- PGS000001.txt.gz
  |       |
  |       |-- Harmonized/
  |       |   |
  |       |   |-- PGS000001_hmPOS_GRCh37.txt.gz
  |       |   |
  |       |   |-- PGS000001_hmPOS_GRCh38.txt.gz
  |       |
  |       |-- archived_versions/
  |
  |-- PGS000002/
  .   |-- ...
  .   |-- ...
  .
  |-- PGS00XXXX/
      |-- ...
      |-- ...

