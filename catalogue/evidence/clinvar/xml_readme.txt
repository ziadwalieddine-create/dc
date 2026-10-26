README file for the path ftp.ncbi.nlm.nih.gov/pub/clinvar/xml/
Last updated September 2, 2025


This directory contains XML for ClinVar accessions. 

As of August 7, 2025, ClinVar provides two types of XML files, one for VCV accessions and one for RCV accessions.
	
Each file corresponds to a monthly comprehensive data release of either VCV or RCV accessions.

VCV XML 
- ClinVarVCVRelease is provided in this directory (ftp.ncbi.nlm.nih.gov/pub/clinvar/xml)

RCV XML 
- ClinVarRCVRelease is provided in the RCV_release subdirectory (ftp.ncbi.nlm.nih.gov/pub/clinvar/xml/RCV_release/)

Each file name is formatted like
	ClinVarVCVRelease_YYYY-MM.xml.gz for VCV XML
	ClinVarRCVRelease_YYYY-MM.xml.gz for RCV XML   
where YYYY and MM are the year and month and day the file was released.

Each XML file is updated the first Thurday of the month. A symbolic link  is established to the most recent file:
	ClinVarVCVRelease_00-latest.xml.gz for VCV XML
        ClinVarRCVRelease_00-latest.xml.gz for RCV XML
	
The xsd's for validating the XML is provided in these paths:
	ftp.ncbi.nlm.nih.gov/pub/clinvar/xsd_public/ for VCV XML
	ftp.ncbi.nlm.nih.gov/pub/clinvar/xsd_public/RCV/ for RCV XML


------------------------------------------------------------
Subdirectories
------------------------------------------------------------
sample_xml
----------------------------------------
The sample_xml subdirectory contains the following examples.

An example of a single ClinVarResult-Set, from the VCV XML and as retrieved by efetch, based on VCV000091629.19:
ftp.ncbi.nlm.nih.gov/pub/clinvar/xml/sample_xml/VCV_XML_VCV000091629.xml

An example of a docsum, as retrieved by esummary, for VCV000091629.19:
ftp.ncbi.nlm.nih.gov/pub/clinvar/xml/sample_xml/document_summary_VCV000091629.xml


----------------------------------------
archive
----------------------------------------
This subdirectory contains VCV XML (current format only) from the monthly releases of previous years.

RCV XML from previous years (current format only) is in ftp.ncbi.nlm.nih.gov/pub/clinvar/xml/RCV_release/archive/

Each archive directory has year-specific subdirectories.


----------------------------------------
weekly_release
----------------------------------------
The weekly release of ClinVar data as VCV XML is generated each week as the web site is updated. All the content of the weekly_release 
subdirectory is removed when the monthly release XML is generated at the first of the month.
	
The files for weekly releases are named similarly to the monthly releases except that the day is included in the date, e.g.

	ClinVarVCVRelease_YYYY-MMDD.xml.gz

where YYYY, MM, and DD are the year, month and day the file was created.

Weekly releases for the RCV XML is available in	ftp.ncbi.nlm.nih.gov/pub/clinvar/xml/RCV_release/weekly_release/


----------------------------------------
VCV_xml_old_format/
----------------------------------------

The old format for VCV XML, ClinVarVariationRelease, is in the VCV_xml_old_format subdirectory (ftp.ncbi.nlm.nih.gov/pub/clinvar/xml/VCV_xml_old_format/). 

- This format is no longer supported as of August 7, 2025.
- This format had a single element for classification, rather than separate elements for germline and somatic classifications.
- File names are formatted like ClinVarVariationRelease_YYYY-MM.xml.gz 
- The symbolic link ClinVarVariationRelease_00-latest.xml.gz links to the last file for the old VCV XML
- The XSD for validation is in ftp.ncbi.nlm.nih.gov/pub/clinvar/xsd_public/VCV_xsd_old_format/ 

----------------------------------------
VCV_xml_old_format/beta
----------------------------------------
Versions of the VCV XML that were generated as part of a beta release. As of April 4, 2019, this directory is no longer updated.

----------------------------------------
RCV_xml_old_format/
----------------------------------------

The old format for RCV XML, ClinVarFullRelease, is provided in the RCV_xml_old)format subdirectory (ftp.ncbi.nlm.nih.gov/pub/clinvar/xml/RCV_xml_old_format/).

- This format is no longer supported as of August 7, 2025.
- This format had a single element for classification, rather than separate elements for germline and somatic classifications.
- File names are formatted like ClinVarFullRelease_YYYY-MM.xml.gz for RCV XML
- The symbolic link ClinVarFullRelease_00-latest.xml.gz links to the last file for the old RCV XML
- The XSD for validation is in ftp.ncbi.nlm.nih.gov/pub/clinvar/xsd_public/RCV_xsd_old_format/ for old RCV XML

Please direct any questions or comments to clinvar@ncbi.nlm.nih.gov
