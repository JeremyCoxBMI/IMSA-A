# IMSA-A
IMSA+A is an extension of Arron-IMSA program for metagenomic taxonomic classificiation using RNAseq reads. 
# Installation
IMSA+A should be installed into it's own directory, independent of IMSA.
See Detailed Directions.
## July 2017 Update
* Expanded tips and tricks in directions
* Include example results to validate pipeline
## December 2016 update
* Starting with DEC 2016 update, IMSA code is distributed with IMSA+A.
* Very minor updates for compatability are included; said changes are documented with comments.
* A new script "postprocesscount4acc.py" allows use of accession number data.
New options in systemSettings.py, start with "ACC_" to point to a new accession version to taxon database.
# GI number phase out; run IMSA+A using either
https://www.ncbi.nlm.nih.gov/news/03-02-2016-phase-out-of-GI-numbers/
GI numbers are being phased out (09/2106).  IMSA and IMSA+A counting scripts count using GI numbers.
Any files downloaded here together will function, but updated NCBI databases after 09/2016 will not.
The hand off to accession.version should be relatively painless and require very minor code changes.  However, IMSA+A will need NCBI to release an accession.version database dump in the same way that they do a GI dump.  Until then, IMSA+A is not future compatible.
## 12/2016 update
## 02/2017 update
We found a bug in IMSA paired end processing.  We are working to correct and will release soon.
