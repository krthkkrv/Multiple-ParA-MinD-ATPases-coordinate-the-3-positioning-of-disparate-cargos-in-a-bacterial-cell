Pipline for aquiring cargo specific consensus sequences

(1) BLASTp
	-	Using the ParA-like proteins from H. neapolitanus as the query, a protein BLAST search was done against the RefSeq Selectproteins (Refseq select)database.
	-	Alignmnet .txt file and Description file .csv is saved in alignment_hits folder under each protein

(2) Gene neighborhood analysis
	-	The 100 proteins from BLASTp (saved in query.txt) were then fed into FlaGs analysis (http://130.239.193.227/html/webFlaGs.html)
	-	From these, proteins with similar flanking genes were selected
	-	Analysis results saved under GNA folder

(3) MSA
	-	The filtered proteins were then aligned using NCBI's MSA tool (https://www.ncbi.nlm.nih.gov/tools/cobalt/re_cobalt.cgi)
	-	The MSA .fasta file was then fed into EMBOSS to get a consensus sequence (https://www.bioinformatics.nl/cgi-bin/emboss/cons)
	-	Consensus sequence from EMBOSS saved in *_consenss.txt file
	-	Files saved under MSA folder

(4) Conensus_MSA
	-	All consensus sequence for te 6 ParA-like proteins saved in Consensus.txt file 
	
