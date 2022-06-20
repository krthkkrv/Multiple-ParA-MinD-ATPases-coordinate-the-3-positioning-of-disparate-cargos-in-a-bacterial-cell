Pipeline for data generated in Consensus_based_tBLASTn

(1) tBLASTn
	-	Using the consensus sequences from "Consensus_seq/", tBLASTn search was done against RefSeq Representative genomes (refseq representative genomes) databaset with each of the 21 phyla (proteobacteria as its different classes), with max target sequences as 5000 and E.value threshold at 0.0001
	-	Alignment.txt file and description_file.csv was downloaded for each search
	-	Files for each protein were saved under individual folders names *_tBLASTn_hits/

(2) Getting Assembly accession number and organism names for all representative genomes (to get the 12607 number; used in Counts and ratio stats)
	-	following code using on command-line
	esearch -query '"bacteria"[filter] AND "latest refseq"[filter] AND ("reference genome"[filter] OR "representative genome"[filter])' -db assembly | esummary | xtract -pattern DocumentSummary -def "NA" -element Organism -element AssemblyAccession
	-	File saved as ../AssemblyAccession.csv 

(3) Getting taxonomy information for all representative genomes.
	- 	esearch was used to get taxID for all representative genomes and saved under taxids.txt
	-	Taxonomy information was found using the NCBI Taxonomy Toolkit - TaxonKit on the command line
	taxonkit lineage -t taxids.txt 
	-	File saved as ../Taxonomy_info.xlsx
  
(2) Run ../R/Rscript.r
