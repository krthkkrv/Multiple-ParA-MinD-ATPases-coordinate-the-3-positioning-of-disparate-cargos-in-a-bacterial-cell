pipeline for data generated in Consensus_based_tBLASTn

(1) tBLASTn
	-	Using the consensus sequences from Consensus_seq/ tBLASTn search was done against RefSeq Representative genomes (refseq presentative geomes) database, with each of 21 phyulum (proteinbacteria as its different classes), with max target sequences as 5000 and E.value threshold at 0.0001
	-	alignment .txt file and description file .csv downloaded for each search
	-	files for each protein saved under individual folders names *_tBLASTn_hits/

(2) getting Assembly acession number and organism names for all represntative genomes (to get the 12607 number; used in Counts and ratio stats)
	-	following code using on command-line
esearch -query '"bacteria"[filter] AND "latest refseq"[filter] AND ("reference genome"[filter] OR "representative genome"[filter])' -db assembly | esummary | xtract -pattern DocumentSummary -def "NA" -element Organism -element AssemblyAccession 
	-	file saved as ../AssemblyAccession.csv 
  
(2) Running Rscript.r
