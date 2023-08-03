# Load following packages
library("tidyverse")
library("genoPlotR")
library("reshape2") 
library("readxl")
library("xml2")
library("xlsx")
library("extraoperators")

## Requires following functions
# hit_count()
# esummary_fetch_taxid()
# read_blast_description()
# filter_overlapping_ranges()
# merge_overlapping_ranges_dataset()
# condensed_phylum()


## Running the tblastn inputs
#  List of phylum names
phylum_list <- c( "Acidobacteria","Actinobacteria","Alphaproteobacteria", "Aquificae", "Bacteroidetes", "Betaproteobacteria", "Chlamydiae", "Chlorobi",'Chloroflexi','Cyanobacteria','Deinococcales', 'Deltaproteobacteria', 'Epsilonproteobacteria', 'Fibrobacteres','Firmicutes', 'Fusobacteria', 'Gammaproteobacteria', 'Planctomycetes', "Spirochaetes", 'Thermotogae', 'Zetaproteobacteria')

#  List of query proteins
protein_list <- c("ParA", "ParC", "McdA", "FlhG", "MinD")

#  The main path to the subfolders holding the tblastn output downloads
path_to <- "../Consensus_based_tBLASTn/"

#  Creating a empty list
Consensus_tBLASTn <- NULL

#  Combining the blast output in a phylum specific manner and adding to the empty list
for(i in 13:length(phylum_list)){
  phylum <- phylum_list[i]
  print(phylum)
  names(Consensus_tBLASTn)[[i]] <- phulym
  Consensus_tBLASTn[[i]] <- as.data.frame(condensed_phylum(phylum = phylum, path = path_to, protein = protein_list))
}

## 795 bacterial species did not have any hits. To confirm this, a local database of the reference genome assemblies were created and the consensus sequence for the 5 ParA/MinD family proteins was done using BLAST+ tblastn
## In the terminal

#cat assembly_query_list.txt | while read -r acc ; do
#    esearch -db assembly -query $acc </dev/null \
#        | esummary \
#        | xtract -pattern DocumentSummary -element FtpPath_GenBank \
#        | while read -r url ; do
#            fname=$(echo $url | grep -o 'GCA_.*' | sed 's/$/_genomic.fna.gz/') ;
#            wget "$url/$fname" ;
#       done ;
#   done

#gzip -d GCA_*

#cat GCA_*  > ref_seq_genome_subset.fasta
#makeblastdb -in ref_seq_genome_subset.fasta -out my_reference -parse_seqids -dbtype nucl


# tblastn -db my_reference -evalue 0.0001 -query query.fasta -out no_hits_tblastn_results.txt -outfmt "6 qseqid aasinames ataxids sacc pident length evalue slen sstart send qstart qend" ->
# useing esearch to also get scientific name and genome type information for these hits ->
####
#Back to R

# reading the two files in
no_hits <- read.csv("../Consensus_based_tBLASTn/no_hits_tblastn_results.txt", sep = "\t", header = F)
colnames(no_hits1400) <- c( "start1","end1", "start2","end2", "name1" , "name2", "per_id","aln_len", "mism" ,"gaps", "e_value", "bit_score" )

accession_meta_2 <- read.csv("../Consensus_based_tBLASTn/nucleotide_accession_output_2.txt", sep = "\t", header = F)
colnames(accession_meta_2) <- c("Accession", "name2", "Scientific.Name", "genome_type")

no_hits <- merge(accession_meta_2, no_hits, by = "name2", all.y = T)

# adding the hits from the second round of blast analysis to the Consensus_tBLASTn list.

tmp <- no_hits[no_hits$genome_type=="",]
tmp$genome_type <- c("genomic")

no_hits <- rbind(no_hits[!no_hits$genome_type=="",], tmp) %>% 
  select(-Accession)

ParA_df <- no_hits[no_hits$name1 == "ParA",]
ParC_df <- no_hits[no_hits$name1 == "ParC",]
McdA_df <- no_hits[no_hits$name1 == "McdA",]
FlhG_df <- no_hits[no_hits$name1 == "FLhG",]
MinD_df <- no_hits[no_hits$name1 == "MinD",]

if(!is.null(ParA_df) && !is.null(ParC_df)){
  file_1 <-  merge_overlapping_ranges_dataset(data_1 = ParA_df, data_2 = ParC_df)
}else{
  file_1 <- "ParA"
}

if(!is.null(McdA_df)){
  file_2 <- merge_overlapping_ranges_dataset(data_1 = file_1, data_2 = McdA_df)
}else{
  file_2  <- file_1
}

if(is.null(FlhG_df)){
  file_3  <- file_2
}else{
  file_3 <- merge_overlapping_ranges_dataset(data_1 = file_2, data_2 = FlhG_df)
}

if(!is.null(MinD_df)){
  file_4 <- merge_overlapping_ranges_dataset(data_1 = file_3, data_2 = MinD_df)
}else{
  file_4  <- file_3
}

sp_count <- hit_count(df = file_4, count_clmn = file_4$Scientific.Name) %>% 
  rename("sp_count" = "nmbr_of_hits", "Scientific.Name" = "object")
acc_count <- hit_count(df = file_4, count_clmn = file_4$name2) %>% 
  rename( "acc_count" = "nmbr_of_hits", "name2" = "object")

file_4_raw <- merge(merge(file_4, sp_count, by = "Scientific.Name"),acc_count, by = "name2")

file_4_raw_chr <- file_4_raw[file_4_raw$genome_type=="genomic",]

sp_list <- unique(file_4_raw_chr$Scientific.Name)
to_remove <- NULL
for(i in 1:length(sp_list)){
  sp <- sp_list[i]
  subset_1 <- file_4_raw_chr[file_4_raw_chr$Scientific.Name == sp,]
  sp_c <- unique(subset_1$sp_count)
  acc_c <- unique(subset_1$acc_count)
  if(length(acc_c)==1){
    if(acc_c != sp_c){
      subset_1 <- subset_1[order(sort(subset_1$acc_count), decreasing = T),]
      to_remove <- rbind(to_remove, subset_1[subset_1$name2!=subset_1$name2[1],])
    }
  }else{
    subset_1 <- subset_1[order(sort(subset_1$acc_count), decreasing = T),]
    to_remove <- rbind(to_remove, subset_1[subset_1$name2!=subset_1$name2[1],])
  }
  
}
file_4_filter_chr <- file_4_raw_chr[!file_4_raw_chr$name2 %in% to_remove$name2,]

final <- file_4_raw %>% 
  filter(.,genome_type != "genomic") %>% 
  rbind(.,file_4_filter_chr) %>% 
  select(-sp_count, -acc_count) %>% 
  as.data.frame(.)

final <- distinct(final)

sp_count <- hit_count(df = final, count_clmn = final$Scientific.Name) %>% 
  rename("sp_count" = "nmbr_of_hits", "Scientific.Name" = "object")

final <- merge(final, sp_count, by = "Scientific.Name") %>% 
  select(Scientific.Name, genome_type, name1, name2, everything()) %>% 
  group_by(Scientific.Name, genome_type)


SN_list <- NULL
for (i in 1:21){
  SN_list <- append(SN_list, Consensus_tBLASTn[[i]]$Scientific.Name)
}
SN_list <- unique(SN_list)

Consensus_tBLASTn$no_hits <- final[!final$Scientific.Name%in%SN_list,]


Consensus_tBLASTn$no_hits$direction <- ""


## Generating the hits_ratio_table
#  Creating an empty list
binary_table <- NULL
#  Looping thorough the lists of data frame to generate one final data frame with columns specifying the number of ParA/MinD hits per species, the hits in each genetic element, the number of unique ParA/MinD hits in each species, the ratio of hits/genetic_element, as well as a count of different ParA/MinD-like protein hits found in a species.

for(i in 1:length(Consensus_tBLASTn)){
  df <- as.data.frame(Consensus_tBLASTn[[i]])
  name2_list <- unique(df$Scientific.Name)
  sp.Name <- NULL
  parA <- NULL
  parC <- NULL
  mcdA <- NULL
  minD <-  NULL
  flhG <- NULL
  other_parAs <- NULL
  nmbr_genetic_elements <- NULL
  sp.hit.count <- NULL
  g.acc <- NULL
  plasmid <- NULL
  chromosome <- NULL
  plas.para <- NULL
  chr.para <- NULL
  
  for (i in 1:length(name2_list)){
    sp.Name <- append(sp.Name, name2_list[i])
    sub <- subset(df, df[,"Scientific.Name"] == name2_list[i])
    
    #chromosome.accession
    g.acc <- append(g.acc, ifelse(length(unique(sub[sub$genome_type=="chromosome", 4]))>0,unique(sub[sub$genome_type=="chromosome", 4]),"NA"))
    
    # total.ParA.hits
    sp.hit.count <- append(sp.hit.count, unique(sub$sp_count))
    
    # ParAs.in.plasmid/chromosome
    genome_type <- sub$genome_type
    plas.para <- append(plas.para, sum(genome_type=="plasmid"))
    chr.para <- append(chr.para, sum(genome_type == "chromosome"))
    
    # nmbr.genetic.elements
    nmbr_genetic_elements <- append(nmbr_genetic_elements, length(unique(sub$name2)))
    
    # plasmid/chromosome
    sub_1 <- subset(sub, sub[,"genome_type"] == "plasmid")
    plasmid <-  append(plasmid,length(unique(sub_1$name2)))
    sub_1 <- subset(sub, sub[,"genome_type"] == "chromosome")
    chromosome <- append(chromosome, length(unique(sub_1$name2)))
    
    # ParA/ParC/McdA/MinD/flhG
    name1 <- sub$name1
    parA <- append(parA, sum(name1 == "ParA"))
    parC <- append(parC, sum(name1 == "ParC"))
    mcdA <- append(mcdA, sum(name1 == "McdA"))
    minD <- append(minD, sum(name1 == "MinD"))
    flhG <- append(flhG, sum(name1 == "FlhG"))
    
    # unknown.cargo.ParAs
    n=0
    for(i in 1:length(name1)){
      if(nchar(name1[i])>4){
        n=n+1}
    }
    other_parAs <- append(other_parAs, n)
  }
  # knitting all the lists
  output <- data.frame("Scientific.Name" = sp.Name,
                       "chromosome.accession" = g.acc,
                       "total.ParA.hits" = as.numeric(sp.hit.count),
                       "ParAs.in.plasmid" = plas.para,
                       "ParAs.in.chromosome"= chr.para,
                       "nmbr.genetic.elements" = as.numeric(nmbr_genetic_elements),
                       "plasmid.count" = as.numeric(plasmid),
                       "chromosome.count" = as.numeric(chromosome),
                       "ParA" = parA,
                       "ParC" = parC,
                       "McdA" = mcdA,
                       "MinD" = minD,
                       "FlhG" = flhG,
                       "unknown.cargo.ParAs" = other_parAs
  )
  output <- output %>% 
    as.data.frame()
  
  binary_table <- rbind(binary_table, output)
}

#  calculating the hits/genetic_element ratio
binary_table$ratio <- binary_table$total.ParA.hits-binary_table$nmbr.genetic.elements+1

# calculating the single.cargo.hits
binary_table$single.cargo.hits <- binary_table$ParA+binary_table$ParC+binary_table$McdA+binary_table$MinD+binary_table$FlhG

# Aquiring information on hits, based on assmebly accession numbers in linux
# for g in accessionQuer.txt
# do
# esearch -db assembly -query GCF_005121165.2|efetch -format docsum | xtract -pattern DocumentSummary -element AssemblyAccession Taxid >> all_accession.txt
# done

# converitng file to .csv and importing it
global <- read.csv("../AssemblyAccession.csv", header = F, sep = '\t')
taxonomy <- read_xlsx("../Taxonomy_info.xlsx") %>% 
  as.data.frame()
global <- merge(global, taxonomy[,c(1,3)], all.x= T, by = "taxID")

global <- global[!duplicated(global$AssemblyAccession),]

#  merging in the Assembly accession numbers
binary_table_meta <- merge(global, binary_table, by = "Scientific.Name")

# some scientific name entires had two assmebly accessions. these were manually scrrened thorugh and then added to the final data (supplemtal data in the paper).

--complete--