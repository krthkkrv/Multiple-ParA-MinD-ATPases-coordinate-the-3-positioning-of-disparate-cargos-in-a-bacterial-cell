## Running the blast inputs
#  list of phylum names
phylum_list <- c( "Acidobacteria","Actinobacteria","Alphaproteobacteria", "Aquificae", "Bacteroidetes", "Betaproteobacteria", "Chlamydiae", "Chlorobi",'Chloroflexi','Cyanobacteria','Deinococcales', 'Deltaproteobacteria', 'Epsilonproteobacteria', 'Fibrobacteres','Firmicutes', 'Fusobacteria', 'Gammaproteobacteria', 'Planctomycetes', "Spirochaetes", 'Thermotogae', 'Zetaproteobacteria')

#  list of query proteins
protein_list <- c("ParA", "ParC", "McdA", "FlhG", "MinD")

#  the main path to the subfolders holding the blast output downloads
path <- "Consensus_based_tBLASTn/"

#  Creating a empty list
Consensus_tBLASTn <- NULL

#  combining the blast output in a phylum specific manner and adding to the empty list
for(i in 13:length(phylum_list)){
  phylum <- phylum_list[i]
  print(phylum)
  names(Consensus_tBLASTn)[[i]] <- phulym
  Consensus_tBLASTn[[i]] <- as.data.frame(condensed_phylum(phylum = phylum, path = "Consensus_based_tBLASTn/", protien = protein_list))
}



## Generating the hits_ratio_table
#  creating an empty list
binary_table <- NULL
#  looping thorough the listss of data frame to generate one final data frame with columns specifying the number of parA hits per species, the hits in heat genetic element, the ration of hits/genetic_element, as well as a count of different ParA-like protein hits found in a species.

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
    
    #cheomosome.accession
    g.acc <- append(g.acc, ifelse(length(unique(sub[sub$genome_type=="chromosome", 4]))>0,unique(sub[sub$genome_type=="chromosome", 4]),"NA"))
    
    # total.ParA.hits
    sp.hit.count <- append(sp.hit.count, unique(sub$sp_count))
    
    # ParAs.in.plasmid/chromosome
    genome_type <- sub$genome_type
    plas.para <- append(plas.para, sum(genome_type=="plasmid"))
    chr.para <- append(chr.para, sum(genome_type == "chromosome"))
    
    # nmbr.genetic.elements
    nmbr_genetic_elements <- append(nmbr_genetic_elements, length(unique(sub$name2)))

    #plasmid/chromosome
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
    
    #unknown.cargo.ParAs
    n=0
    for(i in 1:length(name1)){
      if(nchar(name1[i])>4){
        n=n+1}
    }
    other_parAs <- append(other_parAs, n)
  }
  #knitting all the lists
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
binary_table$ratio <- binary_table$total.ParA.hits/binary_table$nmbr.genetic.elements


# cleaning the scientifc.Name column so that they match the naming format in the AssemeblyAccession.csv file
tmp <- binary_table %>% 
    separate(Scientific.Name, c("name1", "name2"), extra = "drop", fill = "right", remove = T)
  
tmp1 <- na.omit(tmp)
tmp1$Scientific.Name <- paste(tmp1$name1, tmp1$name2, sep = " ")
binary_table <- tmp1[,c(17,3,4,5,6,7,9,8,16,10:15)]

# calculating the single.cargo.hits
binary_table$single.cargo.hits <- binary_table$ParA+binary_table$ParC+binary_table$McdA+binary_table$MinD+binary_table$FlhG

# reading in Assembly accession file 
AssemblyAccession <- read.csv("../Input_files/AssemblyAccession.csv")

# merging in the Assembly accession numbers
binary_table_meta <- merge(AssemblyAccession, binary_table, by = "Scientific.Name", all.y = T)

# getting taxid
query <- sort(binary_table$chromosome.accession)[9:11905]
quer_op <- esummary_fetch_taxid(query)
colnames(quer_op)[1] <- "chromosome.accession"

# adding the taxID information to the data
binary_table_meta <- merge(quer_op, binary_table, by = "chromosome.accession", all.y = T)

# reading in the taxonomy file (has info for phylum,class....)
taxonomy <- read_xlsx("Taxonomy_info.xlsx") %>% 
  as.data.frame()

# merging the taxonomy info 
binary_table_meta <- merge(taxonomy[,c(1,3,4)], binary_table_meta, by = "taxId", all.y = T)

#  Wrting out the final output
write.csv(binary_table_meta, "Consensus_based_tBLASTn/ratio_and_binary_info_for_parA_hits_4.csv", row.names = F)

--complete--

## Statistics
#  number/percent of species having atleast 1 ParA hits
hits_ratio_file <- read.xlsx("Consensus_based_tBLASTn/ratio_and_binary_info_for_parA_hits.xlsx", sheetIndex = 1)
total_sum <- dim(hits_ratio_file[hits_ratio_file$total.ParA.hits>=1,])[1] # change the >=1 to any number to filter out the data and count; for example to get number of species with atleast 2 ParA hits, change the number to >=2

#  the total number  of reference genome as on feb 4 2021 is 12607
(total_sum/12607)*100
total_sum

#  number/percent of species having the hits/genomic_element ratio of atleast 1 
hits_ratio_file <- read.xlsx("Consensus_based_tBLASTn/ratio_and_binary_info_for_parA_hits.xlsx", sheetIndex = 1)
total_sum <- dim(hits_ratio_file[hits_ratio_file$ratio>=1,])[1] # change the >=1 to any number to filter out the data and count; for example to get number of species with a hits/genomic_element ratio of atleast 2, change the number to >=2

#  the total number  of reference genome as on feb 4 2021 is 12607
(total_sum/12607)*100
total_sum


