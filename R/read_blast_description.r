# This function reads in the two tBLASTn output files and combines them as one 
# Input is the downloaded text alignment output and the Description alignment output from BLAST
# Output is one data frame that has the scientific name, genome type information combined with the .txt alignment output
# this function yelds the output that is saved as an EXSCEL workbook of multiple sheets for every query BLAST (ParA/ParC/McdA/FlhG/MinD)

read_blast_description <- function(blast_file, description_file){
  df1 <- read_comparison_from_blast(blast_file)
  df2 <- read.csv(description_file) %>% 
    rename("name2"= "Accession") %>% 
    select(1,2,9)
  
  df2$genome_type <- geneome_type(Description = df2$Description)
  
  combine_df <- merge(df2, df1, by = "name2") %>% 
    filter(., Scientific.Name!="") %>%
    select(-Description) %>% 
    group_by(.,Scientific.Name, genome_type)

  
  return(as.data.frame(combine_df))
  
}
