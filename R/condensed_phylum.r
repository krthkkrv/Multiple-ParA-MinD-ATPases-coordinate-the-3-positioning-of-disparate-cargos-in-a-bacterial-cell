# For a given phulym this function merges the hits result form all the 5 ParA-like proteins tBLAStn and filters out the redundant hits.
# Input - phylum name; main path that leads to the folder containing the subfolders (i.e. the query proteins) holding the downloaded tBlastn output files; finally a list of names of the query proteins(in this case it would be a c("ParA", "ParC", "McdA", "FlhG", "MinD"))
# Uses following functions
    # read_blast_description()
    # merge_overlapping_ranges_dataset()
    # filter_overlapping_ranges()
    # hit_count()
# Output is one dataframe showing the ocurrance of hits of ParA-like proteins in the given phylum
condensed_phylum <- function(phylum, path, protein){
    ParA_df <- NULL
    ParC_df <- NULL
    McdA_df <- NULL
    MinD_df <- NULL
    FlhG_df <- NULL
  for(i in 1:length(protein)){
    tryCatch({
    txt_file <- paste(path, protein[i],"/", phylum,"_",protein[i], "-Alignment.txt", sep = "")
    csv_file <- paste(path, protein[i],"/", phylum,"_",protein[i], "-Alignment-Descriptions.csv", sep = "")
    assign(paste(protein[i], "_df", sep = ""), value = read_blast_description(blast_file = txt_file, description_file = csv_file))
    }, error = function(e){cat(protein[i], "data not available", "\n")})
  }
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
  
  sp_count <- hit_count(df = file_4, count_clmn = "Scientific.Name") %>% 
    rename("sp_count" = "nmbr_of_hits", "Scientific.Name" = "object")
  acc_count <- hit_count(df = file_4, count_clmn = "name2") %>% 
    rename( "acc_count" = "nmbr_of_hits", "name2" = "object")
  
  file_4_raw <- merge(acc_count, merge(file_4, sp_count, by = "Scientific.Name"), by = "name2")
  file_4_raw_chr <- file_4_raw[file_4_raw$genome_type=="chromosome",]
  
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
    filter(.,genome_type == "plasmid") %>% 
    rbind(.,file_4_filter_chr) %>% 
    select(-sp_count, -acc_count) %>% 
    as.data.frame(.)
    
  
  sp_count <- hit_count(df = final, count_clmn = "Scientific.Name") %>% 
    rename("sp_count" = "nmbr_of_hits", "Scientific.Name" = "object")
  
  final <- merge(final, sp_count, by = "Scientific.Name") %>% 
    select(Scientific.Name, genome_type, name1, name2, everything()) %>% 
    group_by(Scientific.Name, genome_type)%>% 
    arrange(desc(sp_count, name2, .by_group = T)) 

  return(final)
  
}
