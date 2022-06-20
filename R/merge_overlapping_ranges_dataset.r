#  This function is used to check for overlapping hits for individual ParA-like proteins. Used when combining the individual tBLASTn analysis using the 5 consensus sequence to yield one final ParA hits file.
#  Input is two data frame (df1, df2) generated using the read_blast_description (ex. the hits output for two query (ParA and ParC) for one phylum)
#  Uses filter_overlappping_ranges() function
#  Output is one dataframe that has hits unique and common to both querys (ex. ParA and ParC unique hits and hits overlapping both querys)

merge_overlapping_ranges_dataset <- function(data_1, data_2){
  step_1 <- filter_overlapping_ranges(df1 =data_1, df2 =data_2, get_common = T, get_df2_unique = T)
  common_1 <- step_1[[1]]
  common_2 <- step_1[[2]]

  if(is.null(common_1)){
    output <-  rbind(data_1,data_2) %>% 
             arrange(.,name2, start2)
  }else if(dim(common_1)[1]>0){
      common_1$name1 <- paste(data_1$name1[1], data_2$name1[1],sep = "/")
      unique_df1 <- data_1[!data_1$start2 %in% common_1$start2,]
      unique_df2 <- data_2[!data_2$start2 %in% common_2$start2,]
      output <-  rbind(unique_df1,rbind(unique_df2, common_1)) %>% 
             arrange(.,name2, start2)
  }

  return(output)
  
}
