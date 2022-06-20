#  This function is used to check for overlapping hits for individual ParA-like proteins. Used when combining the individual tBLASTn analysis using the 5 consensus sequence to yield one final ParA hits file.
#  Input is two data frame (df1, df2) generated using the read_blast_description (i.e. the hits output for two query (ParA and ParC) for one phylum)
#  Output is a list of data frames
  #  data frame 1 is hit ranges from df1 that overlapped with df2 (i.e hits from ParA table that overlapped with ParC)
  #  data frame 1 is hit ranges from df2 that overlapped with df1 (i.e hits from ParC table that overlapped with ParA)
filter_overlapping_ranges <- function(df1 , df2, get_common = F, get_df2_unique = T){
  df1 <- df1%>% 
    arrange(.,name2, start2)
  print(paste("univ df:", dim(df1)))
  df2 <- df2[df2$name2 %in% df1$name2,]%>% 
    arrange(.,name2, start2)
  print(paste("filter df:", dim(df2)))

  acc <- unique(df1$name2)
  filtered_list <- NULL
  filtered_list_2 <- NULL

  for(i in 1:length(acc)){
  subset_1 <- df1[df1$name2 == acc[i],]
  subset_2 <- df2[df2$name2 == acc[i],]
  if(dim(subset_1)[1]==0 && dim(subset_2)[1]>0){
    df2_unique <- rbind(df2_unique, subset_2)
  }
  if(dim(subset_1)[1] > 0 && dim(subset_2)[1] > 0){
    pos_1 <- subset_1 %>% 
      select(.,start2, end2) %>% 
      as.matrix()
    pos_2 <- subset_2 %>% 
      select(.,start2, end2) %>% 
      as.matrix()
    for(i in 1:dim(pos_2)[1]){
      tmp <- pos_2[i,]
      for (i in 1:dim(pos_1)[1]){
        if ( tmp[1] %gl% c(pos_1[i,1] , pos_1[i,2]) || tmp[1] %gl% c(pos_1[i,2], pos_1[i,1])){
          filtered_list <- rbind(filtered_list, subset_1[i,])
          filtered_list_2 <- rbind(filtered_list_2, subset_2[subset_2$start2==tmp[1],])
        }else if(tmp[2] %gl% c(pos_1[i,1] , pos_1[i,2]) || tmp[2] %gl% c(pos_1[i,2], pos_1[i,1])){
          filtered_list <- rbind(filtered_list, subset_1[i,])
          filtered_list_2 <- rbind(filtered_list_2, subset_2[subset_2$start2==tmp[1],])
        }else if(tmp[1] <= pos_1[i,1] && tmp[2] >= pos_1[i,2] || tmp[1] >= pos_1[i,1] && tmp[2] <= pos_1[i,2]){
          filtered_list <- rbind(filtered_list, subset_1[i,])
          filtered_list_2 <- rbind(filtered_list_2, subset_2[subset_2$start2==tmp[1],])
        }else if(sum(tmp == pos_1[i,]) >0){
          filtered_list <- rbind(filtered_list, subset_1[i,])
          filtered_list_2 <- rbind(filtered_list_2, subset_2[subset_2$start2==tmp[1],])
        }else{
          next
             }
        }
      }
    }
  }

    my_list <- list("common_1" = filtered_list, "common_2" = filtered_list_2)

  
  if(get_df2_unique == T && get_common ==T){
    return(my_list)
  }else if(get_common ==T){
    return( filtered_list)
  }else if(get_df2_unique ==T){
      return(df2_unique)
  }else{
      print("OUTPUT not selected")
    }

}
