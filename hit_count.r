# To count through a column in a data frame and generate a counts list. (for example coutning the number of ParA hits per species)
# Requires a data frame, and mention the column you want to count
# output is a data frame with two columns; the name of the element counted and the count
hit_count <- function(df, count_clmn){
  name2_list <- unique(df[,count_clmn])
  final_output <- data.frame()
  nmbr_of_hits <- NULL
  object_name <- NULL
    for (i in seq_along(name2_list)) {
      # Sub setting based on the unique list
      sub <- subset(df, df[,count_clmn] == name2_list[i])
      object_name <- append(object_name, name2_list[i])
      # Using dim() function to get the number of repeats - rows
      nmbr_of_hits <-append(nmbr_of_hits, dim(sub)[1])
    }
  final_output <- data.frame("object" = object_name,
                             "nmbr_of_hits" = nmbr_of_hits)
  final_output$nmbr_of_hits <- as.numeric(final_output$nmbr_of_hits)
  final_output$object <- as.character(final_output$object)
  return(final_output)
}
