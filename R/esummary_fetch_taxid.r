# Used to fetch the taxonomy ID for a given list of accession number 
# Input is a list of accession numbers
# Output is  data frame of two collumns; Accession number and taxonomy ID
esummary_fetch_taxid <- function(access.nb){

    N <- length(access.nb)
    a <- 1
    b <- if(N<100){
      N
    }else{
      100
    }
    taxid <- NULL
    Access.nb <- NULL
    repeat{ 
      # get url
      URL <- paste("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=nucleotide&id=", 
                   paste(access.nb[a:b], collapse = ","), sep = "")
      # downlaod the file
      xml2::download_xml(url = URL, file = "try.xml")
      # read in the xmlfile fro parsing
      data <- xmlInternalTreeParse("try.xml")
      top <- xmlRoot(data)
      # get the description nodes
      taxid <- append(taxid, unlist(xpathApply(top,  "//DocSum/Item[@Name = 'TaxId']", xmlValue)))
      # get the accession nodes
      Access.nb <- append(Access.nb, unlist(xpathApply(top, "//DocSum/Item[@Name='Caption']", xmlValue)))
      if(b==N) break
      a <- b+1
      b <- b+100
      if (b>N) {b <- N}
    }
    
    df <- data.frame("name2" = Access.nb,
                     "taxId" = taxid)
    return(df)
}
