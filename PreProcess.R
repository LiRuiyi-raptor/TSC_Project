# Created by Ruiyi Li on April 2, 2020
# This function return the data after filtering genes
#______________________
#Input:
# data: a numberic matrix with the rows are genes and colums are cells
# th : The parameter to filter genens. If the gene is not expressed more than th% cells, it will be removed.

#Output
# data: data after filtering uninformative genes.

#_______________________

find0 <- function(x)
{
  loc <- which(x==0)
  return(length(loc)/length(x))
}

PreProcess <- function(data,th=0.98)
{
  # Gene filter
  gSum <- apply(data,1,find0)
  loc <- which(gSum>=th)
  data <- data[-loc,]
  return(data)
}