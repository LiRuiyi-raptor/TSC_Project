# Created by Ruiyi Li on April 2, 2020
# This function return the right skewed coefficient (RSC) of data.
#______________________
#Input:
# data: a numberic matrix with the rows are genes and colums are cells

#Output
# RSC: a numeric value

#_______________________

CmpRSC <- function(data)
{
  staRe <- apply(data,1,max)
  staRe <- as.matrix(staRe[order(staRe,decreasing = T)])
  re <- as.matrix(quantile(staRe))
  Q1 <- re[2]
  Q3 <- re[4]
  IQR <- Q3-Q1
  th1 <- Q1-1.5*IQR
  th2 <- Q3+1.5*IQR
  loc <- intersect(which(staRe>=th1),which(staRe<=th2))
  length(loc)/nrow(staRe)
  staRe <- staRe[loc]
  meaV <- mean(staRe)
  staReright <- staRe[which(staRe>=meaV)]
  staReright <- staReright-meaV
  RSC <- mean(staReright)/meaV
  return(RSC)
}
