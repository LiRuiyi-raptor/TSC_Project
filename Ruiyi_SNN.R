# Created by Ruiyi Li on January 13, 2021
# This function return the similarity between nodes based on SNN(Shared Neighbor hoods)
#______________________
#Input:
# inputTags: a numberic matrix with rows are genes/transcripts and colums are cells
# K: The number of neighbores used to compute the similarity

#Output
# Edge: similarity between each two nodes(V1=node1,V2=node2,V3=SNN similarity)
#_______________________

mtrx2cols = function(simM){
  lt = lower.tri(simM)  #获取下半角为TRUE，上半角为FALSE的同维度矩阵；
  res = data.frame(row = row(simM,as.factor = T)[lt],  #返回矩阵m1的下半角元素对应的行名
                   col = col(simM,as.factor = T)[lt],  #返回矩阵m1的下半角对应的列名
                   val1 = simM[lt]) #按列依次获取矩阵下半角的元素
  names(res)[1:3] = c("V1","V2","V3") #对后两列重命名，支持多个矩阵合并
  return(res)
}

SNNew <- function(inputTags,k)
{
  inputTags <- t(as.matrix(inputTags))
  numSpl<-dim(inputTags)[1]
  m <- dist(inputTags, diag=TRUE, upper=TRUE)
  x<-as.matrix(m)
  IDX<-t(apply(x,1,order)[1:k,]) # knn list
  simM <- matrix(0,numSpl,numSpl)
  for (i in 1:numSpl)
  {
    interF <- function(x)
    {
      strength <- 0
      shared <- intersect(x, IDX[i,])
      if(length(shared)>0)
      {			
        s<-k-0.5*(match(shared, IDX[i,])+match(shared, x))
        strength<-max(s)
      } 
      return(strength)
    }
    simM[i,] <- apply(IDX,1,interF)
  }
  Edge <- mtrx2cols(simM)
  Edge <- Edge[-which(Edge$V3==0),]
  return(Edge)
}
