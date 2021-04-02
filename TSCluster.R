# Created by Ruiyi Li on Apri 2, 2021
# This function return the predicting results of TSC
#______________________
#Input:
# data: a numberic matrix with rows are genes/transcripts and colums are cells
# all_edge: The similarity between each two cell
# temp: The similarity between each two core cell

#Output
# re_TSC_all: clustering outcomes with two columns that one indicates cell names and the other is the predicted cell labels.
#_______________________


library("igraph")
TSCCluster <- function(data,all_edge,temp)
{
  isoNum <- ncol(data)
  g <- graph_from_data_frame(temp,directed=F)
  while (isoNum>0)
  {
    bad.vs = V(g)[degree(g) == 0]
    g = delete.vertices(g, bad.vs)
    temp <- get.data.frame(g)
    re <- walktrap.community(g,weights=E(g)$weight,step=4)
    pre_l <- as.matrix(re$membership)
    re_temp <- table(pre_l)
    loc <- which(re_temp<=1)
    if(length(loc)>0)
    {
      isoNum <- length(loc)
      isoName <- re$names[match(loc,re$membership)]
      g = delete.vertices(g, isoName)
    }else{
      isoNum=0
    }
  }
  pre_all <- matrix(0,ncol(data),1)
  k <- length(unique(pre_l))
  center <- matrix(0,nrow(data),k)
  for (i in 1:k)
  {
    tempCell <- re$names[which(re$membership==i)]
    loc <- match(tempCell,colnames(data))
    tempData <- data[,loc]
    if(is.null(ncol(tempData)))
    {
      center[,i] <- tempData
    }else{
      center[,i] <- apply(tempData,1,mean)
    }
  }
  
  if (ncol(data)>length(pre_l))
  {
    cells <- union(temp$from,temp$to)
    loc <- match(cells,colnames(data))
    pre_all[loc,1] <- pre_l
    
    leftCell <- as.matrix(colnames(data)[-loc])
    ln <- nrow(leftCell)
    pre_l_add <- matrix(0,ln,1)
    for (i in 1:ln)
    {
      tempNode <- data[,which(colnames(data)==leftCell[i])]
      tempNode <- rbind(tempNode,t(center))
      disE <- as.matrix(cor(t(tempNode)))
      diag(disE) <- 0
      pre_l_add[i,1] <- which(disE[1,]==max(disE[1,]))-1
    }
    loc1 <- match(leftCell,colnames(data))
    pre_all[loc1,1] <- pre_l_add
  }else{
    pre_all <- pre_l
  }
  re_TSC_all <- data.frame(cell_name=colnames(data),cell_label=pre_all)
  re_TSC_core <- data.frame(cell_name=re$names,cell_label=pre_l)
  re_TSC <- list(TSC_all=re_TSC_all,TSC_core=re_TSC_core)
  return(re_TSC)
}