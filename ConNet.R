# Created by Ruiyi Li on Apri 2, 2021
# This function return the graph created by different methods
#______________________
#Input:
# data: a numberic matrix with rows are genes/transcripts and colums are cells
# method：choose one method to cmpute the similiraty between cells: pcc (Pearson correlation coefficient), scc (Spearman correlation coefficient), ed (Euclidean distance), md (Manhattan distance),snn (shared nearest neighbors)
# remain: how many edges to be remain: OneQtr,TwoQtr,ThreeQtr,All

#Output
# all_edge: The similarity between each two cell
# temp : The similarity between each two core cell
#_______________________


source("Ruiyi_SNN.R")
ConNet <- function(data,method,remain)
{
  if (method=="ed")
  {
    ed <- dist(t(as.matrix(data))) 
    ed <- as.matrix(ed)
    g <- graph_from_adjacency_matrix(ed,weighted=TRUE,mode="undirected",diag=F)
    edge_ed <- get.data.frame(g)
    edge_sim <- 1-(edge_ed$weight/max(edge_ed$weight))
    ed.pvalue = pnorm(-abs(scale(edge_sim, center=T,scale=T)))
    all_edge <- data.frame(from=edge_ed$from,to=edge_ed$to,weight=edge_sim,pvalue=ed.pvalue,distance=edge_ed$weight)
  }else if(method=="md"){
    md <- dist(t(as.matrix(data)),method="manhattan") 
    md <- as.matrix(md)
    g <- graph_from_adjacency_matrix(md,weighted=TRUE,mode="undirected",diag=F)
    edge_md <- get.data.frame(g)
    edge_sim <- 1-(edge_md$weight/max(edge_md$weight))
    md.pvalue = pnorm(-abs(scale(edge_sim, center=T,scale=T)))
    all_edge <- data.frame(from=edge_md$from,to=edge_md$to,weight=edge_sim,pvalue=md.pvalue,distance=edge_md$weight)
  }else if (method=="pcc"){
    pcc <- cor(as.matrix(data))
    pcc <- as.matrix(pcc)
    g <- graph_from_adjacency_matrix(pcc,weighted=TRUE,mode="undirected",diag=F)
    edge_pc <- get.data.frame(g)
    if (length(which(edge_pc$weight<0))>0)
    {
      edge_pc <- edge_pc[-(which(edge_pc$weight<0)),]
    }
    pc.pvalue = pnorm(-abs(scale(edge_pc$weight, center=T,scale=T)))
    all_edge <- data.frame(from=edge_pc$from,to=edge_pc$to,weight=edge_pc$weight,pvalue=pc.pvalue)
  }else if(method=="scc"){
    scc <- cor(as.matrix(data),method="spearman")
    scc <- as.matrix(scc)
    g <- graph_from_adjacency_matrix(scc,weighted=TRUE,mode="undirected",diag=F)
    edge_sc <- get.data.frame(g)
    loc <- which(edge_sc$weight<0)
    if (length(loc)>0)
    {
      edge_sc <- edge_sc[-loc,]
    }
    sc.pvalue = pnorm(-abs(scale(edge_sc$weight, center=T,scale=T)))
    all_edge <- data.frame(from=edge_sc$from,to=edge_sc$to,weight=edge_sc$weight,pvalue=sc.pvalue)
  }else {
    re <- SNNew(inputTags=data,k=5)
    map <- data.frame(cellID=c(1:ncol(data)),cellName=colnames(data))
    re <- merge(re,map,by.x="V1",by.y="cellID")
    re <- merge(re,map,by.x="V2",by.y="cellID")
    all_edge <- data.frame(from=re$cellName.x,to=re$cellName.y,weight=re$V3)
  }
  
  # Delete Edges
  if(method!="snn")
  {
    lx <- summary(all_edge$weight)
    if (remain=="OneQtr")
    {
      temp <- all_edge[all_edge$weight>=lx[5],]
    }else if(remain=="TwoQtr"){
      temp <- all_edge[all_edge$weight>=lx[3],]
    }else if(remain=="ThreeQtr"){
      temp <- all_edge[all_edge$weight>=lx[2],]
    }else if(remain=="All"){
      temp <- all_edge
    }
    
  }else{
    temp <- all_edge
  }
  
  return(list(all_edge=all_edge,temp=temp))
}