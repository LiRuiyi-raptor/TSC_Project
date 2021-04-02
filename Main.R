setwd(".../CaFew") 
source("PreProcess.R")
source("CmpRSC.R")
source("Ruiyi_SNN.R")
source("ConNet.R")
source("TSCluster.R")

data <- read.csv(paste0("/Users/ruiyi/Documents/Project/MixMethod/data/",su,"_data.csv"),sep="")

data <- read.csv("...",sep="")
th=0.98
method="scc"
remain = "OneQtr"
data <- PreProcess(data)
RSC <- CmpRSC(data)
if (RSC>=0.8)
{
  data <- log2(data+1)
}
re <- ConNet(data,method,remain)
all_edge <- re$all_edge
temp <- re$temp
pre <- TSCCluster(data,all_edge,temp)
pre_core <- pre$TSC_core # The clustering results for core cells
pre_all <- pre$TSC_all # The clustering results for all cells