rm(list = ls())
library(pcalg)
library(Rgraphviz)

source("R/slicegm.R") 
source("R/ace_ai.R")
source("R/miscfunctions.R")

data.samp <- read.table("../../Data/chalearn/pair0008.txt",sep = "")
colnames(data.samp) <- c("X","Y")
p <- ncol(data.samp)
n.samp <- nrow(data.samp)

DAG.true.am <- matrix(data = c(0,1,0,0),nrow = 2, byrow = TRUE)
colnames(DAG.true.am)<-c("X", "Y")
DAG.true <- graphAM(adjMat = DAG.true.am, edgemode = "directed")

CPDAG.est.am <- matrix(data = c(0,1,1,0),nrow = 2,  byrow = TRUE)
colnames(CPDAG.est.am)<-c("X", "Y")
CPDAG.est <- graphAM(adjMat = CPDAG.est.am, edgemode = "directed")

#Determine the causal edge direction from an estimated CPDAG
slice.res <- slicegm(CPDAG.est.am,data.samp) 
shd.rates <-mycomparegraphs3(DAG.true,slice.res$G.est,CPDAG.est)
print(sprintf("From an estimated CPDAG: SHD=%.2f TotEdg=%i",shd.rates$shd.norm,shd.rates$tot.edg))
par(mfcol=c(1,2))
plot(DAG.true)
plot(slice.res$G.est)