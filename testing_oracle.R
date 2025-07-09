# The input adjacency matrix must has the format (the edgemark-code refers to the row index)
# 0, 1 for tail and arrowhead 
# 
# amat[a,b]=0 and amat[b,a]=1 implies a --> b
# amat[a,b]=1 and amat[b,a]=0 implies a <-- b
# amat[a,b]=0 and amat[b,a]=0 implies a     b
# amat[a,b]=1 and amat[b,a]=1 implies a <-> b


rm(list = ls())
library(pcalg)
library(Rgraphviz)

source("R/slicegm-oracle.R") 
source("R/ace_ai.R")
source("R/miscfunctions.R")
source("R/randomGuessing.R")

set.seed(143)
p <- 6
DAG.true <- synt.mdl(p,pconn = 0.35)
CPDAG.true <- dag2cpdag(DAG.true) 


## generate a set of observations
n.samp <- 1000
## simulate type = c('normal', 'cauchy', 'subsupergaussian') data from the true DAG and the
errM <-disturb(type = 'cauchy',p = p,nsamp = n.samp, p1 = 0, p2=1)
#data.samp <- rmvDAG(n.samp, DAG.true)
data.samp <- rmvDAG(n.samp, DAG.true,errMat = errM,use.node.names = TRUE)

## generate a causal effects matrix row(interv) column(observ)
ceff.mat <- matrix(0,nrow = p,ncol = p)
for(x in 1:(p-1)){
  for(y in (x+1):p){
    ceff.mat[x,y] <- causalEffect(DAG.true,y,x)
    ceff.mat[y,x] <- causalEffect(DAG.true,x,y)
  }
}


## Determine the causal edge direction from a given true CPDAG
CPDAG.true.am <- as(CPDAG.true,"matrix")
slice.res.orc <- slicegm.oracle(CPDAG.true.am,data.samp,ceff.mat) 
shd.rates <-mycomparegraphs3(DAG.true,slice.res.orc$G.est,CPDAG.true)
print(sprintf("From a true CPDAG(orc): SHD=%.2f TotEdg=%i",shd.rates$shd.norm,shd.rates$tot.edg))

# Determine the causal edge directions by a random guessing coin flip
randGuess.res <- randomGuessing(CPDAG.true.am) 
shd.rates <-mycomparegraphs3(DAG.true,randGuess.res$G.est,CPDAG.true)
print(sprintf("From a true CPDAG(randGuess): SHD=%.2f TotEdg=%i",shd.rates$shd.norm,shd.rates$tot.edg))



##Determine the causal edge direction from an estimated CPDAG
#Estimate a CPDAG with PC algorithm
#suffStat <- list(C = cor(data.samp), n = n.samp)
#CPDAG.est <- pc(suffStat, indepTest = gaussCItest, alpha = 0.01, labels = colnames(data.samp), skel.method = "stable")
#CPDAG.est.am <- as(CPDAG.est@graph,"matrix")
#slice.res.est <- slicegm(CPDAG.est.am,data.samp) 
#shd.rates <-mycomparegraphs3(DAG.true,slice.res.est$G.est,CPDAG.est)
#print(sprintf("From an estimated CPDAG(NonShuf): SHD=%.2f TotEdg=%i",shd.rates$shd.norm,shd.rates$tot.edg))


par(mfrow=c(2,2))
plot(DAG.true,main="TrueModel")
plot(CPDAG.true,main="CPDAG")
plot(randGuess.res$G.est,main="RandomGuessing")
plot(slice.res.orc$G.est,main="SLICE-orc")




