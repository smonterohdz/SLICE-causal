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

source("R/slicegm.R") 
source("R/ace_ai.R")
source("R/miscfunctions.R")

set.seed(124)
p <- 5
DAG.true <- synt.mdl(p,pconn = 0.5)
CPDAG.true <- dag2cpdag(DAG.true) 


## generate a set of observations
n.samp <- 1000
## simulate type = c('normal', 'cauchy', 'subsupergaussian') data from the true DAG and the
errM <-disturb(type = 'normal',p = p,nsamp = n.samp, p1 = 0, p2=1)
#data.samp <- rmvDAG(n.samp, DAG.true)
data.samp <- rmvDAG(n.samp, DAG.true,errMat = errM,use.node.names = TRUE)

## shuffle ordr in columns
idx.shuffle <- sample(1:p)
data.samp.sh <- data.samp[,idx.shuffle]

## estimate a SKEL with skeleton (pc algorithm 1st part)
suffStat <- list(C = cor(data.samp), n = n.samp)
skel.est <- skeleton(suffStat, indepTest = gaussCItest, alpha = 0.01, labels = colnames(data.samp), method = "stable")
suffStat.sh <- list(C = cor(data.samp.sh), n = n.samp)
skel.est.sh <- skeleton(suffStat.sh, indepTest = gaussCItest, alpha = 0.01, labels = colnames(data.samp.sh), method = "stable")

## estimate a CPDAG with PC algorithm
suffStat <- list(C = cor(data.samp), n = n.samp)
CPDAG.est <- pc(suffStat, indepTest = gaussCItest, alpha = 0.01, labels = colnames(data.samp), skel.method = "stable")
suffStat.sh <- list(C = cor(data.samp.sh), n = n.samp)
CPDAG.est.sh <- pc(suffStat.sh, indepTest = gaussCItest, alpha = 0.01, labels = colnames(data.samp.sh), skel.method = "stable")

## (NonShuf)
#Determine the causal edge direction from an estimated CPDAG
CPDAG.est.am <- as(CPDAG.est@graph,"matrix")
slice.res <- slicegm(CPDAG.est.am,data.samp) 
shd.rates <-mycomparegraphs3(DAG.true,slice.res$G.est,CPDAG.est)
print(sprintf("From an estimated CPDAG(NonShuf): SHD=%.2f TotEdg=%i",shd.rates$shd.norm,shd.rates$tot.edg))
#Determine the causal edge direction from an estimated skeleton
skel.am <- as(skel.est@graph,"matrix")
slice.res.skel <- slicegm(skel.am,data.samp) 
shd.rates <-mycomparegraphs3(DAG.true,slice.res.skel$G.est,skel.est)
print(sprintf("From an estimated SKEL(NonShuf): SHD=%.2f TotEdg=%i",shd.rates$shd.norm,shd.rates$tot.edg))

## (Shuf) 
# shufling true DAG
DAG.am <- as(DAG.true,"matrix")
DAG.am <- DAG.am[idx.shuffle,idx.shuffle]
DAG.gam <- graphAM(adjMat = DAG.am,edgemode = "directed")
DAG.true.sh <- as(DAG.gam,"graphNEL")

#Determine the causal edge direction from an estimated CPDAG
CPDAG.est.am.sh <- as(CPDAG.est.sh@graph,"matrix")
slice.res.sh <- slicegm(CPDAG.est.am.sh,data.samp.sh)
shd.rates <-mycomparegraphs3(DAG.true.sh,slice.res.sh$G.est,CPDAG.est.sh)
print(sprintf("From an estimated CPDAG(Shuf): SHD=%.2f TotEdg=%i",shd.rates$shd.norm,shd.rates$tot.edg))
#Determine the causal edge direction from an estimated skeleton
skel.am.sh <- as(skel.est.sh@graph,"matrix")
slice.res.sh.skel <- slicegm(skel.am.sh,data.samp.sh) 
shd.rates <-mycomparegraphs3(DAG.true.sh,slice.res.sh.skel$G.est,skel.est.sh)
print(sprintf("From an estimated SKEL(Shuf): SHD=%.2f TotEdg=%i",shd.rates$shd.norm,shd.rates$tot.edg))

par(mfrow=c(4,3))
plot(DAG.true,main="TrueModel(NonShuf)")
plot(CPDAG.est,main="CPDAG(NonShuf)")
plot(slice.res$G.est,main="SLICE(NonShuf)")
plot(DAG.true.sh,main="TrueModel(Shuf)")
plot(CPDAG.est.sh,main="CPDAG(Shuf)")
plot(slice.res.sh$G.est,main="SLICE(Shuf)")
plot(DAG.true.sh,main="TrueModel(Shuf)")
plot(skel.est.sh,main="Skel(Shuf)")
plot(slice.res.sh.skel$G.est,main="SLICE(Shuf)")
plot(DAG.true,main="TrueModel(NonShuf)")
plot(skel.est,main="Skel(NonShuf)")
plot(slice.res.skel$G.est,main="SLICE(NonShuf)")





## Determine the causal edge direction from a given true CPDAG
CPDAG.true.am <- as(CPDAG.true,"matrix")
slice.res <- slicegm(CPDAG.true.am,data.samp.sh) 
par(mfrow=c(1,3))
plot(DAG.true,main="TrueModel")
plot(CPDAG.true,main="CPDAGtrue")
plot(slice.res$G.est,main="SLICEoutcome")
shd.rates <-mycomparegraphs3(DAG.true,slice.res$G.est,CPDAG.true)
print(sprintf("From a true CPDAG: SHD=%.2f TotEdg=%i",shd.rates$shd.norm,shd.rates$tot.edg))

