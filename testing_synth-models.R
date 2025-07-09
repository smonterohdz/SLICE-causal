# First run:
# source("inaoe/ACE-IntervComp/code/slicegm/startupSLICE.R",chdir=TRUE)
# source("inaoe/ACE-IntervComp/code/otherCausalDiscoveryMethods/GDSEEV/startups/startupGDS.R",chdir=TRUE)
#rm(list = ls())
library(pcalg)
library(base)


#library(Rgraphviz)
#library(ggplot2)
#library(dplyr)


so <- Sys.info()[['sysname']]
if (so == "Linux") {
  setwd("/home/samuel/Dropbox/inaoe/Reuniones-Seminarios/ConectividadSintetica/Data/")
} else{
  setwd("c:/Users/Gesture/Dropbox/inaoe/Reuniones-Seminarios/ConectividadSintetica/Data/")
}

## Results
path.to.files <- "../results/"
## Output filename format AAAAMMDD-results-EXPERIMENT-PC_XXX
out.file <- paste0(path.to.files,"20190520-results-synthMod-LAP_SAM.RData")
df.results <- data.frame(model=NULL,beta=NULL,lambda=NULL,p=NULL,t.links=NULL,
                         u.links=NULL,nsamp=NULL,shd.val=NULL,shd.norm=NULL,
                         telap=NULL,algo=NULL,errDist=NULL)
## Vars Args to other methods
pars <- list(regr.pars = list())

for(imdl in 1:6){
  #Connectivity network (DAG and CPDAG)
  input.dag <- paste0("model",imdl,"_am.csv")
  DAG.true.am <- as.matrix(read.table(input.dag,sep = ",",header = TRUE))
  rownames(DAG.true.am)<-colnames(DAG.true.am)
  DAG.true <- as(graphAM(adjMat = DAG.true.am, edgemode = "directed",values=list(weight=1)),"graphNEL")
  CPDAG.true <- dag2cpdag(DAG.true)
  #CPDAG.true@graphData$edgemode<-"undirected"
  CPDAG.true.am <-wgtMatrix(CPDAG.true,transpose = FALSE)
  #jpeg(filename = paste0("model",imdl,"_synth.jpg"),bg = "white",width = 800,height = 600)
  #par(mfrow=c(3,4))
  #plot(DAG.true,main="TRUE DAG")
  #plot(CPDAG.true,main="CPDAG")
  
  for(beta in c(0.1,0.5,1.0)){
    for(lambda in c(0.1,0.5,0.9)){
      tim_ch <- paste0("[",format(Sys.time(),format = "%H:%M:%S"),"]")
      cat(tim_ch,"model =",imdl,", beta =",beta, ", lambda = ", lambda, "\r")
      #Synthetic observations
      input.samp <- paste0("model",imdl,"_synth_lambda-",lambda,"_beta-",beta,".csv")
      data.samp <- as.matrix(read.table(input.samp,sep = ",",header = TRUE))
      nsamp <- nrow(data.samp)
      p <- ncol(data.samp)
      t.links <- sum(DAG.true.am)
      #u.links <- nrow(get.undiredges(CPDAG.true.am)$e.list)
      n.samp <- nrow(data.samp)
      
      #==========
      # SLICE
      ptm <- proc.time() 
      slice.res <- slicegm(CPDAG.true.am,data.samp) 
      timeSLICE <- proc.time() - ptm
      my.shd <-mycomparegraphs3(DAG.true,slice.res$G.est,CPDAG.true)
      
      tmp.df <- data.frame(model=paste0("model",imdl),beta=beta,lambda=lambda,
                           p=p,t.links=t.links,u.links=my.shd$tot.edg,nsamp=nsamp,
                           shd.val=my.shd$shd.val,shd.norm=my.shd$shd.norm,
                           telap=timeSLICE['elapsed'],algo="SLICE",errDist="fnirs")
      df.results <- rbind(df.results,tmp.df)
      #==========
      # LINGAM
      ptm <- proc.time() 
      resLINGAM <- 1L*(lingam(data.samp)$Bpruned!=0)
      timeLINGAM <- proc.time() - ptm
    
      resLINGAM.gNEL <- as(graphAM(adjMat=resLINGAM,edgemode="directed"),"graphNEL")
      my.shd <- mycomparegraphs3(DAG.true,resLINGAM.gNEL,CPDAG.true)
  
      tmp.df <- data.frame(model=paste0("model",imdl),beta=beta,lambda=lambda,
                           p=p,t.links=t.links,u.links=my.shd$tot.edg,nsamp=nsamp,
                           shd.val=my.shd$shd.val,shd.norm=my.shd$shd.norm,
                           telap=timeLINGAM['elapsed'],algo="LiNGAM",errDist="fnirs")
      df.results <- rbind(df.results,tmp.df)
  
      #==========
      # GDS
      ptm <- proc.time() 
      resGDS <- GDS(data.samp, scoreName = "SEMSEV", pars, check = "checkUntilFirstMinK", output = FALSE,startAt = "emptyGraph")$Adj
      timeGDS <- proc.time() - ptm
      
      resGDS.gNEL <- as(graphAM(adjMat=resGDS,edgemode="directed"),"graphNEL")
      my.shd <- mycomparegraphs3(DAG.true,resGDS.gNEL,CPDAG.true)
      
      tmp.df <- data.frame(model=paste0("model",imdl),beta=beta,lambda=lambda,
                           p=p,t.links=t.links,u.links=my.shd$tot.edg,nsamp=nsamp,
                           shd.val=my.shd$shd.val,shd.norm=my.shd$shd.norm,
                           telap=timeGDS['elapsed'],algo="GDS",errDist="fnirs")
      df.results <- rbind(df.results,tmp.df)
    }
  }
  #dev.off()
}
save.image(out.file)