#' Structure Learning via Intervals of Causal Effects (SLICE)
#'
#' Resolve undefined causal edges froma an equivalence class 
#' represented by a completed partially directed acyclic graphCPDAG
#' 
#' @param G.am An adjacency matrix of the initial CPDAG containig one or more undefined causal edges to be resolved. The edgemark-code in am refers to the row index.
#' @param V.data A matrix with \code{"m"} observations (rows) of the \code{"p"} variables (columns)
#'
#' @return A data frame with the learnt structure (\code{"G.est"}), the type of outcome (\code{"G.outcome"}), 
#' a data frame with information regarding the computed intervals of causal effects (\code{"G.df_ice"}),
#' the number of iteration elapsed (\code{"iters"}).
#' @author Samuel Montero-Hernandez, \email{samuel@@inaoep.mx}
#' @keywords graphical models
#'
#' @examples
#' source("slicegm.R") 
#' source("ace_ai-v2.R")
#' source("miscfunctions.R")
#' 
#' set.seed(101)
#' p <- 7
#' DAG.true <- randomDAG(p, prob = 0.5) 
#' CPDAG.true <- dag2cpdag(DAG.true) 
#' 
#' ## generate a set of observations
#' n.samp <- 10000
#' data.samp <- rmvDAG(n.samp, DAG.true)
#' 
#' ## estimate a CPDAG with PC algorithm
#' suffStat <- list(C = cor(data.samp), n = n.samp)
#' CPDAG.est <- pc(suffStat, indepTest = gaussCItest, alpha = 0.01, p=p)
#' 
#' par(mfrow=c(1,3))
#' plot(DAG.true,main="TrueModel")
#' plot(CPDAG.true,main="CPDAGtrue")
#' plot(CPDAG.est,main="CPDAGestim")
#' 
#' 
#' ## Determine the causal edge direction from an estimated CPDAG
#' CPDAG.est.am <- as(CPDAG.est@@graph,"matrix")
#' slice.res <- slicegm(CPDAG.est.am,data.samp) 
#' par(mfrow=c(1,3))
#' plot(DAG.true,main="TrueModel")
#' plot(CPDAG.est,main="CPDAGestim")
#' plot(slice.res$G.est,main="SLICEoutcome")
#' shd.rates <-mycomparegraphs3(DAG.true,slice.res$G.est,CPDAG.est)
#' print(sprintf("From an estimated CPDAG: SHD=%.2f",shd.rates$shd.norm))
#' 
#' ## Determine the causal edge direction from a given true CPDAG
#' CPDAG.true.am <- as(CPDAG.true,"matrix")
#' slice.res <- slicegm(CPDAG.true.am,data.samp) 
#' par(mfrow=c(1,3))
#' plot(DAG.true,main="TrueModel")
#' plot(CPDAG.true,main="CPDAGtrue")
#' plot(slice.res$G.est,main="SLICEoutcome")
#' shd.rates <-mycomparegraphs3(DAG.true,slice.res$G.est,CPDAG.true)
#' print(sprintf("From a true CPDAG: SHD=%.2f",shd.rates$shd.norm))
#' 
#' @export
#' 




randomGuessing <- function(G.am) {
  amGM <- graphAM(adjMat = G.am,edgemode = "directed")
  G.cpdag <- as(amGM,"graphNEL")
  cpdagFlag <- FALSE
  beta <- (-1)
  iters <- 0
  
  # covariance matrix from data
  
  ## list of undefined edges
  undir.res <- get.undiredges(G.am)
  e.list <- undir.res$e.list
  nr.eUndir <- nrow(e.list)
  nodes.lab <- G.cpdag@nodes
  G.est <- G.cpdag
  while (nr.eUndir > 0 & cpdagFlag == FALSE) {
    iters <- iters + 1
    #df.iceL <- df.null
    
    for (i in 1:nr.eUndir) {
      
      x <- e.list[i, 1]
      y <- e.list[i, 2]
	  ## Random orientation
      flip <- abs(rnorm(1,0,1))
      
      if(flip>=0.5){
        G.est <-removeEdge(nodes.lab[y],as.character(nodes.lab[x]),G.est)
        G.est@graphData$edgemode<-"directed"
        G.est <-addEdge(as.character(nodes.lab[x]),as.character(nodes.lab[y]),G.est)
      }else{
        G.est <-removeEdge(nodes.lab[x],as.character(nodes.lab[y]),G.est)
        G.est@graphData$edgemode<-"directed"
        G.est <-addEdge(as.character(nodes.lab[y]),as.character(nodes.lab[x]),G.est)
      }
      
    }    
    ## New undefined edges
    undir.res <- get.undiredges(as(G.est,"matrix"))
    e.list <- undir.res$e.list
    nr.eUndir <- nrow(e.list)
    cpdagFlag <- FALSE#ornt.res$cpdagFlag
    #browser()
  }# en while n.eundir >0 AND cpdagFlag==FALSE
  if (isDirected(G.est) == TRUE) {
    outcome <- "DAG"
  } else{
    outcome <- "CPDAG"
  }

  return(list(
    G.est = G.est,
    G.outcome = outcome#,
    #G.df_ice = G.df_ice,
    #iters = iters
  ))
}
