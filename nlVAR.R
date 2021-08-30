#source("/home/students/afshinh/Thesis/2014_11_12/nlVAR_util.R")
source("C:/Users/Afshin/Documents/Ali/Thesis/ThLasso/2014_11_21/nlVAR_util.R")

#~~~~~~~~~~~~plot.cv~~~~~~~~~~~~
#This funciton plots the cv.profile
plot.tune <- function(tune.profile){
  
  avg <- tune.profile$pen
  sdev <- 0
  if(!is.null(tune.profile$pen.se))
    sdev <- tune.profile$pen.se
  x <- log(tune.profile$parm)
  parm.min <- tune.profile$parm.min
  
  title <- paste("Tuning by", tune.profile$tune.method)
  ylab <- paste("Loss: ", tune.profile$pen.name)
  xlab <- paste("log (",tune.profile$parm.name,")")
  rng <- range(c(avg-sdev, avg+sdev))
  par (mfrow=c(1,1))
  plot(x, avg,
       ylim=rng,
       pch=19, xlab=xlab, ylab=ylab,
       main=title)
  if(!is.null(tune.profile$pen.se))
    arrows(x, avg-sdev, x, avg+sdev, length=0.05, angle=90, code=3)
  arrows(x0=log(parm.min),y0=rng[1],x1=log(parm.min),y1=rng[2],angle=0,lty=2,col=2)
  if(!is.null(tune.profile$parm.min.1se)){
    parm.min.1se <- tune.profile$parm.min.1se
    arrows(x0=log(parm.min.1se),y0=rng[1],x1=log(parm.min.1se),y1=rng[2],angle=0,lty=2,col=3)
  }
}

#~~~~~~~~~~~~~~~~~nlVAR.graph~~~~~~~~~~~~~~~~~~~~~~~~
# given an adjacency matrix or Var.fit, it generates
# an i.graph

nlVAR.graph <- function(g,main=NULL,diag = TRUE, arrow.sz = .4){
  require(igraph)
  if(!is.matrix(g))
    g <- g$est.network
  p <- ncol(g)
  node.names <- dimnames(g)[[1]]
  if(is.null(node.names))
    node.names <- paste("p",1:p,sep="")
  dimnames(g) <- list(node.names,node.names)
  plot(graph.adjacency(adjmatrix=g,mode="directed",diag=diag)
       , main = main,edge.arrow.size = arrow.sz,layout=layout.circle)
}
#~~~~~~~~~~~~~~~~~~~~~~filter.out.edge~~~~~~~~~~~~~~~~~~
filter.out.edge <- function(X.npt, trt.target.map,sig.level = .3,link = log, ctrl.name = NULL){
  exclusion.adj <- get.exclusion.adj(X.npt = X.npt, trt.target.map = trt.target.map ,sig.level = sig.level,link = link, ctrl.name = ctrl.name)
  return(exclusion.adj)
}
#~~~~~~~~~~~~~~~~~~~~~~~eval.network~~~~~~~~~~~~~~~~~~~~
# Given the true and estimated logical adjacency matrices
# it finds the TP, FP, TN, and FN
eval.network <- function(adj.true, adj.est,exclude.diag = F){
  if(!all.equal(dim(adj.true),dim(adj.est))) stop("dimensions of adj.true and adj.est don't match")
  if(!is.logical(adj.true)|| !is.logical(adj.est)) stop("adj.true and adj.est don't match need to be logical arrays")
  if(exclude.diag){
    diag(adj.true) <- F
    diag(adj.est) <- F
  }
  TP <- sum(adj.est[adj.true])
  FP <- sum(adj.est[!adj.true])
  TN <- sum(!adj.est[!adj.true])
  FN <- sum(!adj.est[adj.true])
  return(list(TP=TP,FP=FP,TN=TN,FN=FN))
}
#~~~~~~~~~~~~~~~~~~~~~~nlVAR~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Given an n x p x t data array X, it separates X and Y and expands X.
# It then regresses each of the p columns of Y on the X and computes
# p columns of length npt of the estimates in a npt x p matrix format

nlVAR <- function(X, lambda,trt.target.map = NULL ,smoother = c("nspline","bspline","polynomial"),exclude.adj = NULL, smooth.deg=3, intercept = TRUE, get.adj = F, standardize = F){
  if(is.null(trt.target.map))
    trt.target.map <- list(none = "")
  smoother <- match.arg(smoother)
  dm <- dim(X)
  if(length(dm)!=3) stop("An array of nxpxT is expected")
  n <- dm[1]
  p <- dm[2]
  t <- dm[3]
  X <- scale.XArray(X,n,t)
  dat <- get.XY(X,n,p,t,smooth.deg,smoother = smoother)
  #~~~~get penalty factor based onExclude edge~~
  if(is.null(exclude.adj))
    exclude.adj <- matrix(FALSE,p,p)
  #pfs <- get.penFactors(exclude.adj,t.num = t)
  grp <- rep(1:(p*(t-1)),each = smooth.deg)
  pfs <- get.penFactors(exclude.adj,grp)
  #~~~~~~~~~~~~get.include.obs.mat~~~~~~~~~~~~~~
  include.obs.mat <- get.include.obs.mat(Y.np = dat$Y, trt.target.map = trt.target.map)
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  # Construct standardization parameters
  #-------------------------------------
  stand <- NULL
  if(standardize){
    stand <- get.stand.param(X=dat$XXpanded,smooth.deg = smooth.deg)
  }
  nlVAR.fit <- VAR(X=dat$XXpanded, Y=dat$Y, pfs = pfs,lambda=lambda, smooth.deg=smooth.deg, intercept = intercept, get.adj = get.adj,include.obs.mat = include.obs.mat, stand = stand)

}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~nlVAR.tune~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# This function performs cross validation by calling nlVAR

nlVAR.tune <- function(X
                       ,trt.target.map = NULL
                       ,lmb.range = c(.005,1)
                       ,lmbs.len = 15
                       , cRatio.range = c(.5,100)
                       ,cRatios.len = 20
                       ,gamma.const = .1
                       ,exclude.adj = NULL
                       , smoother = c("nspline","bspline","polynomial")
                       ,tune.method = c("BIC","CV")
                       ,kfold = NULL
                       ,smooth.deg = 3
                       , intercept = TRUE
                       , beta.err = .1
                       ,standardize = F
                       , threshold.method = c("uniform","decaying")
                       ,verbose = F){
  tune.method <- match.arg(tune.method)
  smoother <- match.arg(smoother)
  if(is.null(trt.target.map))
    trt.target.map <- list(none = "")
  if(tune.method == "CV" & is.null(kfold)){
    kfold <- 5
  } 
  dm <- dim(X)
  if(length(dm)!=3) stop("An array of nxpxT is expected")
  n <- dm[1]
  p <- dm[2]
  t <- dm[3]
  #~~~~~~~~~~Scale X and expand~~~~~~
  X <- scale.XArray(X,n,t)
  XY <- get.XY(X,n,p,t,smooth.deg,smoother=smoother)
  XXpanded <- XY$X
  Y <- XY$Y
  #~~~~get penalty factor based onExclude edge~~
  if(is.null(exclude.adj))
    exclude.adj <- matrix(FALSE,p,p)
  #pfs <- get.penFactors(exclude.adj,t.num = t)
  grp <- rep(1:(p*(t-1)),each = smooth.deg)
  pfs <- get.penFactors(exclude.adj,grp)
  
  # Construct standardization parameters
  #-------------------------------------
  stand <- NULL
  if(standardize){
    stand <- get.stand.param(X=XXpanded,smooth.deg = smooth.deg)
  }
  #~~~~~~~~~~~~~~~~~~get.include.obs.mat~~~~~~~~~~~~~~~
  #include.obs.mat <- get.include.obs.mat(Y.np = Y, trt.target.map = trt.target.map)
  #~~~~~~~~~~~~~~~~~~get.variance~~~~~~~~~~~~~~~~~~~~~~
  #sigma2.hat <- get.variance(X.npt=X,include.obs.mat)
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if(verbose)
    cat(smoother, "degree =", smooth.deg,"\n")
  
  lmbs <- sort(exp(seq(from=log(lmb.range[1]),to=log(lmb.range[2]),length.out=lmbs.len)),decreasing = T)
  cRatios <- sort(c(exp(seq(from=log(cRatio.range[1]),to=log(cRatio.range[2]),length.out=cRatios.len)),0),decreasing = T)
  
  #~~~~~~~~~~~~~~~~~~~BIC~~~~~~~~~~~~~~~~~~~~~~~~~
  if(tune.method == "BIC"){
    BIC.tune <- tune.BIC(XXpanded = XXpanded
                         ,Y = Y
                         ,trt.target.map = trt.target.map
                         ,pfs = pfs
                         ,lmbs = lmbs
                         ,cRatios = cRatios
                         ,smooth.deg = smooth.deg
                         ,intercept = intercept
                         ,beta.err = beta.err
                         ,threshold.method = threshold.method
                         ,gamma.const = gamma.const
                         ,stand = stand)
    
    
    tune.profile <- BIC.tune
  }else{
    #~~~~~~~~~~~~~~~~~~~CV~~~~~~~~~~~~~~~~~~~~~~~~
    CV.tune <- tune.CV(XXpanded = XXpanded
                       ,Y = Y
                       ,trt.target.map = trt.target.map
                       ,pfs = pfs
                       ,lmbs = lmbs
                       ,cRatios = cRatios
                       ,kfold = kfold
                       ,smooth.deg = smooth.deg
                       ,intercept = intercept
                       ,threshold.method = threshold.method)
    tune.profile <- CV.tune
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     
  }
  
 return(tune.profile)
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~Th.nlVAR~~~~~~~~~~~~~~~~~~~~~~~
# This function is the main threshold non-linear VAR
# it calls the nlVAR.tune, using this best lambda it then
# does cross varlidation to find the best threhold
# It then thresholds the estimate and returns

Th.nlVAR <- function(X
                     ,trt.target.map = NULL
                     ,lambda.n = NULL
                     ,cRatio.n = NULL
                     ,exclude.adj = NULL
                     ,lmb.range = c(.005,1)
                     ,lmbs.len=15
                     ,cRatio.range = c(.5,100)
                     ,cRatios.len = 20
                     ,gamma.const = 0.1
                     ,smoother = c("nspline","bspline","polynomial")
                     ,smooth.deg=3
                     ,kfold=5
                     ,tune.method = c("BIC","CV")
                     #,parm.pick = c("min","min.1se")
                     ,intercept= TRUE
                     ,return.BIC = FALSE
                     ,verbose=FALSE
                     ,beta.err = .1
                     #,cRatio = NULL
                     ,threshold.method = c("uniform","decaying")
                     , standardize = F){ # note only uniform has theoretical proof 
 
  smoother <- match.arg(smoother)
  threshold.method <- match.arg(threshold.method)
  dm <- dim(X)
  if(length(dm)!=3) stop("An array of nxpxT is expected")
  n <- dm[1]
  p <- dm[2]
  t <- dm[3]
  if(!is.numeric(lmb.range) | length(lmb.range)!=2)
    stop("lmb.range needs to be a numeric of length specifying the start and end of the penalty term range")
  if(!is.numeric(cRatio.range) | length(cRatio.range)!=2)
    stop("cRatio.range needs to be a numeric of length specifying the start and end of the cRatio range")
  
  if(is.null(trt.target.map))
    trt.target.map <- list(none = "")
  
  tune.profile <- NULL
  if(is.null(lambda.n) | is.null(cRatio.n)){
    tune.method <- match.arg(tune.method)
    #parm.pick <- match.arg(parm.pick)
    tune.profile <- nlVAR.tune(X=X
                               ,trt.target.map = trt.target.map
                               ,lmb.range = lmb.range
                               ,lmbs.len=lmbs.len
                               ,cRatio.range = cRatio.range
                               ,cRatios.len = cRatios.len
                               ,gamma.const = gamma.const
                               ,smoother = smoother
                               ,exclude.adj = exclude.adj
                               ,tune.method= tune.method
                               ,kfold=kfold
                               ,smooth.deg=smooth.deg
                               ,intercept= intercept
                               ,beta.err = beta.err
                               ,threshold.method = threshold.method
                               ,standardize = standardize
                               ,verbose=verbose) 
    
    lambda.n <- tune.profile$lambda.min
    cRatio.n <- tune.profile$cRatio.min
    
    #if(tune.method == "CV")
      #lambda.n <- c(tuneLmb.profile$parm.min,tuneLmb.profile$parm.min.1se)[c("min","min.1se") == parm.pick]
  }
  
  #FITTING THE WITH lambda.n
  nlVAR.fit <- nlVAR(X=X,trt.target.map=trt.target.map ,lambda=lambda.n,smoother = smoother ,exclude.adj = exclude.adj,smooth.deg = smooth.deg,intercept = intercept, get.adj = T, standardize = standardize)
  
  b0 <- NULL
  if(intercept)
    b0<- nlVAR.fit$beta.mat[1,]

  #Thresholding
  #tunecRatio.profile <- NULL
  tau <- cRatio.n * lambda.n
  
  
  #+++++++++
  #tmp0 <- nlVAR.fit
  #tmp1 <- tune.profile
  #tmp2 <- tau
  #tmp3 <- as.numeric(cRatio.n!=0) * beta.err
  #tmp4 <- threshold.method
  #tmp5 <- cRatio.n
  #tmp6 <- lambda.n
  #save(list = c("tmp0","tmp1","tmp2","tmp3","tmp4","tmp5", "tmp6"),
       #file = "C:/Users/Afshin/Documents/Ali/Thesis/ThLasso/2014_11_21/output_dump/debug.RData")
  #+++++++++
  
  Th.adj <- threshold(A.pptd=nlVAR.fit$adj.mat, tau=tau, beta.err=as.numeric(cRatio.n!=0) * beta.err, threshold.method = threshold.method )$Th.A
  
  # refit
  grp <- rep(1:(p*(t-1)),each=smooth.deg)
  Th.coef <- Adj.to.coef(adj=Th.adj,beta0=b0,smooth.deg=smooth.deg)
  include.cf.mat <- Th.coef[1:ncol(nlVAR.fit$input$XX)+ intercept,] !=0
  pfs.th <-sqrt(smooth.deg)*(as.matrix(aggregate(.~grp,data = data.frame(grp=grp,include.cf.mat),FUN = sum)[,-1])!=0)
  pfs.th[pfs.th==0] <- 10000*sqrt(smooth.deg)
  
  refit <- VAR(X = nlVAR.fit$input$XX
               ,Y = nlVAR.fit$input$YY
               ,lambda = lambda.n
               ,pfs = pfs.th
               ,smooth.deg = smooth.deg
               ,intercept = intercept
               ,get.adj = T
               ,include.obs.mat = nlVAR.fit$input$include.obs.mat
               ,stand = nlVAR.fit$input$stand)

  #cat("gamma const", gamma.const,"\n")
  # Compute BIC
  BIC <- list()
  if(return.BIC){
    BIC.parm <- list(lambda = lambda.n,
                     gamma.const = gamma.const,
                     grp = grp)
    #print(BIC.parm$gamma.const)
    
    BIC$fit <- mean(get.cost(Y=nlVAR.fit$input$YY
                               ,X = nlVAR.fit$input$XX
                               ,betas = nlVAR.fit$beta.mat
                               ,intercept = intercept
                               ,include.obs.mat = nlVAR.fit$input$include.obs.mat
                               ,include.cf.mat = nlVAR.fit$input$include.cf.mat
                               ,cost.type = "BIC"
                               ,BIC.parm = BIC.parm
                               ,SGL = standardize))
    
    BIC$refit <- mean(get.cost(Y=nlVAR.fit$input$YY
                         ,X = nlVAR.fit$input$XX
                         ,betas = refit$beta.mat
                         ,intercept = intercept
                         ,include.obs.mat = refit$input$include.obs.mat
                         ,include.cf.mat = refit$input$include.cf.mat
                         ,cost.type = "BIC"
                         ,BIC.parm = BIC.parm
                         ,SGL = standardize))
  }
  
  
  est.network <- apply(refit$adj.mat !=0,MARGIN = 1:2,FUN = any)
  dimnames(est.network) <- list(dimnames(X)[[2]], dimnames(X)[[2]])
  return(list(Th.adj = Th.adj,
              nlVAR.fit =  nlVAR.fit,
              nlVAR.refit = refit,
              tune.profile = tune.profile,
              est.network = est.network,
              tau = tau,
              BIC = BIC))
}


#~~~~~~~~~~~~~~~~~~~~~~BMA~~~~~~~~~~~~~~~~~~~~~~~~~
# This function receiving an optimial lambda, lambda.n
# perform Bayesian Model Averaging by calling Th.nlVAR

BMA.nlVAR <- function(X
                      ,trt.target.map = NULL
                      ,lambda.n = NULL
                      ,cRatio.n = NULL
                      ,exclude.adj = NULL
                      ,lmb.range = c(0.005,1)
                      ,lmbs.len = 15
                      ,cRatio.range = c(0.5,100)
                      ,cRatios.len = 20
                      ,gamma.const = 0.1
                      ,BMA.threshold = .8
                      ,smoother = c("nspline","bspline","polynomial")
                      ,smooth.deg=3
                      ,intercept = TRUE
                      ,beta.err = .1
                      ,standardize = F){
  
  smoother <- match.arg(smoother)
  dm <- dim(X)
  if(length(dm)!=3) stop("An array of nxpxT is expected")
  n <- dm[1]
  p <- dm[2]
  t <- dm[3]
  if(!is.numeric(lambda.n)| length(lambda.n)!=1)
    stop("lambda.n needs to be signle positive value")
  
  if(BMA.threshold<0 | BMA.threshold>1)
    stop("BMA.threshold is in [0 1]")
  
  if(is.null(trt.target.map))
    trt.target.map <- list(none = "")
  
  
  t.nets <- array(NA,dim = c(p,p,t - 1),
                  dimnames = list( dimnames(X)[[2]],dimnames(X)[[2]], paste(1:(t -1)))   )
  
  t.BICs <- rep(NA, t-1)
  
  for (i in 2:t){
    
    fit<-Th.nlVAR(X = X[,,1:i]
                  ,trt.target.map = trt.target.map
                  ,lambda.n = lambda.n
                  ,cRatio.n = cRatio.n
                  ,exclude.adj = exclude.adj
                  ,lmb.range = lmb.range
                  ,lmbs.len = lmbs.len
                  ,cRatio.range = cRatio.range
                  ,cRatios.len = cRatios.len
                  ,gamma.const = gamma.const
                  ,tune.method = "BIC"
                  ,smoother = smoother
                  ,smooth.deg = smooth.deg
                  ,intercept = intercept
                  ,beta.err = beta.err
                  ,threshold.method = "uniform"
                  ,standardize = standardize)
    
    
    t.BICs[i-1] <- min(fit$tune.profile$pen)[1]
    t.nets[,,i-1] <- fit$est.network
    #nlVAR.graph(trueNetwork ==1,diag = F, main = "True")
    #nlVAR.graph(fit$est.network,diag = F, main = paste("BIC =",round(t.BICs[i-1],2) )) 
  }
  
  M <- max(-t.BICs/2)
  log.W <- -t.BICs/2 - (M + log(sum(exp(-t.BICs/2 - M) )) )
  
  W.agg.net <- matrix(0, p, p,dimnames = list(dimnames(X)[[2]],dimnames(X)[[2]])  )
  
  for(i in 1: (t-1)){
    W.agg.net <- W.agg.net + exp(log.W[i]) * t.nets[,,i]
  }
  
  BMA.net <- BMA.threshold< W.agg.net
  return(BMA.net)
}
