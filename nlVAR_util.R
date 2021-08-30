read.sim.data <- function(dir,folder.name = "sim",dat.fl.name,num.sim, true.net.name = NULL){
  
  dat <- NULL
  for(i in 1:num.sim){
    folder <- paste(folder.name,i,"/",sep="")
    fl.path<- paste(dir,folder,dat.fl.name,sep="")
    dat <- rbind(dat,read.table(file=fl.path,header=T,sep="\t",blank.lines.skip=T))
  }
  
  p <- ncol(dat) - 1
  t <- length(unique(dat$Time))
  n <- nrow(dat)/t
  t.points <- dat$Time[1:t]
  nodes.name <- colnames(dat)[-1]
  
  #n;p;t
  dat <- as.matrix(dat)
  X <- array(NA,dim = c(n,p,t))
  dimnames(X) <- list(NULL,nodes.name,as.character(t.points))
  
  for(time in 1:t){
    X[,,time] <- dat[dat[,"Time"] == t.points[time],-1]
  }
  rm(dat)
  #~~~~~~~~~~~~~~~~~get true network~~~~~~~~~~~~~~~~~~
  true.adj <- NULL
  if(!is.null(true.net.name)){
    true.net <- read.table(paste(dir,true.net.name,sep=""),header = F,stringsAsFactors = F)
    true.adj <- matrix(0,nrow=p,ncol=p, dimnames = list(nodes.name,nodes.name))
    for(i in 1:nrow(true.net)){
      node.from <- true.net[i,1]
      node.to <- true.net[i,2]
      true.adj[node.from, node.to] <- 1
    } 
  }
  
  return(list(X = X, true.adj = true.adj ))
}

# Given a matrix X, this function centers ecah column
# and returned the centered matrix along with the column mean
# as a list
#------------------------------------------------------------
centerX <- function(X){
  col.mean <- apply(X,2,mean)
  X <- sweep(X,2,col.mean)
  return(list(X=X,col.mean=col.mean))
}

# This function block standardize the design matrix X
# by group id and size using QR decomposition
# It then returns the block standardized matrix along
# with the block inverse R matrices needed for transforming
# the coefficients back
#----------------------------------------------------
blk.standardize <- function(X,group.id,spline.deg){
  require(Matrix)
  
  p <- ncol(X)/spline.deg
  if(p!=floor(p))
    stop(paste("col numb is",p,"number of columns of x should be a integer multiple of spline.deg"))
  
  
  # initializging
  #--------------
  Qs <- array(dim = dim(X),dimnames = dimnames(X))
  R.inv <- matrix(0, spline.deg*p , spline.deg*p)
  
  
  for (gp in 1:length(unique(group.id))){
    rows <- cols <- group.id == gp
    decomp <- qr(X[, cols])
    if (decomp$rank < length(sum(cols))) 
      stop("Block belonging to columns ", paste(cols, collapse = ", "), 
           " has not full rank! \n")
    R.inv[rows,cols] <- solve(qr.R(decomp) * 1/sqrt(nrow(X)))
    Qs[, cols] <- qr.Q(decomp) * sqrt(nrow(X))
  }
  
  return(list(Qs = Qs, R.inv = R.inv))
}

#Given X and the smooth.deg, this function returns a list
# containing the colmean of X, block standardized X and the
# block diagonal matrix consiting of inverses of R (of QR decomposition)
#-----------------------------------------------------------------------
get.stand.param <- function(X,smooth.deg){
  grp <- rep(1:(ncol(X)/smooth.deg),each=smooth.deg) 
  stand <- list()
  
  # Center
  #-------
  cent <- centerX(X = X)
  X <- cent$X
  stand$col.mean <- cent$col.mean
  
  # Standardize
  #------------
  blk.st.X <- blk.standardize(X = X,group.id = grp,spline.deg = smooth.deg)
  stand$X <- blk.st.X$Qs
  stand$R.inv <- blk.st.X$R.inv
  
  return(stand)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~fit.STgglasso~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# This funciton receives the design matrix (expanded) and an n-vector Yi perform standardized
# group lasso. It standardizes the design matrix as needed
# Note include.obs.Yi serves to exlude betas (pfs used to do) pfs [,i] == min(pfs[,i]) would be a value
# Or if we want to refit on non-zeor betas beta!=0 or beta[-1] !=0 (intercept) would be good value

fit.STgglasso <- function(X, Yi, intercept,include.obs.Yi,include.cf.Yi = NULL , stand, smooth.deg,lambda){
  require(gglasso)
  
  #+++++++++++++++++++
  #cat("length include.cf.Yi",length(include.cf.Yi),"\n")
  #cat("ncol X",ncol(X),"\n")
  #cat("ncol stand$X", ncol(stand$X),"\n")
  #+++++++++++++++++++
  if(is.null(include.cf.Yi)) 
    include.cf.Yi = rep(T,ncol(X))
  if(length(include.cf.Yi)!=ncol(X))
    stop("Length of include.cf.Yi should equal ncol(X). No intercept element should be included!")
  
  grp <- rep(1:(sum(include.cf.Yi)/smooth.deg),each=smooth.deg)
  
  # Standardize as needed
  #----------------------
  recalculate.stand <- any(!include.obs.Yi) | any(!include.cf.Yi) | is.null(stand) #any(!include.cf.Yi) doesn't need recalc just exclusion and update!!! improve later
  if(recalculate.stand){
    stand <- get.stand.param(X = X[include.obs.Yi, include.cf.Yi],smooth.deg = smooth.deg)
  }
  
  # Center Yi
  #----------
  ybar.i <- mean(Yi[include.obs.Yi])  # This should be close to zero as y is scaled already
  yc.i <- Yi[include.obs.Yi,drop=F] - ybar.i # This also subtracting zero, these two steps are for throughness and effectively unnecessary
  
  
  # Fit gglasso
  #------------
  mod <- gglasso(x=stand$X,  #Note stand$X already has the include.obs and include.cf info factored in
                 y=yc.i,
                 group=grp,
                 loss="ls",
                 lambda=lambda/sqrt(length(yc.i)),
                 #pf = pfs[,i,drop=F], # Note: pf is not longer used to weigh lambda differently instead X columns are excluded
                 intercept =intercept)
  cf <- coef(mod)
  cf[1] <- ybar.i - cf [-1]%*% apply(stand$X[include.obs.Yi,,drop=F],MARGIN = 2,mean) # This should effectively be zero as X is column centered already
  
  # Transfrom back the Solns
  #-------------------------
  cf[-1] <- stand$R.inv %*% cf[-1]
  cf[1] <- cf[1,] - sum(cf[-1] * stand$col.mean)
  
  # Place within beta vector
  #-------------------------
  beta <- rep(0,ncol(X)+ intercept)
  beta[c(intercept,include.cf.Yi)[1:length(beta)+ !intercept]] <- cf
  beta <- as.matrix(beta)
  
  return(beta)
}


#~~~~~~~~~~~~~~~~~~~~~~~~~Computes the exclusion adjacency matrix~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# It receives  
# X.npt: labeled input data array of n x p x t
# trt.target.map: a named (by treatment) list of length equal number of unique treatments each element of which
#                 corresponds to node/s being targeted by the treatment. Note the names need to match those used
#                 in X.npt. Also note there needs to be at least 3 rows of control (as defined by ctrl.name or "control")
#                 in X.npt.
# sig.level: the p.value below which no edge is excluded
# link: A function used to transform data prior to background subtraction
# ctrl.name: If the control is not specifically named as such, it could be identified here

get.exclusion.adj <- function(X.npt, trt.target.map,sig.level = .3,link = log, ctrl.name = NULL){
  require(ROC)
  require(abind)
  
  dm <- dim(X)
  if(length(dm)!=3) stop("An array of nxpxT is expected")
  n <- dm[1]
  p <- dm[2]
  t <- dm[3]
  if(!is.null(ctrl.name)){
    tmp <- names(trt.target.map) 
    tmp[tmp == ctrl.name] <- "control"
    names(trt.target.map) <- tmp
    tmp <- unlist(dimnames(X.npt)[1])
    tmp[tmp == ctrl.name] <- "control"
    dimnames(X.npt) <- list(tmp,unlist(dimnames(X.npt)[2]),unlist(dimnames(X.npt)[3]))
  }
  obs.names <- unname(unlist(dimnames(X.npt)[1]))
  node.names <- unname(unlist(dimnames(X)[2]))
  time.points <- as.numeric(unname((unlist(dimnames(X)[3]))))
  
  if(! all(obs.names %in% names(trt.target.map)))
    stop("row names of the X.npt need to be from the names in trt.target.map")
  
  if(sum(obs.names == "control")<3)
    stop("tr.node.list names need to contrain at least 3 \"control\" ") 
  
  #~~~~~~~~~~~~~~~~~~~compute AUC~~~~~~~~~~~~~~~~~~
  ctrls.means <- array(apply(X.npt[obs.names == "control",,,drop=F],MARGIN = c(2:3),FUN = mean),dim = c(1,p,t))
  lnk.SMN <- sweep(link(X.npt),MARGIN = c(2,3),STATS = link(ctrls.means),"-")
  
  lnk.auc <- apply(aperm(lnk.SMN,perm = c(1,3,2)),MARGIN = c(1,3),
                   FUN = function(y) 
                     trapezint(x = time.points,y,a = time.points[1],b = time.points[t]))
  #~~~~~~~~~~~~~~~~~~~p.value mat~~~~~~~~~~~~~~~~~~~~
  rep.names <- unique(obs.names)
  rep.list <- sapply(rep.names, FUN = function(rep.name) which(obs.names == rep.name),simplify = F,USE.NAMES = T)
  pval.mat <- matrix(NA,nrow = length(rep.names), ncol = length(node.names),
                     dimnames = list(rep.names, node.names) )
  
  ctrl.rows <-  rep.list[["control"]]
  for(rep.name in rep.names){
    
    rep.rows <- rep.list[[rep.name]]
    for(node.name in node.names){
      pval.mat[rep.name,node.name] <- t.test(lnk.auc[rep.rows,node.name],lnk.auc[ctrl.rows,node.name],alternative = "two.sided",var.equal = FALSE)$p.value
    }
  }
  pvalMat.adjusted <- matrix(p.adjust(pval.mat[rep.names != "control",]), nrow = nrow(pval.mat) -1, ncol = ncol(pval.mat))
  dimnames(pvalMat.adjusted) <- list(rep.names[rep.names != "control"],node.names)
  
  #~~~~~~~~~~~~~~~~~~~Excluding ineffective trt~~~~~~~~~~
  # if the targeted node not significant --> renders trt ineffective
  # by setting all corresponding p-values to 0 --> none will be
  # excluded
  
  for(trt in rownames(pvalMat.adjusted)){
    if(!all(trt.target.map[[trt]] %in% node.names))
      next
    if(! all(pvalMat.adjusted[trt,trt.target.map[[trt]]] < sig.level))
      pvalMat.adjusted[trt,] <- 0
  }
  
  #~~~~~~~~~~~~~~~~~~~exclusion mat~~~~~~~~~~~~~~~~~~~~~~
  all.trt <- names(trt.target.map)
  target.trt.map <- lapply(node.names,FUN = function(node)  
    all.trt[unlist(lapply(trt.target.map,function(trt)any(node %in% trt) ))])
  names(target.trt.map) <- node.names
  
  exclusion.adj <- matrix(FALSE,p,p,dimnames = list(node.names, node.names))
  
  sigMat <- pvalMat.adjusted < sig.level
  for(node.name in node.names){
    exclusion.adj[node.name,] <-!apply (sigMat[target.trt.map[[node.name]], ,drop=F],MARGIN = 2,
                                        FUN = function(col){
                                          if(length(col)==0)
                                            TRUE
                                          else
                                            any(col)
                                        }) 
  }
  
  return(exclusion.adj)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~Gets the penalty factors for gglasso~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# This function receives a logical matrix in the dimension of
# the pxp adjacency matrix and generates the appropriate pentaly
# factor matrix for the gglasso
# test case:
#           grp <- c(c(1,1,1,2,2,5,4,4,4,4,4),10+c(1,1,1,2,2,5,4,4,4,4,4),20+c(1,1,1,2,2,5,4,4,4,4,4))
#           exclude.adj <- matrix(F,4,4)
#           exclude.adj[2,3] <- TRUE

get.penFactors <- function(exclude.adj,grp,large.weight = 1000){
  if(!is.matrix(exclude.adj) | !is.logical(exclude.adj[1])|ncol(exclude.adj)!=nrow(exclude.adj)) 
    stop("exclude.adj needs to be a pxp logical matrix")
  p <- nrow(exclude.adj)
  t <- length(unique(grp))/p
  grp.count <- table(grp)[as.factor(grp[!duplicated(grp)])]
  pf.yi <- sqrt(grp.count)
  pfs <- matrix(pf.yi, nrow = p * t, ncol = p,byrow = FALSE)
 
  #stack t exclude.adj on top of one another
  stk.ex.adj <- do.call("rbind", rep(list(exclude.adj),t))
  pfs[stk.ex.adj] <- large.weight

  return(pfs)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Expands the X using Bspline or polynomial~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Xpand <- function(X.npt,smooth.deg, smoother){
  require(splines)
  n <- nrow(X.npt)
  p <- ncol(X.npt)
  k <- smooth.deg
  
  d.names <- dimnames(X.npt)
  if(is.null(d.names))
    stop("Each of n observations, p nodes and t timepoints in the data matrix X need to be named. Repeated names are ok!")
  
  #~~~~~~~~~~~~~~expansion function~~~~~~~~~~~~
  if(smoother == "nspline"){
    xfun <- function(x,degree){
      return(ns(x = x,df = degree))
    }
  }
  if(smoother == "bspline"){
    xfun <- function(x,degree){
      return(bs(x = x,degree = degree))
    }
  }
  if(smoother == "polynomial"){
    xfun <- function(x, degree){
      xx <- NULL
      for(i in 1:degree)
        xx <- cbind(xx,x^i)
      return(xx)
    } 
  }
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  Xpanded <- matrix(NA, ncol = k * p, nrow = n)
  
  for (j in 1:p) {
    #Xpanded[,p * (1:k - 1) + j ] <-
    #  svd(apply(bs(X.npt[, j]), 2, scale, scale = F))$u * sqrt(n - 1)
    Xpanded[,1:k + (j-1)*k ] <- xfun(x = X.npt[, j],degree= smooth.deg)
  }
  colnames(Xpanded) <- paste(rep(d.names[[2]],each = k),".",1:smooth.deg,sep="")
  rownames(Xpanded) <- d.names[[1]]
  return(Xpanded)
}

#~~~~~~~~~~~~~~~~~~~~~~~Gets the random kfold indices~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Given the length of indices, n, and the number of folds, kfold,
# it generates kfolds disjoint sets of indices the union of which
# is 1:n
get.kfold.index <- function(n,kfold){
  if(n<kfold) stop(paste("n needs to be at least as large as", kfold))
  split.factor <- c(ceiling(1:(n-n%%kfold)/ ((n-n%%kfold)/kfold)), sample(1:kfold, size = n%%kfold,replace = TRUE))
  splt <- split(x=sample(1:n),f=split.factor)
  return(splt)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Computes the Cost (BIC or RSS)~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Given an nxq matrix X, an nxp matrix Y, an qxp beta matrix, and an nxp boolean include.obs.mat
# it will compute a p-vector of cost values either BIC or RSS
# X: nxq design matrix
# Y: nxp response matrix measured
# betas: A qxp beta matrix
# intercept: If TRUE, the 1st row of beta matrix needs to be the intercept
# include.obs.mat: nxp boolean matrix, ith column of which determins which row of observations
# include.cf.mat: pxp boolean matrix, ith column of which determins which cols of X (and corresponding betas) were included in fitting
# should be included in computing rss for the ith column of Y.
# cost.type: RSS or BIC of the fit
# BIC.parm: It needs to be a list containing regression pentaly, lambda, the E-BIC constant, gamma.const, and a q-vector, grp,
#           specifying the grouping of betas for the nxq design matrix X
# SGL: If true it uses standardized group lasso dof

get.cost <- function(Y, X, betas, intercept,include.obs.mat =NULL, include.cf.mat = NULL,cost.type = c("RSS", "BIC"), BIC.parm = NULL, SGL = !is.null(BIC.parm$stand)){

  #+++++
  #cat(BIC.parm$gamma.const,"\n")
  #+++++
  
  # Checking the input
  #-------------------
  cost.type <- match.arg(cost.type)
  if(cost.type == "BIC" & is.null(BIC.parm))
    stop("BIC.parm needs to be a list containing regression pentaly, lambda, the E-BIC constant, gamma.const, and a q-vector, grp,
         specifying the grouping of betas for the nxq design matrix X")
  
  # Cast formatting and initializing
  #---------------------------------
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  betas <- as.matrix(betas)
  
  if(nrow(betas)!=ncol(X)+intercept)
    stop("length of betas should be ncol(X) + intercept (where intercept is 0/1)")
  
  if(is.null(include.obs.mat)) include.obs.mat <- Y==Y
  include.obs.mat <- as.matrix(include.obs.mat)
  
  # Check include.cf.mat format and cast into correct length
  if(is.null(include.cf.mat)){
    include.cf.mat <- (betas = betas)[intercept + 1:ncol(X)]
  }
  if(nrow(include.cf.mat) != ncol(X))
    stop("include.cf.mat should equal ncol(X). It shouldn't include element for intercept!")
  
  include.cf.mat <- rbind(intercept,include.cf.mat)
  X <- cbind(1,X)
  if(!intercept){
    betas <- rbind(0,betas)
  }
  
  RSS.per.y <- DOF.per.y <- rep(NA,ncol(Y))

  # Compute RSS and BIC
  #--------------------
  for(i in 1:ncol(Y)){
    include.obs <- include.obs.mat[,i]
    include.cf <- include.cf.mat[,i]
    
    if(intercept < sum(include.cf)){
      #+++++++++++++++++++
      #cat("***********************\n")
      #cat("sum(include.cf)",sum(include.cf),"\n")
      #cat("dim(X[include.obs,include.cf,drop=F])", nrow(X[include.obs,include.cf,drop=F]),",")
      #cat(ncol(X[include.obs,include.cf]),"\n")
      #cat("dim(betas[include.cf,i,drop=F])", nrow(betas[include.cf,i,drop=F]),",")
      #cat(ncol(betas[include.cf,i,drop=F]),"\n")
      #+++++++++++++++++++++++++++++++++++++++++++
      RSS.per.y[i] <- sum(as.matrix((Y[include.obs,i,drop=F] - X[include.obs,include.cf,drop=F] %*% betas[include.cf,i,drop=F])^2))
      
      #++++++
      #cat("sum(include.cf)",sum(include.cf),"\n")
      #cat("indices", which(include.cf)+ intercept,"\n")
      #cat("Length beta.hat",length(betas[ which(include.cf)+ intercept ,i]),"\n")
      #cat("length grp",length(grp),"\n")
      #cat("length group.id",length(BIC.parm$grp[include.cf]),"\n")
      #++++++
      
      # BIC
      #----
      if(cost.type == "BIC"){
        DOF.per.y[i] <- intercept + get.GL.DOF(X = X[include.obs, c(F,include.cf[-1]),drop=F]
                                               , beta.hat = betas[c(F,include.cf[-1]) ,i,drop=F]
                                               , lambda = BIC.parm$lambda
                                               , group.id = BIC.parm$grp[include.cf[-1],drop=F]
                                               , SGL = SGL)
      }
    }
    else{
      RSS.per.y[i] <- sum((Y[include.obs,i] - intercept * mean(Y[include.obs,i]))^2)
        #sum((sweep(x = Y[include.obs,i], MARGIN = 2, intercept * apply(Y[include.obs,i],2,mean)))^2)
      DOF.per.y[i] <- 1 * intercept
    }
  }
  
  # Set the cost to RSS
  #--------------------
  cost <- RSS.per.y
  
  # Set the cost to BIC
  #--------------------
  if(cost.type == "BIC"){
    reg.obs.num <- apply(include.obs.mat,MARGIN = 2,FUN = sum)
    #NOTE TWO DIFFERENT FORMULATION OF BIC:
    # 1- ZOU, HASTIE, TIBSHIRANI: ON THE DEGRESS OF FREEDOM OF THE LASSO, note for this sigma2.hat is needed!
    # 2- SCHWARZ 1978: ESTIMATING THE DIMENSION OF A MODEL, specialized for iid normal noise
    # 3- CHEN 2008: EXTENDED BAYESIAN INFORMAITON CRITERIA FOR MODEL SELECITON WITH LARGE MODEL SPACES.
    
    #1)
    #BICs <- rss.s/(reg.obs.num * sigma2.hat) + log(reg.obs.num)* DOF/reg.obs.num
    
    #2)
    BIC.per.y  <- reg.obs.num * log(RSS.per.y/reg.obs.num) + log(reg.obs.num)* DOF.per.y
    
    #3)
    if(all(intercept < apply(include.cf.mat[-1,],2,sum))){
      # note this gamma.const> 1-1/(2k) where k is defined in p=O(n^k) e.g. p=O(n) then gamma.const > .5
      BIC.per.y <- BIC.per.y + 2* DOF.per.y * BIC.parm$gamma.const * log(apply(include.cf.mat[-1,],2,sum)) # NOTE: if gamma.const is zero BIC is Schwartz
      
    }
    cost <- BIC.per.y
  }
  #++++++++++++++++++++
  #if(any(BIC.per.y == -Inf)){
  #  cat("RSS: ")
  #  print(RSS.per.y)
  #  cat("log(apply(include.cf.mat[-1,],2,sum)): ")
  #  print(log(apply(include.cf.mat[-1,],2,sum)))
  #}
  #print("cost is")
  #print(cost)
  
  #++++++++++++++++++++
  
  return(cost)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~Thresholds the Adj matrix~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
threshold<- function(A.pptd, tau, beta.err, threshold.method){
  dm <- dim(A.pptd)
  p <- dm[2]
  t <- dm[3] #Note this is T-1
  smooth.deg <- dm[4]
  
  A.norm <- array(0,dim=c(dm[1:3],1))
  A.pptd.2 <- A.pptd^2
  
  for(i in 1:smooth.deg){
    A.norm[,,,1] <- A.norm[,,,1,drop=FALSE] + A.pptd.2[,,,i,drop=FALSE]
  }
  
  A.norm <- sqrt(A.norm) # Note: A.norm is ||A|| and not ||A||^2
  norm0.Th <- rep(p^2 * beta.err/(t),t)  #This is the version as given in the paper
  if(threshold.method == "decaying"){
    norm0.Th <- p^2 * beta.err /(t+1 - t:1) #This is similar to the version implemented in Thlasso code
    #cat(threshold.method,"\n")
  }
    
  Th.A.norm <- array(0,dim = dim(A.norm))  
  VAR.order <- 0
  
  ThbyLMNT.A.norm <- A.norm *(abs(A.norm)> tau) # only for order calc
  for(i in t:1){
    norm0 <- sum(A.norm[,,i,,drop=FALSE] !=0)
    ThbyLMNT.norm0 <- sum(ThbyLMNT.A.norm[,,i,,drop=FALSE] !=0) # only for order calc
    if(norm0 < norm0.Th[i]){
      next # next effectively sets the whole A(t) --> 0
    }else{
      #VAR.order <- t+1-i
      if(norm0.Th[i] <= ThbyLMNT.norm0) # only for order calc
        VAR.order <- t+1-i
      #Th.A.norm[,,i,] <- A.norm[,,i,,drop=FALSE]*(abs(A.norm[,,i,,drop=FALSE]) > tau)
      Th.A.norm[,,i,] <- ThbyLMNT.A.norm [,,i,,drop=FALSE]
    }
  }
  
  Th.A.filter <- Th.A.norm !=0
  Th.A <- array(0,dim=dm)
  for(i in 1:smooth.deg){
    Th.A[,,,i] <- A.pptd[,,,i,drop=FALSE]* Th.A.filter
  }
  #~~~~~~~~Format Th.A.norm~~~~~~~~~~~
  Th.A.norm <- array(Th.A.norm,dim= dm[1:3])
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  return(list(Th.A = Th.A, Th.A.norm = Th.A.norm, VAR.order = VAR.order ))
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Converting the betas to Adj matrix~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Receives a  (smooth.deg*p*t) x p matrix of betas and reshape it into an
# array of p x p x t x smooth.deg adjacency matrix
# Note: the value of p w/ and w/o intercept differs by 1
#       however the returned adj never include intercept!

coef.to.Adj <- function(betas,intercept, smooth.deg){
  dm <- dim(betas)
  if(length(dm)!= 2) 
    stop("Incompatible betas dimension, betas need to be (p*smooth.deg*t) x p or (p*smooth.deg*t + 1) x p if intercept is TRUE")
  p <- dm[2]
  Px <- p*smooth.deg 
  t <- (nrow(betas) - intercept)/Px 
  if(!any(dm == c(Px * t + intercept,p)))
    stop("Incompatible betas dimension, betas need to be (p*smooth.deg*t) x p or (p*smooth.deg*t + 1) x p if intercept is TRUE")
  if(intercept){
    b0 <- betas[1,]
    betas <- betas[-1,]
  }
  Adj <- aperm(array(betas,dim = c(smooth.deg,p,t,p)),perm = c(2,4,3,1))
  return(Adj)
}
#~~~~~~~~~~~~Converts adjacancy matrix back to ceof matrix ~~~~~~~~~~~~~
Adj.to.coef <- function(adj,beta0 = NULL,smooth.deg){
  dm<- dim(adj)
  p<- dm[2]
  t<- dm[3]
  #adj
  coef <-rbind(beta0,
               array(aperm(adj,perm=c(4,1,3,2)),dim=c(smooth.deg*p*t,p))) 
  return(coef)
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~scale x~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
scale.XArray <- function(X,n,t){
  for (i in 1:t){
    X[,,i] = scale( X[,,i] )
  }
  X = sqrt(n/(n-1))* X
  return(X)
}
#~~~~~~~~~~~~~~~~~~~~get.XY~~~~~~~~~~~~~~~~~~
get.XY <- function(X,n,p,t,smooth.deg,smoother){
  d.names <- dimnames(X)
  if(is.null(d.names))
    stop("Each of n observations, p nodes and t timepoints in the data matrix X needs to be named. Repeated names are ok!")
  rnames <- d.names[[1]]
  cnames <- rep(d.names[[2]],t)
  matX <- matrix(X,nrow=n,ncol=p*t,dimnames = list(rnames,cnames))
  Y <- matX[,(1:p) + (t-1)*p,drop= FALSE]
  XX <- matX[, -((1:p) + (t-1)*p), drop = FALSE]
  XXpanded <- Xpand(X.npt=XX,smooth.deg=smooth.deg, smoother = smoother) 
  return(list(XXpanded=XXpanded, Y=Y))
}
# Given the data matrix Y.np and and the mapping for treatment to target node:
# trt.taget.map, it returns a nxp matrix each columns of which corresponds to a node
# such that in the column corresponding to node A the elements of the column are TRUE 
# if in the corresponding observation node is not targeted by the treatment.

get.include.obs.mat <- function(Y.np, trt.target.map){
  
  obs.names <- rownames(Y.np)
  node.names <- colnames(Y.np)
  
  #~~getting the target.trt.mapt~~~
  all.trt <- names(trt.target.map)
  target.trt.map <- lapply(node.names,FUN = function(node)  
    all.trt[unlist(lapply(trt.target.map,function(trt)any(node %in% trt) ))])
  names(target.trt.map) <- node.names
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  include.obs.mat <- Y.np == Y.np
  
  #exclude.obs.list <- list()
  #length(exclude.obs.list) <- length(node.names)
  #names(exclude.obs.list) <- node.names
  
  for(node in node.names){
    trts <- target.trt.map[[node]]
    #exclude.obs.list[[node]] <- which(!is.na(pmatch(obs.names, trts,duplicates.ok = TRUE)))
    exclude.obs<- which(!is.na(pmatch(obs.names, trts,duplicates.ok = TRUE)))
    include.obs.mat[exclude.obs, node ] <- FALSE
  }
  return(include.obs.mat)
}

#~~~~~~~~~~~~~~~~~~~~~~~~CVL.var~~~~~~~~~~~~~~~~~~~~~~~~~
# This function receives an nxp matrix X and an n-vector
# y and estimates the variance of y using the CV method
# This funciton is called by get.variance
# Note see Reid et al. 2014 "A study of Error Variance Estimation in Lasso Reression"

CVL.var <- function(X,y){
  n <- nrow(X)
  mod <- cv.glmnet(x = X,y = y, nfolds = 5)
  cf <- coef(mod)
  sigma2.hat <- sum((y - cbind(1,X) %*% cf )^2)/(n - sum(cf!=0))
  return(sigma2.hat)
}


#~~~~~~~~~~~~~~~~~~~~~~~~CVL.var~~~~~~~~~~~~~~~~~~~~~~~~~
# This function receives an nxp matrix X and an n-vector
# y and estimates the variance of y using the RCV method
# This funciton is called by get.variance
# Note see Reid et al. 2014 "A study of Error Variance Estimation in Lasso Reression"

RCV.var <- function(X,y){
  n <- nrow(X)
  X.divisions <- get.kfold.index(n=n,kfold=2)
  
  X1 <- X[X.divisions[[1]], ]
  y1 <- y[X.divisions[[1]]]
  n1 <- length(X.divisions[[1]])
  
  X2 <- X[X.divisions[[2]], ]
  y2 <- y[X.divisions[[2]]]
  n2 <- length(X.divisions[[2]])
  
  cf.1 <- coef(cv.glmnet(X1, y1, nfolds = 5))
  index.1 <- which(cf.1 !=0)
  
  cf.2 <- coef(cv.glmnet(X2, y2, nfolds = 5))
  index.2 <- which(cf.2 !=0)
  
  #~~~~~~~~~~refit across~~~~~~~~~~~~
  #var index1 data2
  dat.2 <- cbind(1,X2)[, index.1, drop=F]
  refit.1 <- lm(y2 ~ 0 + dat.2)
  
  #var index2 data1
  dat.1 <- cbind(1,X1)[, index.2, drop=F]
  refit.2 <- lm(y1 ~ 0 + dat.1)
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  res.1 <- refit.1$residuals
  res.2 <- refit.2$residuals
  
  v1 <- sum(res.1^2/(n2 - sum(cf.1 != 0)))
  v2 <- sum(res.2^2/(n1 - sum(cf.2 != 0)))
  
  sigma2.hat <- (mean(v1,v2))
  return(sigma2.hat)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~get.variance~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Estimates variance of the last time point in the VAR
# X.npt: the orginal data array
# include.obs.mat: the observations to inclue for each of the p regressions
# Returns for each node (column of Y) an estimated variance (\hat{sigma^2})
# Note see Reid et al. 2014 "A study of Error Variance Estimation in Lasso Reression"

get.variance <- function(X.npt,include.obs.mat, scale = TRUE,est.method = c("RCV","CVL")){
  require(glmnet)
  dm <- dim(X.npt)
  if(length(dm)!=3) stop("An array of nxpxT is expected")
  n <- dm[1]
  p <- dm[2]
  t <- dm[3]
  #~~~~~~~~~~~~~~Scale X ~~~~~~~~~~~~
  if(scale)
    X.npt <- scale.XArray(X.npt,n,t)
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  matX <- matrix(X.npt,nrow=n,ncol=p*t)
  Y <- matX[,(1:p) + (t-1)*p,drop= FALSE]
  XX <- matX[, -((1:p) + (t-1)*p), drop = FALSE]
  sigma2.hat <- rep(NA,p)
  
  #~~~~set estimation function~~~
  est.method <- match.arg(est.method)
  est.var <- RCV.var
  if(est.method == "CVL")
    est.var <- CVL.var
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  for(i in 1:p){
    include.obs <- include.obs.mat[,i]
    #mod <- cv.glmnet(x = XX[include.obs,],y = Y[include.obs,i],nfolds = 5)
    #cf <- coef(mod)
    #sigma2.hat[i] <- sum((Y[include.obs,i] - cbind(1,XX[include.obs,]) %*% cf )^2)/(sum(include.obs)- sum(cf!=0))
    sigma2.hat[i] <- est.var(X = XX[include.obs,], y = Y[include.obs,i])
  }
  return(sigma2.hat)
}


#~~~~~~~~~~~~~~~~get.GL.DOF~~~~~~~~~~~~~~~~~~~~~~~
# Given the design matrix X and a vector beta.hat estimated
# via group lasso and the corresponding group.id (which
# needs to be the same length as betas in which the jth element
# is the group id of the jth beta) and output the DOF
# note the we need to have X[,active.set].... term be invertible!
# IF SGL is TRUE, it uses the estimator for standardized group lasso
# Otherwise it uses the estimator for the group lasso

#Note this expects all groups to be the same size and in sequence (not interdispersed)

get.GL.DOF <- function(X, beta.hat, lambda, group.id, SGL = F){
  require(Matrix)
  if(length(beta.hat) != length(group.id))
    stop("Length of beta.hat needs to be equal to length of group.id")
  tb.group <- table(group.id)
  if(!all(tb.group ==tb.group[1]))
    stop("all groups need to be the same size")
  if(any(c(group.id,group.id[length(group.id)])-c(1,group.id) < 0))
    stop("all groups need to be in blocks (not inter-dispersed) and in increasing order")
  
  #Cast group.id into sequential order
  #-----------------------------------
  group.id <- as.numeric(factor(group.id,labels=1:length(unique(group.id))))
  
  active.set <- which(beta.hat != 0)
  if(length(active.set) == 0) # DOF is zero if all betas are zero
    return(0)
  groups.active <- group.id[active.set]
  bHat.active <- beta.hat[active.set]
  sz.active <- length(bHat.active)
  sz.group <- tb.group[1]
  num.group <- sz.active/sz.group
  
  
  # If the solver was regular group lasso
  #--------------------------------------
  if(!SGL){
    block.diagonizer <- bdiag( rep(list(matrix(1,sz.group,sz.group)),num.group))
    
    # Normalizing operator
    norm.blocks <- tapply(beta.hat^2,INDEX = group.id,FUN = sum,simplify = TRUE)[groups.active]
    normalizer.mat <- matrix(1/sqrt(norm.blocks),nrow = sz.group*num.group,ncol=sz.group * num.group,byrow = TRUE)
    
    # Projection block diagonal operator
    P_bet <- bHat.active%*%t(bHat.active) * block.diagonizer * normalizer.mat^2
    Id <- diag(1,nrow = sz.group*num.group, ncol = sz.group *num.group)
    D <- normalizer.mat * (Id - P_bet)
  }
  
  # If the solver was standardized group lasso
  #-------------------------------------------
  else{
    D <- matrix(0, nrow = sz.active, ncol = sz.active)
    for(gr in 1:num.group){
      blk.ndx <- groups.active %in% unique(groups.active)[gr]
      group.b <- which(group.id %in% groups.active[blk.ndx])
      X.gr <- X[,group.b]
      beta.gr <- beta.hat[group.b]
      y.hat.gr <- X.gr %*% beta.gr
      
      nrm2.y.hat.gr <- as.numeric(t(y.hat.gr) %*% y.hat.gr)
      D[blk.ndx, blk.ndx] <- t(X.gr) %*% X.gr / sqrt(nrm2.y.hat.gr)  - 
        t(X.gr) %*% y.hat.gr %*% t(y.hat.gr) %*% X.gr / (nrm2.y.hat.gr)^1.5
      
    }
  }
  
  
  # cov(y.hat, y)
  #--------------
  inv.term <- tryCatch(expr = solve(t(X[,active.set,drop=F]) %*% X[,active.set,drop=F] + (lambda*sqrt(sz.group)) * D )
                       ,error = function(exception){
                         return(exception)})
  
  if(inherits(inv.term,"simpleError"))
    stop(paste(inv.term$message,"\n  ***Singularity usually occurs when regression was done with too smal of a lambda.If system is singular try larger lambda***"))
  
  cov.yHat.y <- X[,active.set,drop=F] %*%  inv.term %*% t(X[,active.set,drop=F])
  GL.DOF <-sum(diag(cov.yHat.y))
  return(GL.DOF)
}




#~~~~~~~~~~~~~~~~~~~~~~~VAR~~~~~~~~~~~~~~~~~~~~~~~~~~
# The core of nlVAR after X is scaled, expanded and separated from Y
# It regresses each of the p columns of Y on the X, excluding observations
# in which the corresponding node wasn't targeted and was free to vary.
# it computes p columns of length npt of the estimates in a npt x p matrix format
# it returns the beta.mat as (smooth.deg*p*t) x p matrix of betas
# Note: if get.adj = TRUE, it returns converts  the beta.mat into 
# array of p x p x t x smooth.deg adjacency matrix
# Note: the value of p w/ and w/o intercept differs by 1

VAR <- function(X, Y,lambda, pfs, smooth.deg=3, intercept = TRUE, get.adj = F, include.obs.mat, stand = NULL){
  require(gglasso)
  p <- ncol(Y)
  n <- nrow(Y)
  obs.names <- rownames(Y)
  node.names <- colnames(Y)
  
  #grp <- rep(1:(ncol(X)/smooth.deg),each=smooth.deg)
  
  num.betas <- 1 * intercept + ncol(X)
  beta.mat <- matrix(0,nrow = num.betas, ncol = p) # Row num = edge origin ------------> Col num = edge destination
  
  #Change pfs to include.cf
  #------------------------
  include.cf.mat <- matrix(rep(pfs == sqrt(smooth.deg),each = smooth.deg),nrow = smooth.deg*nrow(pfs),ncol=ncol(pfs))
  
  # If not standardized
  #----------------------
  if(is.null(stand)){
    for(i in 1:p){
      include.obs <- include.obs.mat[,i]
      include.cf <- include.cf.mat[,i]
      
      # Ensure there is at least one group
      if(smooth.deg <= sum(include.cf)){
        grp <- rep(1:(sum(include.cf)/smooth.deg),each=smooth.deg)

        mod <- gglasso(x=X[include.obs,include.cf,drop=F],
                       y=Y[include.obs,i,drop=F],
                       group=grp,
                       loss="ls",
                       lambda=lambda/sum(include.obs), # scales by the effective number of obs
                       #pf = pfs[,i,drop=F],
                       intercept =intercept)
        
        cf <- coef(mod)
        #+++++++++++
        #cat("dim(beta.mat)",dim(beta.mat),"\n")
        #cat("sum(include.cf)",sum(include.cf),"\n")
        #cat("length(cf[1 + which(include.cf)])", length(cf[1 + which(include.cf)]),"\n")
        #cat("beta.mat[intercept + which(include.cf),i])",length(beta.mat[intercept + which(include.cf),i]),"\n")
        #+++++++++++
        beta.mat[intercept + which(include.cf),i] <- cf[-1]
        if(intercept)
          beta.mat[1,i] <- cf[1]
      }
      else{
        if(intercept)
          beta.mat[1,i] <- mean(Y[include.obs,i,drop=F])
      }
    }# Over p columns of Y 
  }
  # If standardized
  #----------------
  else{
    for(i in 1:p){
      # Ensure there is at least one group
      if(smooth.deg <= sum(include.cf.mat[,i])){
        beta.mat[,i] <- fit.STgglasso(X = X
                                      ,Yi = Y[,i]
                                      ,intercept = intercept
                                      ,include.obs.Yi = include.obs.mat[,i]
                                      ,include.cf.Yi = include.cf.mat[,i]
                                      ,stand = stand
                                      ,smooth.deg = smooth.deg
                                      ,lambda = lambda)
      }
      else{
        if(intercept)
          beta.mat[1,i] <- mean(Y[include.obs.mat[,i],i,drop=F])
      }

    }
  }
  
  input <- list(XX= X, YY = Y, include.obs.mat = include.obs.mat, include.cf.mat = include.cf.mat, stand = stand)
  if(!get.adj)
    return(list(beta.mat=beta.mat, input = input))
  
  adj.mat <- coef.to.Adj(betas = beta.mat,intercept=intercept,smooth.deg=smooth.deg)
  return(list(beta.mat=beta.mat, adj.mat = adj.mat, input = input))
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~tune.CV~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# This function tunes by Cross Validation, the thresholding parameter, cRatio, and the regression penalty, lambda, simultaneously
# XXpanded: Expanded design matrix
# Y: nxp matrix of response variables
# trt.target.map: a named (by treatment) list of length equal number of unique treatments each element of which
#                 corresponds to node/s being targeted by the treatment. Note the names need to match those used
#                 in X.npt. Also note there needs to be at least 3 rows of control (as defined by ctrl.name or "control")
#                 in X.npt. 
# pfs: a (p*t) X p matrix of regression penalty factors. This is used to weight penalties differntly for differnt groups
#      of coefficients, but is mainly used to force certain groups of betas to zero (by weighting them by a large number)
# lmbs: a vector of group lasso regression penalties, lambdas. One of the two tunning parameters
# cRatios: a vector of constant by which lambda will be scaled to form the thresholding parameter tau. One of the two tunning
#          parameters
# kfold: Number of cross validation folds
# smooth.deg: The degree of smoother used for expansion of the design matrix
# intercept: Whether to fit model with non-zero intercept
# threshold.method: Either "uniform" or "decaying". Used to determine the thresholding scheme. See threshold

tune.CV <- function(XXpanded, Y,trt.target.map, pfs,lmbs,cRatios, kfold, smooth.deg, intercept, threshold.method){
  n <- nrow(Y)
  lmbs.len <- length(lmbs)
  cRatio.len <- length(cRatios)
  X.divisions <- get.kfold.index(n=n,kfold=kfold)
  
  rss.lcf <- array(NA,c(lmbs.len,cRatio.len,kfold)) # l is lambda, c is cRatio and f s fold
  
  for(i in 1:lmbs.len){
    for( fold in 1:kfold){
      
      X1 <- XXpanded[(1:n)[-X.divisions[[fold]]],,drop=F]
      Y1 <- Y[(1:n)[-X.divisions[[fold]]],,drop=F]
      X2 <- XXpanded[X.divisions[[fold]],,drop=F]
      Y2 <- Y[X.divisions[[fold]],,drop=F]
      #~~~~~~~~~~~make new include map~~~~~~~~~~~~~~
      fold.include.obs.mat1 <-get.include.obs.mat(Y.np = Y1, trt.target.map)
      fold.include.obs.mat2 <-get.include.obs.mat(Y.np = Y2, trt.target.map)
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      nlVAR.fit <- VAR(X=X1,
                       Y=Y1,
                       pfs = pfs,
                       lambda=lmbs[i],
                       smooth.deg = smooth.deg, 
                       intercept = intercept,
                       get.adj=T,
                       include.obs.mat = fold.include.obs.mat1)
      
      b0 <- NULL
      if(intercept)
        b0<- nlVAR.fit$beta.mat[1,]
      
      for(k in 1:cRatio.len){
        tau <- cRatios[k] * lmbs[i]
        Th.adj <- threshold(A.pptd=nlVAR.fit$adj.mat, tau=tau, beta.err=beta.err,threshold.method = threshold.method )$Th.A
        Th.coef <- Adj.to.coef(adj=Th.adj,beta0=b0,smooth.deg=smooth.deg)
        rss.lcf[i,k,fold] <- sum(get.cost(Y=Y2
                                          ,X=X2
                                          ,betas=Th.coef
                                          , intercept = intercept
                                          ,include.obs.mat = fold.include.obs.mat2 
                                          ,cost.type = "RSS"
                                          ,BIC.parm = NULL))
        
      }#cRatios
      
    }#fold
    
  }#lambda
  
  rss.lc <-  apply(rss.lcf,MARGIN = c(1,2),FUN = mean)
  rss.lc.se <- apply(rss.lcf,MARGIN = c(1,2),FUN = sd)
  
  rss.min <- min(rss.lc)
  rss.min.index <- which(rss.lc == rss.min,arr.ind = FALSE)[1] # Gets the first if there are multiple
  col.index <- ceiling(rss.min.index/lmbs.len)
  row.index <- rss.min.index - (col.index - 1)* lmbs.len
  
  lambda.min <- lmbs[row.index]
  cRatio.min <- cRatios[col.index]
  
  #bound.1se <- get.withinSD.indices(mat = rss.lc,sd.mat = rss.lc.se)
  #rss.min.1se.element <- max(which(bound.1se))
  #rss.min.1se.index <- c(rss.min.1se.element %% lmbs.len + (rss.min.1se.element %% lmbs.len ==0)*lmbs.len, rss.min.1se.element / lmbs.len + 1)
  #rss.min.1se <- rss.lc[rss.min.1se.index]
  #lambda.min.1se <- lmbs[rss.min.1se.index[1]]
  #cRatio.min.1se <- cRatios[rss.min.1se.index[2]] 
  
  tune.profile <- list(pen = rss.lc
                       ,pen.se = rss.lc.se
                       ,lambda.min = lambda.min
                       #,lambda.min.1se = lambda.min.1se
                       ,cRatio.min = cRatio.min
                       #,cRatio.min.1se = cRatio.min.1se
                       ,tune.method = "CV")
  
  return(tune.profile)
  
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~tune.BIC~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# This function tunes by BIC, the thresholding parameter, cRatio, and the regression penalty, lambda, simultaneously
# XXpanded: Expanded design matrix
# Y: nxp matrix of response variables
# trt.target.map: a named (by treatment) list of length equal number of unique treatments each element of which
#                 corresponds to node/s being targeted by the treatment. Note the names need to match those used
#                 in X.npt. Also note there needs to be at least 3 rows of control (as defined by ctrl.name or "control")
#                 in X.npt. 
# pfs: a (p*t) X p matrix of regression penalty factors. This is used to weight penalties differntly for differnt groups
#      of coefficients, but is mainly used to force certain groups of betas to zero (by weighting them by a large number)
# lmbs: a vector of group lasso regression penalties, lambdas. One of the two tunning parameters
# cRatios: a vector of constant by which lambda will be scaled to form the thresholding parameter tau. One of the two tunning
#          parameters
# smooth.deg: The degree of smoother used for expansion of the design matrix
# intercept: Whether to fit model with non-zero intercept
# threshold.method: Either "uniform" or "decaying". Used to determine the thresholding scheme. See threshold
# gamma.const: a value between 0 and 1 used to compute EBIC. When zero it will Schwartz BIC 

tune.BIC <- function(XXpanded, Y,trt.target.map, pfs,lmbs,cRatios, smooth.deg, intercept, threshold.method, gamma.const,stand,beta.err = .1){
  n <- nrow(Y)
  lmbs.len <- length(lmbs)
  cRatio.len <- length(cRatios)
  beta.errs <- as.numeric(cRatios!=0) * beta.err
  include.obs.mat <-get.include.obs.mat(Y.np = Y, trt.target.map = trt.target.map)
  grp <- rep(1:((ncol(XXpanded))/smooth.deg),each=smooth.deg)
  
  BIC.lc <- array(NA,c(lmbs.len,cRatio.len)) # l is lambda, c is cRatio
  BIC.parm <- list(lambda = NA, gamma.const = gamma.const, grp = rep(1:(ncol(XXpanded)/smooth.deg),each = smooth.deg), stand = stand)
  
  for(i in 1:lmbs.len){
    
    nlVAR.fit <- VAR(X=XXpanded,
                     Y=Y,
                     pfs = pfs,
                     lambda=lmbs[i],
                     smooth.deg = smooth.deg,
                     intercept = intercept,
                     get.adj= TRUE,  
                     include.obs.mat = include.obs.mat,
                     stand = stand)
    
    
    b0 <- NULL
    if(intercept)
      b0<- nlVAR.fit$beta.mat[1,]
    
    for(k in 1:cRatio.len){
      tau <- cRatios[k] * lmbs[i]
      Th.adj <- threshold(A.pptd=nlVAR.fit$adj.mat, tau=tau, beta.err=beta.errs[k],threshold.method = threshold.method )$Th.A
      Th.coef <- Adj.to.coef(adj=Th.adj,beta0=b0,smooth.deg=smooth.deg)
      
      # Refitting
      #-----------
      
      include.cf.mat <- Th.coef[1:ncol(XXpanded)+ intercept,] !=0
       
      #collapase by group to make pfs
      pfs.th <- sqrt(smooth.deg)*(as.matrix(aggregate(.~grp,data = data.frame(grp=grp,include.cf.mat),FUN = sum)[,-1])!=0)
      pfs.th [pfs.th==0] <- 10000*sqrt(smooth.deg)
      
      nlVAR.refit <- VAR(X=XXpanded,
                       Y=Y,
                       pfs = pfs.th,
                       lambda=lmbs[i],
                       smooth.deg = smooth.deg,
                       intercept = intercept,
                       get.adj= TRUE,  
                       include.obs.mat = include.obs.mat,
                       stand = stand)
      
      
      # Compute BIC for the fit
      #------------------------
      BIC.parm$lambda = lmbs[i]
      BIC.lc[i,k] <- mean(get.cost(Y=Y
                                   ,X = XXpanded
                                   ,betas = nlVAR.refit$beta.mat
                                   ,intercept = intercept
                                   ,include.obs.mat = include.obs.mat
                                   ,include.cf.mat = include.cf.mat
                                   ,cost.type = "BIC"
                                   ,BIC.parm = BIC.parm))
      #++++++++++++++++++++++++++++++++++++
      #cat(i,",",k,": ",   BIC.lc[i,k],"\n")
      #++++++++++++++++++++++++++++++++++++
      
    }#cRatios
    
  }#lambda
  
  
  
  BIC.min <- min(BIC.lc)
  BIC.min.index <- which(BIC.lc == BIC.min,arr.ind = FALSE)[1] # Gets the first if there are multiple
  col.index <- ceiling(BIC.min.index/lmbs.len)
  row.index <- BIC.min.index - (col.index - 1)* lmbs.len
  
  lambda.min <- lmbs[row.index]
  cRatio.min <- cRatios[col.index]
  
  #bound.1se <- get.withinSD.indices(mat = rss.lc,sd.mat = rss.lc.se)
  #rss.min.1se.element <- max(which(bound.1se))
  #rss.min.1se.index <- c(rss.min.1se.element %% lmbs.len + (rss.min.1se.element %% lmbs.len ==0)*lmbs.len, rss.min.1se.element / lmbs.len + 1)
  #rss.min.1se <- rss.lc[rss.min.1se.index]
  #lambda.min.1se <- lmbs[rss.min.1se.index[1]]
  #cRatio.min.1se <- cRatios[rss.min.1se.index[2]] 
  
  tune.profile <- list(pen = BIC.lc
                       #,pen.se = BIC.lc.se
                       ,lambda.min = lambda.min
                       #,lambda.min.1se = lambda.min.1se
                       ,cRatio.min = cRatio.min
                       #,cRatio.min.1se = cRatio.min.1se
                       ,tune.method = "BIC"
                       #,nlVAR.refit = nlVAR.refit#####++++++++++++++delete nlVAR.refit
                       )
  
  return(tune.profile)
  
}