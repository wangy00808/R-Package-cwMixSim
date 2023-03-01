# transformation of Yeo and Johnson, Biometrika (2000)
# this function transforms vector y

YJ.trans <- function(y, lambda){
  
  ind1 <- y < 0
  ind2 <- lambda == 0
  ind3 <- lambda == 2
  
  index <- ind1 & ind3
  if (sum(index) != 0) y[index] <- -log(1 - y[index])
  index <- !ind1 & ind2
  if (sum(index) != 0) y[index] <- log(1 + y[index])
  index <- ind1 & !ind3
  if (sum(index) != 0) y[index] <- (1 - (1 - y[index])^(2-lambda)) / (2 - lambda)
  index <- !ind1 & !ind2
  if (sum(index) != 0) y[index] <- ((y[index] + 1)^lambda - 1) / lambda
  
  return(y)
  
}



# inverse Yeo-Johnson transformation (from normality to skewness)
# provides inverse transformation for vector Ty
# the transform has two conditions

YJ.invtrans <- function(Ty, lambda){
  
  ind1 <- Ty < 0
  ind2 <- lambda == 0
  ind3 <- lambda == 2
  
  if (sum((Ty > 0) & (lambda <= -1 / Ty)) > 0){
    stop("Issue with inverse transformation...\n")
  }
  if (sum((Ty < 0) & (lambda >= 2 - 1 / Ty)) > 0){
    stop("Issue with inverse transformation...\n")
  }    
  
  index <- ind1 & ind3
  if (sum(index) != 0) Ty[index] <- 1 - exp(-Ty[index])
  index <- !ind1 & ind2
  if (sum(index) != 0) Ty[index] <- exp(Ty[index]) - 1
  index <- ind1 & !ind3
  if (sum(index) != 0) Ty[index] <- 1 - (1 - Ty[index] * (2 - lambda))^(1 / (2 - lambda))
  index <- !ind1 & !ind2
  if (sum(index) != 0) Ty[index] <- (Ty[index] * lambda + 1)^(1 / lambda) - 1
  
  return(Ty)
  
}





# function for calculating skewness as defined in Mardia

YJ.skew <- function(lambda, TY)
{
  p = dim(TY[[1]])[2]
  skew = NULL
  K <- length(TY)
  for(j in 1:K)
  {
    Y <- TY[[j]]
    for( i in c(1:p))
       {
          Y[,i] <- YJ.invtrans(TY[[j]][,i], lambda[i])
       }

    skew[j] <-  mardia(Y, na.rm = TRUE, plot = FALSE)$b1p
  }

  return(skew)
}




# function for finding roots
# m - multiplier variable

skew.root <- function(m, lambda, Gamma0, Ifbar, TY)
{
  
  gamma <- YJ.skew((m*lambda), TY)
  #cat(max(gamma), "\n")
  if(Ifbar) {return(mean(gamma) - Gamma0)} else {return(max(gamma) - Gamma0)}
  
}





# Simulate skewed mixture given levels of skewness and overlap
#library(mvtnorm)
#library(MixSim)
#library(rootSolve)
#library(psych)
#library(doParallel)
#library(foreach)

MixSim.skew = function(BarOmega = NULL, MaxOmega = NULL, BarGamma = NULL, MaxGamma = NULL, K, p, sph = FALSE, hom = FALSE,
                    ecc = 0.90, PiLow = 1.0, int = c(0.0, 1.0), nM = 1000, resN = 100, tol=1e-06, eps = 1e-06, lim = 1e06){
   
  if(is.null(BarGamma)&&is.null(MaxGamma))
     {BarGamma = 0
     MaxGamma = NULL
   }
  if((!is.null(BarGamma)) && (!is.null(MaxGamma))) stop("BarGamma and MaxGamma can't be both specified...\n")
   
  if(is.null(BarGamma)) 
  {Gamma = MaxGamma
  Ifbar = FALSE
  }else{
    Gamma = BarGamma
    Ifbar = TRUE
  }
  iter = 0
  n = nM
  
  
  ex = MixSim(BarOmega = BarOmega, MaxOmega = MaxOmega, K = K, p = p, sph = sph, hom = hom, ecc = ecc, PiLow = PiLow, 
              int = int, resN = resN, eps = eps, lim = lim)
  if(is.null(ex)) stop("...\n")
  Mu = ex$Mu
  S  = ex$S
  pi = ex$Pi
  TY1 = simdataset(n, pi, Mu, S)
  TY = NULL
  for(k in c(1:(nrow(Mu)))) TY = append(TY, list(TY1$X[which(TY1$id ==k), ]))
  
  if(Gamma!=0)
  { repeat
     {
      iter = iter +1
      lambda=runif(p,0,1)
      m=tryCatch(uniroot(f = skew.root, lambda = lambda, Gamma0 = Gamma, Ifbar = Ifbar, TY = TY, interval = c(0,4), extendInt = "yes", tol = tol)
                          , error=function(e){})
      if(!is.null(m)) {break}
      if(iter > resN) {break}
              #print(iter)
     }
  }
   
  if(Gamma == 0){
    new.lambda = numeric(K)
    gamma = 0
  }
  
  if(Gamma!=0){   
     if(is.null(m)) {
         stop("the desired skewness has not been reached, try a different seed or change the value for BarGamma or MaxGamma\n")
     }else{
       new.lambda=lambda*(m$root)
       #Y = TY
       #for(k in c(1:K))
       #  {for( i in c(1:p))
        #    { Y[[k]][,i] <- YJ.invtrans(TY[[k]][,i], new.lambda[i])
        #    }
       #}
       gamma = YJ.skew(new.lambda, TY)
     }
  }     

   
   
  return(list(lambda = new.lambda, Pi=ex$Pi, Mu=Mu, S2=S, OmegaMap=ex$OmegaMap, BarOmega = ex$BarOmega, MaxOmega = ex$MaxOmega, rcMax = ex$rcMax, Gamma = gamma, BarGamma=mean(gamma), MaxGamma= max(gamma)))
}



summary.MixSim.skew <- function(object, ...){
  K <- length(object$Pi)
  OmegaMap <- object$OmegaMap
  colnames(OmegaMap) <- paste("k.", 1:K, sep = "")
  rownames(OmegaMap) <- paste("k.", 1:K, sep = "")
  cat("OmegaMap: \n")
  print(OmegaMap)
  cat("\nrcMax:", object$rcMax, "\n")
  invisible()
}



# pdf for original data based on Yeo-Johnson transformation
# y is a vector

YJ.pdf <- function(y, lambda, Mu, S){
  
  p = length(lambda)
  Ty = y
  Jacob = y
  if(!is.null(dim(y))&!is.null(lambda))
  {
    for( i in c(1:p))
    {
      Ty[,i] <- YJ.trans(y[,i], lambda[i])
      ind1 <- y[,i] < 0
      ind2 <- lambda[i] == 0
      ind3 <- lambda[i] == 2
      index <- (ind1 & ind3)|(!ind1 & ind2)
      if (sum(index) != 0) Jacob[index, i] <- 1/(1 + abs(y[index,i]))
      
      index <- (ind1 & !ind3)|(!ind1 & !ind2)
      if (sum(index) != 0) Jacob[index, i] <- (1 + abs(y[index,i]))^(sign(y[index,i]) * (lambda[i]-1))
      
    }
    J=apply(Jacob, 1, prod)
  }
  
  if(is.null(dim(y))&!is.null(lambda)){
    
    for( i in c(1:p))
    {
      Ty[i] <- YJ.trans(y[i], lambda[i])
      ind1 <- y[i] < 0
      ind2 <- lambda[i] == 0
      ind3 <- lambda[i] == 2
      index <- (ind1 & ind3)|(!ind1 & ind2)
      if (sum(index) != 0) Jacob[i] <- 1/(1 + abs(y[i]))
      
      index <- (ind1 & !ind3)|(!ind1 & !ind2)
      if (sum(index) != 0) Jacob[i] <- (1 + abs(y[i]))^(sign(y[i]) * (lambda[i]-1))
      
    }
    
    J=prod(Jacob)
  }
  
  if(is.null(lambda)) J=1
  
  if(max(lambda) == 0){
    res <- dmvnorm(y, Mu, S) 
  }else{
    res <- dmvnorm(Ty, Mu, S) * abs(J)
  }
  
  return(res)
}


######### simulate skew data given object

simdataset.skew= function(n, Pi, Mu, S2, lambda){
  
  TY1 = simdataset(n, Pi, Mu, S2)
  K = length(Pi)
  p = ncol(TY1$X)
  
  TY = NULL
  for(k in c(1:(nrow(Mu)))) {TY = append(TY, list(TY1$X[which(TY1$id ==k), ]))}
  Y = TY
  
  for(k in c(1:K))
   {for( i in c(1:p))
     { Y[[k]][,i] <- YJ.invtrans(TY[[k]][,i], lambda[i])
     }
  }
  

  Y=do.call(rbind,Y)
  if(max(lambda)==0) Y = TY1$X
  data = list(X = Y, id = TY1$id)
  return(data)
}


### pdf of the skewed mixture
pdf.skew = function(y, lambda, Mu, S, Pi)
{
  k = dim(Mu)[1]
  pdf = sapply(c(1:k), function(x) YJ.pdf(y, lambda, Mu[x,], S[,,x]))
  
  return(rowSums(pdf*Pi))
}


### pdf of the skewed mixture for the contour plot
contour.skew = function(x, y, opt)
{   
    i = NULL
    n = length(y)
    
    #cl<-makePSOCKcluster(c(rep("localhost",4)),outfile='',homogeneous=FALSE,port=11001)
    #f = foreach(i = 1: n, .combine = cbind) %dopar% { pdf.skew(cbind(x, rep(y[i], n)), Pi = opt$Pi, Mu = opt$Mu, S = opt$S, lambda = opt$lambda)}
    #stopCluster(cl)
    
    cl = makeCluster(4)
    clusterSetRNGStream(cl, 123)
    clusterExport(cl, varlist=c("n", "x", "y", "opt", "pdf.skew","YJ.pdf","YJ.trans","dmvnorm"), envir=environment())
    f = parLapply(cl, 1:n, function(i) {pdf.skew(cbind(x, y = rep(y[i], n)), Pi = opt$Pi, Mu = opt$Mu, S = opt$S2, lambda = opt$lambda)})
    stopCluster(cl)
    f = matrix(as.numeric(unlist(f)), n, n)
    
    
    
    return(f)
}


### 2d contour plot of the skewed mixture
contour2d.skew = function(opt, xlim =NULL, ylim = NULL, mar = c(3.5,3.5,3.5,3.5), nlevels= 100, lwd = 0.2, col = "red", yaxt = NULL, xaxt =NULL){
  
  #data = opt$data
  
  
  if(is.null(xlim)||is.null(ylim)) {
    data = simdataset.skew(n=1000, Pi = opt$Pi, Mu = opt$Mu, S2= opt$S2, lambda = opt$lambda)$X
    colmax = apply(data,2,max)+0.1
    colmin = apply(data,2,min)-0.1}
  
  if(is.null(xlim)) 
  {
    xlim[1] = colmin[1]
    xlim[2] = colmax[1]
  }
  if(is.null(ylim)) 
  {ylim[1] = colmin[2]
  ylim[2] = colmax[2] 
  }
  
  x = seq(xlim[1], xlim[2], length.out = 500)
  y = seq(ylim[1], ylim[2], length.out = 500)
  z = contour.skew(x, y, opt)
  par(mar=mar)
  contour(x, y, z, nlevels =  nlevels, lwd = lwd, col = col, drawlabels = FALSE, yaxt = yaxt, xaxt = xaxt )
  #centers = cbind(YJ.invtrans(opt$Mu[,1], opt$lambda[1]), YJ.invtrans(opt$Mu[,2], opt$lambda[2]))
  #points(centers, pch = 20, lwd = 0.5)
  
}




