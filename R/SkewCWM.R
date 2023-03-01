#library(mvtnorm)
#library(MixSim)
#library(rootSolve)
#library(psych)

#library(doParallel)

#library(bazar)

MixSim.skewcwm = function(BarOmega=NULL, MaxOmega = NULL, BarGamma = NULL, MaxGamma = NULL,K, p, Pi = NULL, Mu.x_lim = c(-5, 5), Beta_lim = c(-5,5), S2.y_lim = c(0,1), S2.x_lim = c(0,1),N = 500, nM = 200, resN = 100, tol=1e-06)
{
  Parallel = FALSE
  n.cores = 4
  if((!is.null(BarGamma)) && (!is.null(MaxGamma))) stop("BarGamma and MaxGamma can't be both specified...\n")
  if(is.null(BarGamma)) 
  {Gamma0 = MaxGamma
  Ifbar = FALSE
  }else{
    Gamma0 = BarGamma
    Ifbar = TRUE
  }

  if(is.null(BarGamma)&&is.null(MaxGamma)) {
    Gamma0 = 0
  }

  
  ex =  MixSim.cwm(BarOmega = BarOmega, MaxOmega = MaxOmega, K = K, p = p, Pi = Pi, Mu.x_lim = Mu.x_lim, Beta_lim = Beta_lim, S2.y_lim = S2.y_lim, S2.x_lim = S2.x_lim,N =N, resN = resN, tol = tol)
  Mu.x = ex$Mu.x
  Beta = ex$Beta
  S2.x  = ex$S2.x
  S2.y = ex$S2.y
  data1 = simdataset.cwm(n =nM, Pi = ex$Pi, Mu.x = ex$Mu.x, S2.x = ex$S2.x, Beta = ex$Beta, S2.y = ex$S2.y)
  X = data1$X
  id = data1$id
  
  TY = list(c(data1$Y[id==1]))
  for(k in 2:K) TY = append(TY, list(c(data1$Y[id==k])))
  Pi = ex$Pi
  f = NULL

  niter = 0
  if(Gamma0!=0 ){
   repeat
    { niter = niter+1
      lambda=runif(1,0,1)
      m=tryCatch(uniroot(f = skew.root.cwm, lambda = lambda, Gamma0 = Gamma0,  Ifbar = Ifbar, TY = TY, interval = c(0.1, 1), extendInt = "yes", tol = tol)
                 , error=function(e){
                 })
      if(!is.null(m)) {
        new.lambda=lambda*(m$root)
        Gamma = YJ.skew.cwm(new.lambda, TY)
        break}
      if(niter >resN) break
    }
    
    if(is.null(m)) stop("Error: the desired overlap has not been reached in reN simulations...\n Increase the number of simulations allowed (option resN) or change the value of overlap or skewness...\n")
    #Y = TY
    #for(j in c(1:K)) Y[[j]] <- YJ.invtrans(TY[[j]], new.lambda)
    }else{
    #Y = TY
    new.lambda = 1
    Gamma = rep(0,K)
  }
  
  
  
  return(list(lambda = new.lambda, Pi = Pi, Mu.x = Mu.x, S2.x = S2.x, Beta = Beta, S2.y = S2.y, OmegaMap=ex$OmegaMap, BarOmega = ex$BarOmega, MaxOmega = ex$MaxOmega, rcMax = ex$rcMax, Gamma = Gamma, BarGamma=mean(Gamma), MaxGamma= max(Gamma)))
}



# function for calculating skewness as defined in Mardia

YJ.skew.cwm <- function(lambda, TY)
{

  k = length(TY)
  Gamma = NULL
  Y <- TY
  for( j in c(1:k))
  { Y[[j]] <- YJ.invtrans(TY[[j]], lambda)
    Gamma[j]<-  mardia(Y[[j]],na.rm = TRUE, plot = FALSE)$b1p 
  }

  return(Gamma)
}


# function for finding roots
# m - multiplier variable

skew.root.cwm <- function(m, lambda, Gamma0, Ifbar, TY)
{
  
  gamma <- YJ.skew.cwm((m*lambda), TY)
  #if(Ifbar) cat(mean(gamma), "\n") else cat(max(gamma), "\n")
  if(Ifbar) {return(mean(gamma) - Gamma0)} else {return(max(gamma) - Gamma0)}
  
}

simdataset.skewcwm = function(n, Pi, Mu.x, S2.x, Beta, S2.y,lambda)
{
  K = length(Pi)
  data = simdataset.cwm(n, Pi = Pi, Mu.x = Mu.x, S2.x = S2.x, Beta = Beta, S2.y = S2.y)
  Y = data$Y
  id = data$id
  for(j in c(1:K)) Y[id ==j] <- YJ.invtrans(data$Y[id ==j], lambda)
  # data
  data = list(X = data$X, Y = Y, id = id)
  
  return(data)
}

# plot =====================================================================
# pdf for original data based on Yeo-Johnson transformation
# y is a vector
# transformation of Yeo and Johnson, Biometrika (2000)
# this function transforms vector y

#YJ.trans

Jacob <- function(y, lambda)
{
  J = y
  k = ncol(y)
  if(!is.null(lambda))
  {
      ind1 <- y < 0
      ind2 <- lambda == 0
      ind3 <- lambda == 2
      index <- (ind1 & ind3)|(!ind1 & ind2)
      if (sum(index) != 0) J[index] <- 1/(1 + abs(y[index]))
      index <- (ind1 & !ind3)|(!ind1 & !ind2)
      if (sum(index) != 0) J[index] <- (1 + abs(y[index]))^(sign(y[index]) * (lambda-1))
  }
  if(is.null(lambda)) J = 1
  return(J)
}

####### simulate skewed cwm data=============================



#####plot=====================================================
# find the pdf of mixture of regression with point input
pdf.skewcwm = function(x, y, Pi, Mu.x, S2.x, Beta, S2.y, lambda)
{
  k = length(Pi)
  ty = YJ.trans(y, lambda)
  J = Jacob(y, lambda)
  
  if(!is.matrix(x)){
    f = sapply(c(1:k), function(i) Pi[i]*dnorm(ty, cbind(1,x)%*%Beta[i,], S2.y[i]^0.5)*dnorm(x, Mu.x[i,], S2.x[i]^0.5)*abs(J) )
  }else{
    f = sapply(c(1:k), function(i) Pi[i]*dnorm(ty, cbind(1,x)%*%Beta[i,], S2.y[i]^0.5)*dmvnorm(x, Mu.x[i,], S2.x[,,i])*abs(J) )
  }
  return(rowSums(f))
}




contour_skew.cwm = function(x,y,opt)
{
  i = NULL
  N = length(y)
  #cl<-makePSOCKcluster(c(rep("localhost",4)),outfile='',homogeneous=FALSE,port=11001)
  #f=foreach(i = 1: n, .combine = cbind) %dopar% { pdf.skewcwm (x, y = rep(y[i], n), Pi = opt$Pi, Mu.x = opt$Mu.x, S2.x = opt$S2.x, Beta = opt$Beta, S2.y = opt$S2.y, lambda = opt$lambda)}
  #stopCluster(cl)
  
  cl = makeCluster(4)
  clusterSetRNGStream(cl, 123)
  clusterExport(cl, varlist=c("N", "x", "y", "opt", "pdf.skewcwm","YJ.trans","Jacob"), envir=environment())
  f = parLapply(cl, 1:N, function(i) {pdf.skewcwm(x, y = rep(y[i], N), Pi = opt$Pi, Mu.x = opt$Mu.x, S2.x = opt$S2.x, Beta = opt$Beta, S2.y = opt$S2.y, lambda = opt$lambda)})
  stopCluster(cl)
  
  f = matrix(as.numeric(unlist(f)), N, N)
  
  
  return(f)
}



contour2d.skewcwm = function(opt, xlim =NULL, ylim = NULL, mar = c(3.5,3.5,3.5,3.5), nlevels= 100, lwd = 0.2,  col = "red", yaxt = NULL, xaxt = NULL){
  
  
  if(is.null(xlim)||is.null(ylim)) {
    data = simdataset.skewcwm(n=200, Pi = opt$Pi, Mu.x = opt$Mu.x, S2.x = opt$S2.x, Beta = opt$Beta, S2.y = opt$S2.y, lambda = opt$lambda)
    data = cbind(data$X, data$Y)
    colmax = apply(data,2,max)+10
    colmin = apply(data,2,min)-10}
  
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
  z = contour_skew.cwm(x, y, opt)
  par(mar=mar)
  contour(x, y, z, nlevels =  nlevels, lwd = lwd, col = col, drawlabels = FALSE, yaxt = yaxt, xaxt = xaxt )
  #centers = cbind(YJ.invtrans(opt$Mu[,1], opt$new.lambda[1]), YJ.invtrans(opt$Mu[,2], opt$new.lambda[2]))
  #points(centers, pch = 20, lwd = 0.5)
  
}




