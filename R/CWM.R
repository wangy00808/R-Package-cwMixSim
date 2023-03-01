
#library(mvtnorm)
#library(MixSim)
#library(rootSolve)
#library(psych)
#library(doParallel)
#library(foreach)
#library(bazar)
#library(parallel)
#library(graphics)

#Generate a positive definite matrix
Posdef = function(N, ev = runif(N, 0, 10)) 
{ 
  Z <- matrix(ncol=N, rnorm(N^2)) 
  decomp <- qr(Z) 
  Q <- qr.Q(decomp) 
  R <- qr.R(decomp) 
  d <- diag(R) 
  ph <- d / abs(d) 
  O <- Q %*% diag(ph) 
  Z <- t(O) %*% diag(ev) %*% O 
  return(Z) 
} 


# measure the overlap of the CWM model
overlap.cwm = function(Pi, S2.y, N, NN, dX, Mu_y, Parallel, n.cores)
{
  K = length(Pi)
  OmegaMap = diag(1, K, K)
  if(min(S2.y)>1e-04 && max(S2.y)<1e+04)
  {  for(i in 1:K)
     {
     if(Parallel){
       cl = makeCluster(n.cores)
       clusterSetRNGStream(cl, 123)
       clusterExport(cl, varlist=c("i", "Pi", "dX", "Mu_y", "K", "S2.y", "overlap"), envir=environment())
       result = parLapply(cl, 1:NN[i], function(j) {overlap(Pi*dX[[i]][j,], matrix(Mu_y[[i]][j, ],K,1), array(S2.y, c(1, 1, K)))$OmegaMap})
       stopCluster(cl)
       result = array(as.numeric(unlist(result)), dim=c(K, K, NN[i]))
     }else{
       result = array(0, dim = c(K,K,NN[i]))
       for(j in 1:NN[i]) result[,,j] = overlap(Pi*dX[[i]][j,], matrix(Mu_y[[i]][j, ],K,1), array(S2.y, c(1, 1, K)))$OmegaMap
     }
     OmegaMap[i, ] = apply(result, c(1,2), mean)[i, ]
  }
  }
  
  BarOmega = (sum(OmegaMap) - sum(diag(OmegaMap)))/((K^2-K)/2)
  MaxOmega = max((OmegaMap +t(OmegaMap))-diag(2,K))
  #print(BarOmega)
  return(list(OmegaMap = OmegaMap, BarOmega = BarOmega, MaxOmega = MaxOmega))
}

#overlap.cwm = function(Pi, S2.y, N, NN, dX, Mu_y){
#  K = length(Pi)
#  OmegaMap = diag(K)
#  #Mu_y = do.call(rbind,Mu_y)
#  #dX = do.call(rbind, dX)
#  for(k1 in 1:(K-1)){
#    for(k2 in ((k1+1):K)){
#         if(S2.y[k1] == S2.y[k2]){
#          phi_x = S2.y[k1]^0.5/(2*(Mu_y[[k1]][,k1] - Mu_y[[k1]][,k2]))*(-(Mu_y[[k1]][,k1] - Mu_y[[k1]][,k2])^2/S2.y[[k1]][k1]+log((Pi[k2]^2*dX[[k1]][,k2]^2)/(Pi[k1]^2*dX[[k1]][,k1]^2)))
#          misp = pnorm(phi_x)
#          OmegaMap[k1,k2] = mean(ifelse(Mu_y[[k1]][,k1] > Mu_y[[k1]][,k2], misp, 1-misp))
#          
#          phi_x = S2.y[k2]^0.5/(2*(Mu_y[[k2]][,k2] - Mu_y[[k2]][,k1]))*(-(Mu_y[[k2]][,k2] - Mu_y[[k2]][,k1])^2/S2.y[k2]+log((Pi[k1]^2*dX[[k2]][,k1]^2)/(Pi[k2]^2*dX[[k2]][,k2]^2)))
#          misp = pnorm(phi_x)
#          OmegaMap[k2,k1] = mean(ifelse(Mu_y[[k2]][,k2] > Mu_y[[k2]][,k1], misp, 1-misp))
#          }
      
#        if(S2.y[k1] != S2.y[k2]){
#          chisq_x = (S2.y[k2]/(S2.y[k1]-S2.y[k2]))* ((Mu_y[[k1]][,k1] - Mu_y[[k1]][,k2])^2/(S2.y[k1] - S2.y[k2]) + log((Pi[k2]^2*dX[[k1]][,k2]^2*S2.y[k1])/(Pi[k1]^2*dX[[k1]][,k1]^2*S2.y[k2])))
#          chisq_c = (Mu_y[[k1]][,k1] - Mu_y[[k1]][,k2])^2*(S2.y[k1]/(S2.y[k1]-S2.y[k2])^2)
#          OmegaMap[k1,k2] = mean(pchisq(chisq_x, df = 1, ncp = chisq_c))
#         if(S2.y[k1] < S2.y[k2])  OmegaMap[k1,k2] = 1-  OmegaMap[k1,k2]
         
#         chisq_x = (S2.y[k1]/(S2.y[k2]-S2.y[k1]))* ((Mu_y[[k2]][,k2] - Mu_y[[k2]][,k1])^2/(S2.y[k2] - S2.y[k1]) + log((Pi[k1]^2*dX[[k2]][,k1]^2*S2.y[k2])/(Pi[k2]^2*dX[[k2]][,k2]^2*S2.y[k1])))
#         chisq_c = (Mu_y[[k2]][,k2] - Mu_y[[k2]][,k1])^2*(S2.y[k2]/(S2.y[k2]-S2.y[k1])^2)
#         OmegaMap[k2,k1] = mean(pchisq(chisq_x, df = 1, ncp = chisq_c))
#         if(S2.y[k2] < S2.y[k1])  OmegaMap[k2,k1] = 1-  OmegaMap[k2,k1]
#        }
#    }
#  }
#  BarOmega = (sum(OmegaMap) - sum(diag(OmegaMap)))/((K^2-K)/2)
#  MaxOmega = max((OmegaMap +t(OmegaMap))-diag(2,K))
 
#  return(list(OmegaMap = OmegaMap, BarOmega = BarOmega, MaxOmega = MaxOmega))
#}


reg.root1 = function(a, MaxOmega, Pi, S2.y, N, NN,dX, Mu_y, Parallel = TRUE, n.cores = 4)
{
  f = overlap.cwm(Pi,  a*S2.y, N, NN,dX, Mu_y, Parallel = Parallel, n.cores = 4)$MaxOmega - MaxOmega
  #cat("a is ",a, "\n")
  return(f)
}


reg.root2 = function(b, ij, BarOmega, Pi, S2.y, N, NN, dX, Mu_y, Parallel = TRUE, n.cores = 4)
{
  if(!is.null(ij)) S2.y[-ij] = b*S2.y[-ij] else S2.y = b*S2.y
  #cat("b is ",b, "\n")
  f = overlap.cwm(Pi, S2.y, N, NN, dX, Mu_y, Parallel = Parallel, n.cores = 4)$BarOmega - BarOmega
  return(f)
}

# the initialization 
intlz = function(K, p, NN, Mu.x_lim, Beta_lim, S2.y_lim, S2.x_lim){

  p_x = p - 1
  # the mean and variance of x
  Mu.x = matrix(runif(p_x*K, Mu.x_lim[1], Mu.x_lim[2]), K, p_x)
  if(p_x==1) {S2.x =runif(K, 0, S2.x_lim[2])
  }else {S2.x =replicate(K, Posdef(p_x, ev = runif(p_x, 0, S2.x_lim[2])))}
  # the variance and coefficients of y
  S2.y = runif(K, 0, S2.y_lim[2])
  Beta = matrix(runif(K*(p_x+1), Beta_lim[1], Beta_lim[2]), K, p_x+1)
  
  # simulate X, X is a list
  X = list(matrix(0,NN[1],p_x))
  for(i in 1:(K-1)) X = append(X, list(matrix(0,NN[i+1],p_x)))
  
  if(p_x > 1) {X = lapply(1:K, function(i) rmvnorm(NN[i], Mu.x[i, ], S2.x[,,i]))
   }else {X = lapply(1:K, function(i) rnorm(NN[i], Mu.x[i, ], S2.x[i]^0.5))}
  
  # calculate the mean of y given each x, and the densities of x
  dX = list(matrix(0,NN[1],K))
  for(k in 2:K)  dX = append(dX, list(matrix(0,NN[k],K)))
  Mu_y = dX
  for(k in 1:K)
  { if(p_x > 1)
    {dX[[k]] = sapply(1:K, function(j) dmvnorm(X[[k]], Mu.x[j,], S2.x[,,j]))
    Mu_y[[k]] = sapply(1:K, function(j) cbind(1, X[[k]])%*%Beta[j, ])
    }else{
      dX[[k]] = sapply(1:K, function(j) dnorm(X[[k]], c(Mu.x)[j], c(S2.x)[j]^0.5))
      Mu_y[[k]] = sapply(1:K, function(j) cbind(1, X[[k]])%*%Beta[j, ])
    }
  }
  #simulate y
  Y = lapply(1:K, function(j) rnorm(NN[j], Mu_y[[j]][,j], S2.y[j]^0.5))  
  
  return(list(X = X, Y = Y, Mu.x = Mu.x, S2.x = S2.x, Beta = Beta, S2.y = S2.y, dX = dX, Mu_y = Mu_y))
}



# main function: simulate CWM given barOmega and MaxOmega===========
MixSim.cwm = function(BarOmega = NULL, MaxOmega = NULL, K, p, Pi = NULL, Mu.x_lim = c(-5, 5), Beta_lim = c(-5,5), S2.y_lim = c(0,1), S2.x_lim = c(0,1), N = 500, resN = 10,  tol=1e-06){
 
   #options(warn=-1)
  Parallel = FALSE
  n.cores = 4
  if(!(is.null(BarOmega)||is.null(MaxOmega))) if(K==2&&(BarOmega!=MaxOmega)) stop("BarOmega and MaxOmega should equal for a two-dimension data")
 # if(BarOmega>1|MaxOmega>1) stop("average and maximum overlaps are values between 0 to 1")
  if(!is.null(Pi) && ((sum(Pi)!= 1|length(Pi)>K))) stop("incorrect values of mixing proportions")
  if(!(is.null(BarOmega)||is.null(MaxOmega))) if(MaxOmega < BarOmega | MaxOmega > (BarOmega * K *(K - 1) / 2)) stop("Error: incorrect values of average and maximum overlaps for k>2...
  Both conditions should hold:
  1. MaxOverlap > AverOverlap
  2. MaxOverlap < AverOverlap * K (K - 1) / 2")
  if(!(is.null(BarOmega)||is.null(MaxOmega))) if(K>2&&(BarOmega==MaxOmega))  stop("Error: incorrect values of average and maximum overlaps for k>2...
  Both conditions should hold:
  1. MaxOverlap > AverOverlap
  2. MaxOverlap < AverOverlap * K (K - 1) / 2")
  

  if(is.null(Pi)) 
  {Pi = rep(1/K, K)
   Pi = round(Pi,round(log10(N)))
   Pi[K] = 1 - sum(Pi[1:(K-1)])}
  # the number of observation in each category
  NN = Pi*N
  NN = round(NN,0)
  NN[K] = N - sum(NN[1:(K-1)])
 
  ii = 0
  prsn = 0.005
  
  repeat
  { ii = ii + 1
     #initalization
     opt1 = intlz(K,p,NN, Mu.x_lim, Beta_lim, S2.y_lim, S2.x_lim)
     S2.y = opt1$S2.y
     # solve sigma when the max_omega meet
     m1 = NULL
     ij = NULL
     if(!is.null(MaxOmega)) 
       {m1 = tryCatch(uniroot(f = reg.root1, MaxOmega = MaxOmega, Pi = Pi, S2.y = S2.y, 
                         N = N, NN = NN, dX = opt1$dX, Mu_y = opt1$Mu_y, Parallel = Parallel, n.cores = n.cores, interval = c(1e-6, 1), extendInt = "yes",  tol = tol), error=function(e){})
     }
     prsn1 = 1
       # if m1 is not Null
     if(!is.null(m1)){
         #update the variances of Y
         S2.y = m1$root*S2.y
         f = overlap.cwm(Pi, S2.y = S2.y, N, NN, dX = opt1$dX, Mu_y = opt1$Mu_y, Parallel = Parallel, n.cores = n.cores)
         if(MaxOmega!=0) prsn1 = abs((f$MaxOmega-MaxOmega)/MaxOmega)
         if(MaxOmega ==0) prsn1 = abs((f$MaxOmega-MaxOmega)/0.1)
         # Find the pair of clusters which has the largest overlap
         ij = which(f$OmegaMap+t(f$OmegaMap)==f$MaxOmega, TRUE)[1,]
     }
     # solve sigma when the bar_omega meet
     m2 = NULL
     if((K > 2&&!is.null(m1)&&!is.null(MaxOmega)&&!is.null(BarOmega)&&prsn1<prsn)||(is.null(MaxOmega)&&!is.null(BarOmega))){
       m2 = tryCatch(uniroot(f = reg.root2,ij, BarOmega = BarOmega, Pi = Pi, S2.y = S2.y, 
                             N = N, NN = NN,  dX = opt1$dX, Mu_y = opt1$Mu_y,  Parallel = Parallel, n.cores = n.cores, interval = c(1e-6, 1), extendInt = "yes", tol = tol), error=function(e){})
       if(!is.null(m2)){
          if(!is.null(ij)) S2.y[-ij] = m2$root*S2.y[-ij] else S2.y = m2$root*S2.y
          f = overlap.cwm(Pi, S2.y = S2.y, N, NN, dX = opt1$dX, Mu_y = opt1$Mu_y,  Parallel = Parallel, n.cores = n.cores)
          if(BarOmega!=0) if(is.null(MaxOmega)) prsn1 = abs((f$BarOmega-BarOmega)/BarOmega) else prsn1 = max(abs((f$BarOmega-BarOmega)/BarOmega), prsn1)
          if(BarOmega==0) if(is.null(MaxOmega)) prsn1 = abs((f$BarOmega-BarOmega)/0.1) else prsn1 = max(abs((f$BarOmega-BarOmega)/0.1), prsn1)
          ij = which(f$OmegaMap+t(f$OmegaMap)==f$MaxOmega, TRUE)[1,]
       }
     }

     if(!is.null(MaxOmega)&&!is.null(BarOmega))
         { if(K == 2 &&!is.null(m1)&&prsn1<prsn) break
           if(K > 2 &&!is.null(m2)&&prsn1<prsn) break}
     if(is.null(MaxOmega)&&!is.null(BarOmega))  if(!is.null(m2)&&prsn1<prsn) break
     if(!is.null(MaxOmega)&&is.null(BarOmega))  if(!is.null(m1)&&prsn1<prsn) break
     if(is.null(MaxOmega)&&is.null(BarOmega)) break
     #if(abs(f$MaxOmega - MaxOmega)/MaxOmega < eps )  break
     if(ii == resN)    break
  }
  
  
  if(!is.null(MaxOmega)&&is.null(BarOmega)&&(is.null(m1)||prsn1>prsn)) stop("Error: the desired overlap has not been reached in reN simulations...\n Increase the number of simulations allowed (option resN) or change the value of overlap...\n")
  if(is.null(MaxOmega)&&!is.null(BarOmega)&&(is.null(m2)||prsn1>prsn)) stop("Error: the desired overlap has not been reached in reN simulations...\n Increase the number of simulations allowed (option resN) or change the value of overlap...\n")
  if(!is.null(BarOmega)&&!is.null(MaxOmega)) if(is.null(m1)||(K>2&&is.null(m2))||prsn1>prsn) stop("Error: the desired overlap has not been reached in reN simulations...\n Increase the number of simulations allowed (option resN) or change the value of overlap...\n")
  if(is.null(MaxOmega)&&is.null(BarOmega)){ 
     f = overlap.cwm(Pi, S2.y = S2.y, N, NN, dX = opt1$dX, Mu_y = opt1$Mu_y, Parallel = Parallel, n.cores = n.cores)
     ij = which(f$OmegaMap+t(f$OmegaMap)==f$MaxOmega, TRUE)[1,]}

  
  return(list(Pi = Pi, Mu.x = opt1$Mu.x, S2.x = opt1$S2.x, Beta = opt1$Beta, S2.y = S2.y, OmegaMap=f$OmegaMap, BarOmega = f$BarOmega, MaxOmega = f$MaxOmega, rcMax = ij))
}



# simulate the data follows CWM distribution  ==============================================
simdataset.cwm = function(n, Pi, Mu.x, S2.x, Beta, S2.y)
{
  p_x = ncol(Mu.x)
  K = length(Pi)
  # the number of observations in each category
  
  NN = c(rmultinom(prob = Pi, size = n, 1))
  #simulate X
  X = list(matrix(0,NN[1],p_x))
  for(i in 1:(K-1)) X = append(X, list(matrix(0,NN[i+1],p_x)))

  if(p_x > 1) {X = lapply(1:K, function(i) rmvnorm(NN[i], Mu.x[i, ], S2.x[,,i]))
  }else {X = lapply(1:K, function(i) rnorm(NN[i], Mu.x[i, ], S2.x[i]^0.5))}
  
  #simulate y
  Mu_y = lapply(1:K, function(j) cbind(1, X[[j]])%*%Beta[j, ])
  Y = lapply(1:K, function(j) rnorm(NN[j], Mu_y[[j]], S2.y[j]^0.5))

  # data
  if(p_x == 1)
  { data = list(X = unlist(X), Y = unlist(Y), id = rep(1:K,NN))
  }else{
    data = list(X = do.call(rbind,X), Y = unlist(Y), id = rep(1:K,NN))
  }
  return(data)
}

#####plot========================================================================
# find the pdf of cluster weighted model
pdf.cwm = function(x, y, Pi, Mu.x, S2.x, Beta, S2.y)
{
  K = length(Pi)
  p_x = ncol(Mu.x)
  
  if(p_x == 1){
    f = sapply(c(1:K), function(i) Pi[i]*dnorm(y, cbind(1,x)%*%Beta[i,], S2.y[i]^0.5)*dnorm(x, Mu.x[i,], S2.x[i]^0.5) )
  }else{
    f = sapply(c(1:K), function(i) Pi[i]*dnorm(y, cbind(1,x)%*%Beta[i,], S2.y[i]^0.5)*dmvnorm(x, Mu.x[i,], S2.x[,,i]) )
  }
  
  return(rowSums(f))
}



# calculate the pdf for the contour plot 
contour.reg = function(x,y,opt)
{
  N = length(y)
  
  cl = makeCluster(4)
  clusterSetRNGStream(cl, 123)
  clusterExport(cl, varlist=c("N", "x", "y", "opt", "pdf.cwm"), envir=environment())
  f = parLapply(cl, 1:N, function(i) {pdf.cwm(x, y = rep(y[i], N), Pi = opt$Pi, Mu.x = opt$Mu.x, S2.x = opt$S2.x, Beta = opt$Beta, S2.y = opt$S2.y)})
  stopCluster(cl)
  f = matrix(as.numeric(unlist(f)), N, N)
  return(f)
}


# 2d contour plot for the CWM=======
contour2d.cwm = function(opt, xlim =NULL, ylim = NULL, mar = c(3.5,3.5,3.5,3.5), nlevels= 100, lwd = 0.2,  col = "red", yaxt = NULL, xaxt = NULL){
  
  #data = cbind(opt$X, opt$Y)
  #data = opt$data
 
  if(is.null(xlim)||is.null(ylim)) {
    data = simdataset.cwm(n=200, Pi = opt$Pi, Mu.x = opt$Mu.x, S2.x = opt$S2.x, Beta = opt$Beta, S2.y = opt$S2.y)
    data = cbind(data$X, data$Y)
    colmax = apply(data,2,max)+5
    colmin = apply(data,2,min)-5}
  
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
  
  z = contour.reg(x, y, opt)
  par(mar=mar)
  contour(x, y, z, nlevels =  nlevels, lwd = lwd, col = col, drawlabels = FALSE, yaxt = yaxt, xaxt = xaxt )
  #centers = cbind(YJ.invtrans(opt$mu[,1], opt$new.lambda[1]), YJ.invtrans(opt$mu[,2], opt$new.lambda[2]))
  #points(centers, pch = 20, lwd = 0.5)
  
}


