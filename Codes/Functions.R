#Note that for the following functions, it is assumed that
#X is an np x pd matrix produced from reshape_X()
#Y is an np-dimensional vector where Y = (Y_1^T,...,Y_n^T)^T,
#W is a p x p x (K+1) consisting of the set of similarity matrices {W_k: k=0,...,K}.

library(Matrix)
library(emulator)
library(matrixcalc)
library(MASS)

#To estimate the joint mean and correlation regression model
estimate_joint = function(X,Y,W,fam,init_beta,init_phi,init_alpha,stepsize=1,Lambda0=NA,mu=0.05,epsilon = 1e-5,xi = 1e-8,iterMAX = 1000, xi2 = 1e-8, iterMAX2=1000){
  p = dim(W)[1]; K = dim(W)[3] 
  
  n = dim(X)[1] / p
  d = dim(X)[2] / p;
  
  alpha_iter = 0
  alphap_iter = 0
  beta_time = 0
  alpha_time = 0
  alphap_time = 0
  
  

  fam1 = fam()
  
  #Set up initial values
  beta0 = init_beta
  phi0 = init_phi
  alpha0 = init_alpha
  rho0 = alpha0 / alpha0[1]
  
  
  #Beta Step
  iter3 = 1
  start.time = Sys.time()
  f = GEE(beta0,phi0,alpha0,W,Y,X,fam)
  beta1 = solve(a=f$jacobian, b = -stepsize *  f$fval) + beta0
  end.time = Sys.time()
  beta_time = beta_time + (end.time - start.time)
  
  #Phi Step
  ETA = X %*% beta1
  MU = fam1$linkinv(as.vector(ETA))
  YMINUSMU = Y-MU
  YMINUSMU = matrix(YMINUSMU, nrow = p, ncol = n, byrow = FALSE)
  if(fam1$family == 'neg_binomial'){
    phi1 = NULL
    Ymat = matrix(Y,nrow = p, ncol = n, byrow = FALSE)
    mumat = matrix(MU,nrow = p, ncol = n, byrow = FALSE)
    for(species_j in 1:p){
      phi1 = c(phi1, 1/theta.mm(Ymat[species_j,],mumat[species_j,],n-d))
    }
  }else{
    phi1 = apply(((YMINUSMU)^2)/matrix(fam1$variance(MU), nrow = p , ncol = n,byrow = FALSE),1,sum)
    phi1 = phi1/(n-d)
  }
  
  #Alpha step
  if(fam1$family == 'neg_binomial'){
    hMU_half = (matrix(fam1$variance(MU,rep(phi1,n)), nrow = p , ncol = n,byrow = FALSE))^(-1/2)
  }else{
    hMU_half = (phi1*matrix(fam1$variance(MU), nrow = p , ncol = n,byrow = FALSE))^(-1/2)
  }
  
  RESIDUAL = (YMINUSMU) * hMU_half
  
  
  #OLS
  start.time = Sys.time()
  result_OLS = cov_OLS_func(W,RESIDUAL)
  OLS_AA = result_OLS$AA ; hatalpha_OLS = result_OLS$hatalpha 
  end.time = Sys.time()
  alpha_time = alpha_time + (end.time - start.time)
  
  alpha_iter = alpha_iter + 1
  
  
  
  #OLS+
  start.time = Sys.time()
  result_OLSp = cov_OLSp_func(W,RESIDUAL,OLS_AA,hatalpha_OLS,Lambda0,mu,epsilon,xi,iterMAX); 
  hatalpha_OLSp = result_OLSp$hatalpha ; hatSigma_OLSp = result_OLSp$hatSigma ; iter_OLSp = result_OLSp$iter
  end.time = Sys.time()
  alphap_time = alphap_time + (end.time - start.time)
  
  alphap_iter = alphap_iter + iter_OLSp
  
  
  alpha1 = hatalpha_OLSp
  Sigma1 = hatSigma_OLSp
  
  rho1 = alpha1 / alpha1[1]
  R1 = Sigma1/ alpha1[1]
  
  
  
  while( ( (norm(beta1 - beta0, type = "2") > xi2 )| (norm(rho1 - rho0, type = "2") > xi2)  )  & iter3 < iterMAX2 ){
    iter3 = iter3 + 1
    beta0 = beta1
    phi0 = phi1
    alpha0 = alpha1
    rho0 = rho1
    
    
    
    #Beta Step
    start.time = Sys.time()
    f = GEE(beta0,phi0,alpha0,W,Y,X,fam)
    beta1 = solve(a=f$jacobian, b = -stepsize *  f$fval) + beta0
    end.time = Sys.time()
    beta_time = beta_time + (end.time - start.time)
    
    
    
    #Phi Step
    ETA = X %*% beta1
    MU = fam1$linkinv(as.vector(ETA))
    YMINUSMU = Y-MU
    YMINUSMU = matrix(YMINUSMU, nrow = p, ncol = n, byrow = FALSE)
    if(fam1$family == 'neg_binomial'){
      phi1 = NULL
      Ymat = matrix(Y,nrow = p, ncol = n, byrow = FALSE)
      mumat = matrix(MU,nrow = p, ncol = n, byrow = FALSE)
      for(species_j in 1:p){
        phi1 = c(phi1, 1/theta.mm(Ymat[species_j,],mumat[species_j,],n-d))
      }
    }else{
      #Modification of ALgorithm 2a
      phi1 = apply(((YMINUSMU)^2)/matrix(fam1$variance(MU), nrow = p , ncol = n,byrow = FALSE),1,sum)
      phi1 = phi1/(n-d)
    }
    
    #Alpha step
    if(fam1$family == 'neg_binomial'){
      hMU_half = (matrix(fam1$variance(MU,rep(phi1,n)), nrow = p , ncol = n,byrow = FALSE))^(-1/2)
    }else{
      hMU_half = (phi1*matrix(fam1$variance(MU), nrow = p , ncol = n,byrow = FALSE))^(-1/2)
    }
    
    
    RESIDUAL = (YMINUSMU) * hMU_half
    
    
    #OLS
    start.time = Sys.time()
    result_OLS = cov_OLS_func(W,RESIDUAL)
    OLS_AA = result_OLS$AA ; hatalpha_OLS = result_OLS$hatalpha 
    end.time = Sys.time()
    alpha_time = alpha_time + (end.time - start.time)
    
    alpha_iter = alpha_iter + 1
    
    
    #OLS+
    start.time = Sys.time()
    result_OLSp = cov_OLSp_func(W,RESIDUAL,OLS_AA,hatalpha_OLS,Lambda0,mu,epsilon,xi,iterMAX); 
    hatalpha_OLSp = result_OLSp$hatalpha ; hatSigma_OLSp = result_OLSp$hatSigma ; iter_OLSp = result_OLSp$iter
    
    end.time = Sys.time()
    alphap_time = alphap_time + (end.time - start.time)
    
    alphap_iter = alphap_iter + iter_OLSp
    
    alpha1 = hatalpha_OLSp
    Sigma1 = hatSigma_OLSp
    
    rho1 = alpha1 / alpha1[1]
    R1 = Sigma1 / alpha1[1]
    
    
  }
  
  hatbeta = beta1
  hatphi = phi1
  hatalpha = alpha1
  hatSigma = Sigma1
  hatrho = rho1
  hatR = R1
  
  return(list("hatbeta" = as.vector(hatbeta), "hatphi" = as.vector(hatphi), "hatalpha" = as.vector(hatalpha), "hatSigma" = hatSigma, "hatrho" = as.vector(hatrho), "hatR" = hatR, "iter3" = iter3,'alpha_iter' = alpha_iter, 'beta_time' = beta_time, 'alpha_time' = alpha_time,'alphap_iter'=alphap_iter, 'alphap_time'= alphap_time))
  
}



#To compute the estimating function for beta and its Jacobian matrix
GEE = function(beta,phi,alpha,W,Y,X,fam){
  p = dim(W)[1]; K = dim(W)[3] 
  n = dim(X)[1] / p
  d = dim(X)[2] / p;
  
  fam = fam()
  
  ETA = X %*% beta
  M = fam$linkinv(as.vector(ETA))
  
  
  if(fam$family == 'neg_binomial'){
    A_half = sparseMatrix(i = 1:(n*p), j = 1:(n*p), x = (fam$variance(M,rep(phi,n)))^(-1/2))
  }else{
    A_half = sparseMatrix(i = 1:(n*p), j = 1:(n*p), x = (rep(phi,n)*fam$variance(M))^(-1/2))
  }
  
  
  
  Sigma_inv = solve(linear_comb(alpha,W))
  Sigma_inv = bdiag2(replicate(n,Sigma_inv,simplify = FALSE)) 
  
  aa = fam$mu.eta(as.vector(ETA)) * X
  dd = crossprod(aa,A_half %*% Sigma_inv %*% A_half) 
  fval = dd %*% (Y-M)
  jacobian = -dd %*% aa
  
  return(list('fval'=fval,'jacobian'=jacobian))
}

#Creates a block diagonal matrices from the list Ms of SYMMETRIC matrices
bdiag2 <- function(Ms){
  l <- length(Ms)
  N <- nrow(Ms[[1]])  
  i0 <- rep(1:N, times=N:1)
  s <- rep(seq(0,(l-1)*N,by=N),each=length(i0))
  i <- rep(i0,l) + s
  j0 <- 1:N -> j
  for(k in 1:N){
    j0 <- j0[-1]
    j <- c(j, j0)
  }
  j <- rep(j,l) + s
  idx <- t(upper.tri(Ms[[1]], diag = TRUE))
  x <- unlist(lapply(Ms, "[", idx))
  sparseMatrix(i, j, x = x, symmetric = TRUE)
}

#To compute linear combination of W_k matrices
linear_comb = function(beta,W){
  p = dim(W)[1]; K = dim(W)[3]
  TW = array(data=0, dim = c(p,p))
  for (i in 1:K){
    TW = TW + beta[i] * W[,,i]
  }
  return(TW)
}

#Wrapper function to compute the Least square estimator of alpha, i.e. hatalpha_{LS}
cov_OLS_func = function(W,Y){
  cov_OLS = OLS_est(W = W, Y = Y)
  hatalpha = cov_OLS$hatalpha
  hatSigma = cov_OLS$hatSigma
  AA = cov_OLS$AA
  
  
  return(list('AA'= AA, 'hatSigma' = hatSigma, 'hatalpha' = hatalpha))
}

#To compute the Least square estimator of alpha, i.e. hatalpha_{LS}
OLS_est = function(W,Y){
  p = dim(W)[1]; n = length(Y)/p
  K = dim(W)[3] 
  
  aa = array(0, dim = c(K,K))
  for (k1 in 1:K){
    for (k2 in 1:k1){
      aa[k1,k2] = matrix.trace(W[,,k1] %*% W[,,k2])
    }
  }
  aa[upper.tri(aa)] =  t(aa)[upper.tri(aa)] 
  AA = n * aa
  
  
  BB = NULL
  for(k in 1:K){
    BB = c(BB, matrix.trace(crossprod(Y,W[,,k] %*% Y)) )
  }
  
  
  
  hatalpha = solve(AA)%*% BB  
  
  
  hatSigma = linear_comb(hatalpha,W)
  
  
  
  return(list("AA" = AA, "hatalpha" = hatalpha, "hatSigma" = hatSigma))
  
}

#To compute the constrained estimator of alpha, i.e. hatalpha
cov_OLSp_func = function(W,Y,AA,alpha0,Lambda0,mu,epsilon,xi,iterMAX){
  cov_OLSp = OLSp_est(W,Y,AA,alpha0,Lambda0,mu,epsilon,xi,iterMAX)
  hatalpha = cov_OLSp$hatalpha
  hatSigma = cov_OLSp$hatSigma
  iter = cov_OLSp$iter
  
  
  return(list('hatSigma' = hatSigma, 'hatalpha' = hatalpha ,'iter' = iter))
  
}

OLSp_est = function(W,Y,AA,alpha0,Lambda0,mu,epsilon,xi,iterMAX){
  p = dim(W)[1]; n = length(Y) / p
  K = dim(W)[3] 
  
  hatalphaOLS=alpha0;
  
  hatSigmaOLS = linear_comb(hatalphaOLS,W)
  
  test = (is.positive.definite(hatSigmaOLS) & rankMatrix(hatSigmaOLS) == p) 
  
  iter = 0
  if (test == 1){
    hatalpha = hatalphaOLS
    hatSigma = hatSigmaOLS
  }else{
    hatalphaOLS = (2*mu/ (2*mu + 1)) * hatalphaOLS
    AA = solve(AA) * (1/ (2*mu + 1)) 
    
    iter = iter + 1
    
    #Theta Step
    Sigma0 = hatSigmaOLS
    Theta1 = Sigma0 + mu*Lambda0
    
    E = eigen(Theta1)
    D = E$values; Gamma = E$vectors
    D[D <epsilon] = epsilon
    Theta1= quad.tform(diag(D),Gamma) 
    if (!is.symmetric.matrix(Theta1)){
      Theta1 = (Theta1 + t(Theta1))/2
    }
    
    #Alpha Step
    ThetaL1 = Theta1 - mu * Lambda0
    BB = c()
    for (k in 1:K){
      BB = c(BB,n * matrix.trace( W[,,k] %*% ThetaL1)) 
    }
    
    
    alpha1 = hatalphaOLS + AA %*% BB 
    
    #Lambda Step
    Sigma1 = linear_comb(alpha1,W)
    Lambda1 = Lambda0 - ((Theta1 - Sigma1) / (mu) )
    
    while( (norm(alpha1 - alpha0, type = "2") > xi) & (iter < iterMAX) ){
      iter = iter + 1
      
      alpha0 = alpha1
      Lambda0 = Lambda1
      Sigma0 = Sigma1
      
      #Theta Step
      Theta1 = Sigma0 + mu*Lambda0
      
      E = eigen(Theta1)
      D = E$values; Gamma = E$vectors
      D[D <epsilon] = epsilon
      Theta1= quad.tform(diag(D),Gamma) 
      if (!is.symmetric.matrix(Theta1)){
        Theta1 = (Theta1 + t(Theta1))/2
      }
      
      #Alpha Step
      ThetaL1 = Theta1 - mu * Lambda0
      BB = c()
      for (k in 1:K){
        BB = c(BB,n * matrix.trace( W[,,k] %*% ThetaL1)) 
      }
      
      alpha1 = hatalphaOLS + AA %*% BB 
      
      #Lambda Step
      Sigma1 = linear_comb(alpha1,W)
      Lambda1 = Lambda0 - ((Theta1 - Sigma1) / (mu) )
      
    }
    hatalpha = alpha1;
    hatSigma = Sigma1
  }
  
  
  
  
  return(list("hatalpha"= hatalpha, "hatSigma" = hatSigma, "iter" = iter))
}

#Family object for negative binomial distribution
neg_binomial <- function(link = "log") {
  out <- poisson(link)
  out$family <- "neg_binomial"
  out$variance <- function(mu, phi =0) mu + ((mu^2) * phi)
  out$dev.resids <- function(y, mu, wt) {
    stop("'dev.resids' function should not be called")
  }
  out$aic <- function(y, n, mu, wt, dev) {
    stop("'aic' function should not have been called")
  }
  out$simulate <- function(object, nsim)
    stop("'simulate' function should not have been called")
  return(out)
}

#Reshape the original X matrix (of dimension d * n)  into an np x pd matrix for use of the estimation algorithm, which can be used to save computational time in the simulation study
reshape_X = function(W,X){
  d = dim(X)[1];
  p = dim(W)[1]; K = dim(W)[3]
  
  X = apply(X,2,function(x){list( bdiag(replicate(p,t(x),simplify = FALSE)))})
  X = do.call('rbind',unlist(X,recursive=FALSE))
  
  return(list( 'X' = X, 'W' = W))
}

#Compute the matrix required for the multivariate Delta method
delta_method = function(beta,phi,alpha,rho,W,Y,X,fam){
  p = dim(W)[1]; K = dim(W)[3]
  d = dim(X)[2] / p;
  
  r1 = array(data = 0, dim = c(K-1,K-1))
  diag(r1) = 1/alpha[1]
  
  r2 = -alpha[2:K]/ ((alpha[1])^2)
  
  r3 = cbind(r2,r1)
  
  result1 = array(data = 0 , dim = c(p*d,p*d))
  diag(result1) = 1
  
  result = bdiag(list(result1,r3))
  
  return(as.matrix(result))
}

#Compute the asymptotic covariance matrix for the joint estimator of mean and reparameterized correlation parameter
sandwich_cov_matrix = function(beta,phi,alpha,rho,W,Y,X,fam){
  cov_psi = as.matrix(cov_matrix(beta,phi,alpha,rho,W,Y,X,fam))
  derivative = as.matrix(first_derivative(beta,phi,alpha,rho,W,Y,X,fam))
  
  sandwich = quad.tform(cov_psi,solve(derivative))
  return(sandwich)
}

#Compute the covariance matrix of the estimating equation
cov_matrix = function(beta,phi,alpha,rho,W,Y,X,fam){
  p = dim(W)[1]; K = dim(W)[3] 
  n = dim(X)[1] / p
  d = dim(X)[2] / p;
  
  fam = fam()
  
  ETA = X %*% beta
  M = fam$linkinv(as.vector(ETA))
  if(fam$family == 'neg_binomial'){
    A_half = sparseMatrix(i = 1:(n*p), j = 1:(n*p), x = (fam$variance(M,rep(phi,n)))^(-1/2))
  }else{
    A_half = sparseMatrix(i = 1:(n*p), j = 1:(n*p), x = (rep(phi,n)*fam$variance(M))^(-1/2))
  }
  
  Sigma_inv = solve(linear_comb(alpha,W))
  Sigma_inv = bdiag2(replicate(n,Sigma_inv,simplify = FALSE)) 
  
  aa = fam$mu.eta(as.vector(ETA)) * X 
  dd = crossprod(aa,A_half %*% Sigma_inv %*% A_half)
  V_beta = dd %*% aa
  
  
  
  ###V_alpha  
  Sigma = linear_comb(alpha,W)
  L_Sigma = t(chol(Sigma))
  L_Sigma_Inv = solve(L_Sigma)
  Sigma= bdiag2(replicate(n,Sigma,simplify = FALSE))
  L_Sigma_block = bdiag(replicate(n,L_Sigma,simplify = FALSE))
  L_Sigma_Inv_block = bdiag(replicate(n,L_Sigma_Inv,simplify = FALSE))
  
  aa1 = array(0,dim = c(K,K))
  for (k1 in 1:K){
    for (k2 in 1:k1){
      aa1[k1,k2] = matrix.trace(t(L_Sigma) %*%  W[,,k1] %*% L_Sigma %*% t(L_Sigma)  %*% W[,,k2] %*% L_Sigma )
    }
  }
  aa1[upper.tri(aa1)] =  t(aa1)[upper.tri(aa1)]
  AA1 = n * aa1
  
  ###Estimation for mu4
  vareps = L_Sigma_Inv_block %*% A_half %*% (Y-M)
  mu3 = mean(vareps^3)
  mu4 = mean(vareps^4)
  
  
  
  #computing \mathcal{B} matrix
  
  B_pre = array(data = 0,dim = c(p, K))
  for(k1 in 1:K){
    B_pre[,k1] = diag(crossprod(L_Sigma,W[,,k1] %*% L_Sigma))
  }
  B = do.call('rbind',replicate(n,B_pre,simplify = FALSE))
  
  V_alpha = 2* AA1 + (mu4 -3) * crossprod(B,B)
  
  
  ##V_beta_alpha
  cc = crossprod(aa,A_half)
  ee = crossprod(L_Sigma_Inv_block, B)
  V_beta_alpha = mu3 * cc %*% ee
  
  result1 = cbind(V_beta, V_beta_alpha)
  result2 = cbind(t(V_beta_alpha),V_alpha)
  result = rbind(result1,result2)
  return(result)
}

#Compute the expectation of the first derivative of the estimating equations
first_derivative = function(beta,phi,alpha,rho,W,Y,X,fam){
  p = dim(W)[1]; K = dim(W)[3] 
  n = dim(X)[1] / p
  d = dim(X)[2] / p;
  
  fam = fam()
  
  ETA = X %*% beta
  M = fam$linkinv(as.vector(ETA))
  if(fam$family == 'neg_binomial'){
    A_half = sparseMatrix(i = 1:(n*p), j = 1:(n*p), x = (fam$variance(M,rep(phi,n)))^(-1/2))
  }else{
    A_half = sparseMatrix(i = 1:(n*p), j = 1:(n*p), x = (rep(phi,n)*fam$variance(M))^(-1/2))
  }
  
  Sigma_inv = solve(linear_comb(alpha,W))
  Sigma_inv = bdiag2(replicate(n,Sigma_inv,simplify = FALSE)) 
  
  aa = fam$mu.eta(as.vector(ETA)) * X
  dd = crossprod(aa,A_half %*% Sigma_inv %*% A_half)
  S_beta = -dd %*% aa
  
  aa = array(0, dim = c(K,K)) 
  for (k1 in 1:K){
    for (k2 in 1:k1){
      aa[k1,k2] = matrix.trace(W[,,k1] %*% W[,,k2])
    }
  }
  aa[upper.tri(aa)] =  t(aa)[upper.tri(aa)] #here aa will be converted to dsyMatrix class
  S_alpha = -n * aa
  
  S_beta_alpha = matrix(data=0, nrow = p*d, ncol = K)
  
  
  Sigma = linear_comb(alpha,W)
  Sigma_block = bdiag2(replicate(n,Sigma,simplify = FALSE))
  

  
  if(fam$family == 'gaussian'){
    h_prime = function(x){rep(0,length(x))}
  }else if(fam$family == 'poisson'){
    h_prime = function(x){rep(1,length(x))}
  }else if(fam$family == 'binomial'){
    h_prime = function(x){1 - 2*x}
  }else if(fam$family == 'neg_binomial'){
    h_prime = function(x,phi){1+2*x*phi}
  }
  
  
  temp = array(0,dim=c(K,p))
  for(k in 1:K){
    temp[k,] = diag(Sigma %*% W[,,k])
  }
  
  if(fam$family == 'neg_binomial'){
    temp2 = h_prime(M,rep(phi,n)) / fam$var(M,rep(phi,n))
  }else{
    temp2 = h_prime(M) / fam$var(M)
  }
  
  D = fam$mu.eta(as.vector(ETA)) * X
  S_alpha_beta = 0
  for(i in 1:n){
    Di = D[((i-1)*p+1):(i*p),]
    S_alpha_beta = S_alpha_beta - temp %*% diag(temp2[((i-1)*p+1):(i*p)]) %*% Di
  }
  
  
  S1 = cbind(S_beta, S_beta_alpha)
  S2 = cbind(S_alpha_beta,S_alpha)
  S = rbind(S1,S2)
  return(S)
  
}

#Compute the norm of the estimating equation vector
norm_OLSest_eq = function(beta,phi,alpha,W,Y,X,fam){
  p = dim(W)[1]; K = dim(W)[3] 
  n = dim(X)[1] / p
  d = dim(X)[2] / p;
  
  fam1 = fam()
  
  ETA = X %*% beta
  M = fam1$linkinv(as.vector(ETA))
  if(fam1$family == 'neg_binomial'){
    A_half = sparseMatrix(i = 1:(n*p), j = 1:(n*p), x = (fam1$variance(M,rep(phi,n)))^(-1/2))
  }else{
    A_half = sparseMatrix(i = 1:(n*p), j = 1:(n*p), x = (rep(phi,n)*fam1$variance(M))^(-1/2))
  }
  
  Sigma_inv = solve(linear_comb(alpha,W))
  Sigma_inv = bdiag2(replicate(n,Sigma_inv,simplify = FALSE)) 
  
  aa = fam1$mu.eta(as.vector(ETA)) * X
  dd = crossprod(aa,A_half %*% Sigma_inv %*% A_half)
  fval = dd %*% (Y-M)
  
  
  #Get the residual
  ETA = X %*% beta
  MU = fam1$linkinv(as.vector(ETA))
  YMINUSMU = Y-MU
  if(fam1$family == 'neg_binomial'){
    hMU_half = (matrix(fam1$variance(MU,rep(phi,n)), nrow = p , ncol = n,byrow = FALSE))^(-1/2)
  }else{
    hMU_half = (phi*matrix(fam1$variance(MU), nrow = p , ncol = n,byrow = FALSE))^(-1/2)
  }
  
  
  RESIDUAL = (YMINUSMU) * hMU_half
  
  
  
  aa = array(0, dim = c(K,K)) 
  for (k1 in 1:K){
    for (k2 in 1:k1){
      aa[k1,k2] = matrix.trace(W[,,k1] %*% W[,,k2])
    }
  }
  aa[upper.tri(aa)] =  t(aa)[upper.tri(aa)] 
  AA = n * aa
  
  
  BB = c()
  for(k in 1:K){
    BB = c(BB, matrix.trace(crossprod(RESIDUAL,W[,,k] %*% RESIDUAL)) )
  }
  
  
  est_alpha = BB - AA %*% alpha
  
  result = rbind(fval,est_alpha)
  result_norm = norm(result,type = '2')
  
  return(result_norm)
}