rm(list=ls())

############################################################################################
#Load Data and Pre-preprocessing
############################################################################################

library(ade4)
library(here)
i_am("Applications/aravo.R")
data(aravo)

#Center and scale each covariate to have mean zero and variance one
X_matrix = cbind(1,data.frame(scale(aravo$env[,c('Aspect','Slope','PhysD')])))
X_matrix = t(as.matrix(sapply(X_matrix,as.numeric)))



#Only considered ground beetle species that were detected in more than 10 sites
include.ind = apply(aravo$spe,2,function(x){75 - sum(x == 0)}) > 10
p = sum(include.ind)

Y = t(as.matrix(aravo$spe))
binary_Y = matrix(as.numeric(Y > 0),nrow = dim(Y)[1], byrow = FALSE)
Y_matrix = as.matrix(binary_Y)[include.ind,]

#Compute initial values 
init_beta = NULL
for(j in 1:p){
  glm_Y = Y_matrix[j,]
  glm_X = t(X_matrix[-1,])
  fit = glm(glm_Y ~ glm_X, family = binomial)
  init_beta = c(init_beta,fit$coefficients)
}

#Center and scale quantitative trait to have mean zero and variance one
standardized_traits = scale(aravo$traits[include.ind,],center = T,scale = T)

# Construct similarity matrix for quantitative traits
K = dim(aravo$traits)[2]
W = array(data =  NA,dim = c(p,p,K + 1))

W[,,1] = diag(p)

for(k in 2:(K+1)){
  W_pre = array(data = NA, dim = c(p,p))
  for(j1 in 1:p){
    for(j2 in 1:j1){
      W_pre[j1,j2] = exp(-(standardized_traits[j1,(k-1)] - standardized_traits[j2,(k-1)])^2)
    }
  }
  W_pre[upper.tri(W_pre)] =  t(W_pre)[upper.tri(W_pre)]
  diag(W_pre) = 0
  
  W[,,k] = W_pre
}

############################################################################################
#Fit Joint Mean and Correlation Regression Model
############################################################################################
source(here("Codes","Functions.R"))
X = reshape_X(W,X_matrix)$X
result = estimate_joint(X=X,Y = c(Y_matrix),W=W,fam = binomial,init_beta = init_beta, init_phi = rep(1,p), init_alpha = c(1,rep(0,K)),stepsize=1,Lambda0=array(0,dim = c(p,p)),mu=0.05,epsilon = 1e-5,xi = 1e-8,iterMAX = 1000, xi2 = 1e-8, iterMAX2=500)
est_delta_matrix = delta_method(result$hatbeta,result$hatphi,result$hatalpha,result$hatrho,W,c(Y_matrix),X,binomial)
est_final_cov_matrix = quad.tform(sandwich_cov_matrix(result$hatbeta,result$hatphi,result$hatalpha,result$hatrho,W,c(Y_matrix),X,binomial),est_delta_matrix)
lower_ci = c(result$hatbeta,result$hatrho[2:(K+1)]) - qnorm(0.975) * sqrt(diag(est_final_cov_matrix))
upper_ci = c(result$hatbeta,result$hatrho[2:(K+1)]) + qnorm(0.975) * sqrt(diag(est_final_cov_matrix))

############################################################################################
#Estimation Results
############################################################################################

#Formatting species names
sub_aravo_names = aravo$spe.names[include.ind]
sub_col_names = colnames(aravo$spe)[include.ind]
species_names = rep(0,p)
for(v in 1:p){
  p3 = 0
  p1 = substr(sub_aravo_names[v],start = 1,stop=1)
  p2 = '.'
  for(x in unlist(strsplit(sub_aravo_names[v]," "))){
    if(grepl(strsplit(sub_col_names[v],"[.]")[[1]][2],x)){
      p3 = x
    }
  }
  species_names[v] = paste0(p1,p2,p3)
}
species_names[31] = 'P.mutellinoides'
beta_result = data.frame(species = rep(species_names, each = 4), covariate = rep(c('Intercept', 'Aspect', 'Slope', 'Physical Disturbance'), p),hatbeta = result$hatbeta, beta_lower_CI = lower_ci[1:(4*p)], beta_upper_CI = upper_ci[1:(4*p)])
beta_result

rho_result = data.frame(trait = c('Plant Height', 'MLS', 'LEA', 'Leaf Area', 'Leaf Thickness','SLA','Nitrogen Content', 'Seed Mass'),
                        hatrho = result$hatrho[2:(K+1)],
                        rho_lower_CI = lower_ci[(4*p + 1):(4*p + K)], rho_upper_CI = upper_ci[(4*p + 1):(4*p + K)])
rho_result
