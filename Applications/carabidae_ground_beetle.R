############################################################################################
#Load Data and Pre-preprocessing
############################################################################################

rm(list=ls())
library(MASS)
library(ggplot2)
library(here)
i_am("Applications/carabidae_ground_beetle.R")
beetle_env_file = read.csv(here("Applications","beetle_env.csv"))
beetle_spe_file = read.csv(here("Applications","beetle_spe.csv"))
beetle_trait_file = read.csv(here("Applications","beetle_trait.csv"))

#Center and scale each covariate to have mean zero and variance one
X_matrix = cbind(1,data.frame(scale(beetle_env_file[,c('pH','Elevation','Management')])))
X_matrix = t(as.matrix(sapply(X_matrix,as.numeric)))

beetle_spe = data.frame(beetle_spe_file[,2:88],row.names = beetle_trait_file$SPECIES)


#Only considered ground beetle species that were detected in at least 15 sites
include.ind = apply(beetle_spe,1,function(x){sum(x>0)}) >= 15

#Also, removed 'pter nige' and 'pter rhae' which caused errors in computing the initial values of their beta coefficients
include.ind[c(58,60)] = F
p = sum(include.ind)

Y_matrix = as.matrix(beetle_spe)[include.ind,]

spe_mean = apply(Y_matrix ,1,mean)
spe_var = apply(Y_matrix ,1,var)

spe_mean/spe_var


log_log_plot = ggplot(data = data.frame(spe_mean,spe_var), aes(x = spe_mean, y = spe_var)) + geom_point() + scale_x_log10() + scale_y_log10() + theme_bw() + 
  xlab('Species Mean') + ylab('Species Variance') + ggtitle('Log-log Plot for Sample Variances Against Sample Means of Species Counts') + theme(plot.title = element_text(hjust = 0.5))
log_log_plot

#Compute initial values 
init_beta = NULL
init_phi = NULL
for(j in 1:p){
  print(j)
  glm_Y = Y_matrix[j,]
  glm_X = t(X_matrix[-1,])
  fit = glm.nb(glm_Y ~ glm_X,control = glm.control(maxit = 100))
  init_beta = c(init_beta,fit$coefficients)
}

#Center and scale quantitative trait to have mean zero and variance one
beetle_trait = data.frame(beetle_trait_file[,4:8], row.names = beetle_trait_file$CODE)

quan_trait = beetle_trait[include.ind,'LTL']
quan_trait = scale(quan_trait,center = T,scale = T)

qual_trait = beetle_trait[include.ind,c("CLG","WIN","OVE","BRE")]
qual_trait = sapply(qual_trait,function(x){gsub("\\?","",x)})


K = 5

W = array(data =  NA,dim = c(p,p,K+1))
W[,,1] = diag(p)

# Construct similarity matrix for quantitative trait
W_pre = array(data = NA, dim = c(p,p))
for(j1 in 1:p){
  for(j2 in 1:j1){
    W_pre[j1,j2] = exp(-(quan_trait[j1] - quan_trait[j2])^2)
  }
}
W_pre[upper.tri(W_pre)] =  t(W_pre)[upper.tri(W_pre)]
diag(W_pre) = 0
W[,,2] = W_pre

# Construct similarity matrix for qualitative traits
for(k in 3:6){
  W_pre = array(data = NA, dim = c(p,p))
  for(j1 in 1:p){
    for(j2 in 1:j1){
      if(qual_trait[j1,(k - 2)] == qual_trait[j2,(k - 2)]){
        W_pre[j1,j2] = 1
      }else{
        W_pre[j1,j2] = 0
      }
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
result = estimate_joint(X=X,Y = c(Y_matrix),W=W,fam = neg_binomial,init_beta = init_beta, init_phi = rep(1,p), init_alpha = c(1,rep(0,K)),stepsize=0.5,Lambda0=array(0,dim = c(p,p)),mu=0.05,epsilon = 1e-5,xi = 1e-8,iterMAX = 1000, xi2 = 1e-8, iterMAX2=1500)
est_delta_matrix = delta_method(result$hatbeta,result$hatphi,result$hatalpha,result$hatrho,W,c(Y_matrix),X,neg_binomial)
est_final_cov_matrix = quad.tform(sandwich_cov_matrix(result$hatbeta,result$hatphi,result$hatalpha,result$hatrho,W,c(Y_matrix),X,neg_binomial),est_delta_matrix)
lower_ci = c(result$hatbeta,result$hatrho[2:(K+1)]) - qnorm(0.975) * sqrt(diag(est_final_cov_matrix))
upper_ci = c(result$hatbeta,result$hatrho[2:(K+1)]) + qnorm(0.975) * sqrt(diag(est_final_cov_matrix))

############################################################################################
#Estimation Results
############################################################################################
species_names = paste0(sapply(rownames(beetle_spe)[include.ind],substr,start = 1,stop=1),'.',sapply(strsplit(rownames(beetle_spe)[include.ind], ' '), `[`, 2))
beta_result = data.frame(species = rep(species_names, each = 4), covariate = rep(c('Intercept', 'Soil pH', 'Elevation', 'Land management'), p),hatbeta = result$hatbeta, beta_lower_CI = lower_ci[1:(4*p)], beta_upper_CI = upper_ci[1:(4*p)])
beta_result

rho_result = data.frame(trait = c('Total length', 'Leg color', 'Wing development', 'Overwintering', 'Breeding season'),
                        hatrho = result$hatrho[2:(K+1)],
                        rho_lower_CI = lower_ci[(4*p + 1):(4*p + K)], rho_upper_CI = upper_ci[(4*p + 1):(4*p + K)])
rho_result


