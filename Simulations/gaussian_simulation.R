############################################################################################
#Simulate Data
############################################################################################

rm(list=ls())
library(here)
i_am("Simulations/gaussian_simulation.R")
source(here("Codes","Functions.R"))
X_in = readRDS(here("Simulations","X_matrix.rds"))
W_in = readRDS(here("Simulations","W_matrix.rds"))

#Randomly generate true rho and compute the true correlation matrix
n = 400; p =50; d= 4; K = 5

set.seed(255)
rho = c(1,runif(K,-0.05,0.05))
R0 = linear_comb(rho,W_in)


#Randomly generate true beta and compute the mean of each response Y_{ij}
set.seed(7)
beta50 = rnorm(p*d,mean = 0, sd = 0.5)
BETA = matrix(beta50,nrow = p, byrow = T)
ETA = BETA %*% X_in
MU = gaussian()$linkinv(ETA)



p_vec = c(10,25,50)
n_vec = c(50,100,200,400)
count = 1

#Create a list consisting of the W_matrices for different combinations of n = 50/100/200/400 and p=10/25/50
W_list = list()
for (i in 1:3){
  for (j in 1:4){
    temp  = array(data=NA, dim = c(p_vec[i],p_vec[i],K+1))
    
    for (t in 1:(K+1)){
      temp[,,t] = W_in[1:p_vec[i],1:p_vec[i],t]
    }
    W_list[[count]] = temp
    count = count +1
  }
}


#Set true phi to be ones
phi50 = rep(1,50)
A_half = (phi50*matrix(gaussian()$variance(MU), nrow = p , ncol = n))^(1/2)

#Compute true covariance matrix
true_cov_matrix = array(data = NA, dim=c(p,p,n))
for (i in 1:n){
  true_cov_matrix[,,i] = diag(A_half[,i]) %*% R0 %*% diag(A_half[,i])
  if(!is.symmetric.matrix(true_cov_matrix[,,i])){
    true_cov_matrix[,,i] = (true_cov_matrix[,,i] + t(true_cov_matrix[,,i]))/2
  }
}

B=1000
data_Y = array(data = NA, dim= c(p,B,n))

set.seed(1)
for(b in 1:B){
  for (i in 1:n){
    data_Y[,b,i] = mvrnorm(n=1, mu = MU[,i], Sigma = true_cov_matrix[,,i])
  }
}

saveRDS(data_Y,here("Simulations","gaussian_data_Y.rds"))



############################################################################################
#Run Simulation
############################################################################################
B_converge = 1000
count = 1
for(p in p_vec){
  for(n in n_vec){
    #Family
    fam = gaussian
    
    #File path to pre-generated data
    reuse_path = here("Simulations","gaussian_data_Y.rds")
    
    #File path to save simulation result
    save_path = here('Simulations',paste0('gaussian_n',n,'_p',p,'_sim_result.rds'))
    
    #W_k matrices and the vectorized form of X matrix
    W = W_list[[count]]; X = reshape_X(W_list[[count]],X_in[,1:n])$X
    
    #Tuning parameters
    stepsize=1;Lambda0=array(0,dim = c(p,p));mu=0.05;epsilon = 1e-5;xi = 1e-8;iterMAX = 1000; xi2 = 1e-8; iterMAX2=1000; 
    
    #Chi-sq coverage for the first q response's coefficients
    q = 5
    
    K = dim(W)[3]
    d = dim(X)[2] / p;
    
    #True parameters
    beta = beta50[1:(p*d)]; phi = phi50[1:p]; rho = rho
    
    
    #Set up initial values for the estimation algorithm
    init_beta = beta
    init_phi = rep(1,p)
    init_alpha = c(1,rep(0,K-1))
    
    
    fam2 = fam()
    
    R0 = linear_comb(rho,W)
    
    
    reuse_Y = readRDS(reuse_path)
    data_Y = reuse_Y[1:p,,1:n]
    
    
    ETA = X %*% beta
    MU = fam2$linkinv(as.vector(ETA))
    MU = matrix(MU, nrow = p, ncol = n, byrow = FALSE)
    
    A_half = (phi*matrix(fam2$variance(MU), nrow = p , ncol = n))^(1/2)
    
    covMatrix = array(data = NA, dim=c(p,p,n))
    for (i in 1:n){
      covMatrix[,,i] = diag(A_half[,i]) %*% R0 %*% diag(A_half[,i])
      if(!is.symmetric.matrix(covMatrix[,,i])){
        covMatrix[,,i] = (covMatrix[,,i] + t(covMatrix[,,i]))/2
      }
    }
    
    
    #Set up list of dimnames for the arrays that will be used to store simulation results
    namelist = list(c(rep(sprintf("beta_[%d]",seq(1:p)),each = d), sprintf("phi_[%d]",seq(1:p)), sprintf("rho_[%d]",seq(1:K)-1), 'S-norm' ,'F-norm', 'Iter', 'Alpha_Iter','Alpha_time','Alphap_iter','Alphap_time','Beta_time') , c(1:B))
    alpha_namelist = list(c(sprintf("alpha_[%d]",seq(1:K)-1),'S-norm_Sigma' ,'F-norm_Sigma'), c(1:B))
    namelist_CI = list(c(rep(sprintf("beta_[%d]",seq(1:p)),each = d), sprintf("rho_[%d]",seq(1:(K-1)))),c('Lower_Est','Upper_Est'),c(1:B))
    namelist_chisq = list(1:B,c('Beta_Est','Rho_Est'))
    
    #Set up arrays to store simulation results later
    l = p*d + p + K 
    result = array(NA, dim = c((l+8),B), dimnames = namelist)
    alpha_result = array(NA,dim = c(K+2,B), dimnames = alpha_namelist)
    
    OLSp_CIs = array(NA, dim = c((l-p-1),2,B), dimnames = namelist_CI )
    chisq_stat = array(NA, dim = c(B,2), dimnames = namelist_chisq)
    OLSest_norm = rep(NA,B)
    
    time.OLSp = 0
    b=0
    b_converge=0
    while(b_converge<B_converge){
      b= b+1 
      if(b %% 100 == 0 ){
        cat(paste0('p = ',p,'; n = ',n,'; b = ', b, '\n'))
      }
      tryCatch({
        Y = array(data = NA, dim = c(p,n))
        for(i in 1:n){
          Y[,i] = data_Y[,b,i]
        }
        Y = c(Y)
        
        
        
        start.time.OLSp = Sys.time()
        result_OLSp = estimate_joint(X=X,Y = Y,W=W,fam = fam,init_beta = init_beta, init_phi = init_phi, init_alpha = init_alpha,stepsize=stepsize,Lambda0=Lambda0,mu=mu,epsilon = epsilon,xi = xi,iterMAX = iterMAX, xi2 = xi2, iterMAX2=iterMAX2)
        end.time.OLSp = Sys.time()
        time.OLSp = time.OLSp +  (end.time.OLSp - start.time.OLSp)
        
        
        spectral_norm_hatR = norm(result_OLSp$hatR-R0, type = "2") 
        f_norm_hatR =norm(result_OLSp$hatR-R0, type = "F") / sqrt(p)
        
        spectral_norm_hatSigma = norm(result_OLSp$hatSigma-R0, type = "2")
        f_norm_hatSigma = norm(result_OLSp$hatSigma-R0, type = "F")/sqrt(p)
        
        est_delta_matrix = delta_method(result_OLSp$hatbeta,result_OLSp$hatphi,result_OLSp$hatalpha,result_OLSp$hatrho,W,Y,X,fam)
        est_final_cov_matrix = quad.tform(sandwich_cov_matrix(result_OLSp$hatbeta,result_OLSp$hatphi,result_OLSp$hatalpha,result_OLSp$hatrho,W,Y,X,fam),est_delta_matrix)
        
        
        
        chisq_stat_beta_est = quad.form.inv(est_final_cov_matrix[1:(q*d), 1:(q*d)],result_OLSp$hatbeta[1:(q*d)] - beta[1:(q*d)])
        chisq_stat_rho_est = quad.form.inv(est_final_cov_matrix[(p*d+1):(p*d + K -1), (p*d+1):(p*d + K -1)],result_OLSp$hatrho[2:K] - rho[2:K])
        chisq_stat[b,] = c(chisq_stat_beta_est,chisq_stat_rho_est)
        
        OLSp_CIs[,1,b] = c(result_OLSp$hatbeta,result_OLSp$hatrho[2:K]) - qnorm(0.975) * sqrt(diag(est_final_cov_matrix))
        OLSp_CIs[,2,b] = c(result_OLSp$hatbeta,result_OLSp$hatrho[2:K]) + qnorm(0.975) * sqrt(diag(est_final_cov_matrix))
        
        
        OLSest_norm[b] = norm_OLSest_eq(beta = result_OLSp$hatbeta, phi = result_OLSp$hatphi, alpha = result_OLSp$hatalpha, W = W,Y = Y, X=X,fam = fam)
        
        result[,b] = c(result_OLSp$hatbeta,result_OLSp$hatphi,result_OLSp$hatrho,spectral_norm_hatR,f_norm_hatR, result_OLSp$iter3,result_OLSp$alpha_iter, result_OLSp$alpha_time,result_OLSp$alphap_iter,result_OLSp$alphap_time,result_OLSp$beta_time)
        alpha_result[,b] = c(result_OLSp$hatalpha, spectral_norm_hatSigma,f_norm_hatSigma)
        
        if(!is.na(result_OLSp$iter3)){
          if(result_OLSp$iter3 < iterMAX2){
            b_converge=b_converge+1
          }
        }
      },error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
      
      
    }
    
    #Keep track of replications with either errors or non-convergence of the estimation algorithm
    total_trial = b
    num_error= 0
    num_nonconverge=0
    
    errors = rep(FALSE,B) 
    nonconverge = rep(FALSE,B)
    for(b in 1:B){
      if(is.na(result['Iter',b])){
        errors[b] = TRUE
        if(b <= total_trial){
          num_error = num_error+1
        }
      }else{
        if(result['Iter',b] == iterMAX2){
          nonconverge[b] = TRUE
          if(b <= total_trial){
            num_nonconverge = num_nonconverge+1
          }
        }
      }
    }
    
    
    
    
    full_result = result
    result = full_result[,!(errors | nonconverge)]
    
    alpha_result = alpha_result[,!(errors | nonconverge)]
    
    OLSp_CIs = OLSp_CIs[,,!(errors | nonconverge)]
    chisq_stat = chisq_stat[!(errors | nonconverge),]
    OLSest_norm = OLSest_norm[!(errors | nonconverge)]
    
    
    CI_coverage = rowMeans(c(beta,rho[2:K]) >= OLSp_CIs[,1,] & c(beta,rho[2:K]) <= OLSp_CIs[,2,])
    names(CI_coverage) = c(rep(sprintf("beta[%d]",seq(1:p)),each = d),sprintf("rho[%d]",seq(1:(K-1)) ))
    
    chisq_beta_est_coverage = mean(chisq_stat[,1] <= qchisq(0.95, df= (q*d)))  
    chisq_rho_est_coverage = mean(chisq_stat[,2] <= qchisq(0.95, df= (K-1))) 
    
    
    averaged_result = apply(result,1,FUN = function(x){mean(x,na.rm = T)})
    averaged_alpha_result = apply(alpha_result,1,FUN = function(x){mean(x,na.rm = T)})
    
    averaged_beta_bias = averaged_result[1:(p*d)] - beta
    averaged_phi_bias = averaged_result[(p*d + 1):(p*d + p) ] - phi
    averaged_rho_bias = averaged_result[(p*d +p + 1): (p*d + p +K)] - rho
    
    averaged_alpha_bias = averaged_alpha_result[1:K] - rho
    alpha_SD = apply(alpha_result,1, function(x){sd(x,na.rm=T)})[1:K]
    alpha_RMSE = sqrt(   (averaged_alpha_bias)^2 + alpha_SD^2             )
    averaged_S_norm_Sigma = averaged_alpha_result[K+1]
    averaged_f_norm_Sigma = averaged_alpha_result[K+2]
    
    averaged_S_norm = averaged_result[(p*d + p +K + 1)]
    averaged_f_norm = averaged_result[(p*d + p +K + 2)]
    averaged_iter = averaged_result[(p*d + p +K + 3)]
    averaged_alpha_iter = averaged_result[(p*d + p +K + 4)]
    averaged_alpha_time = averaged_result[(p*d + p +K + 5)]
    averaged_alphap_iter = averaged_result[(p*d + p +K + 6)]
    averaged_alphap_time = averaged_result[(p*d + p +K + 7)]
    averaged_beta_time = averaged_result[(p*d + p +K + 8)]
    
    
    nonconverge_ratio = num_nonconverge/total_trial
    error_ratio = num_error/total_trial
    total_replication = total_trial - num_nonconverge - num_error
    
    ESD = apply(result,1, function(x){sd(x,na.rm=T)})[1:l]
    RMSE = sqrt ( (averaged_result[1:l] - c(beta,phi,rho))^2 + ESD ^2  )
    
    summary_result = c(averaged_beta_bias,averaged_phi_bias,averaged_rho_bias,RMSE,averaged_S_norm, averaged_f_norm, averaged_iter, averaged_alpha_iter, averaged_alpha_time,averaged_alphap_iter,averaged_alphap_time, averaged_beta_time,c(nonconverge_ratio,nonconverge_ratio),c(num_error,num_error) )
    names(summary_result) = c(rep(sprintf("Bias(beta[%d])",seq(1:p)),each = d), sprintf("Bias(phi[%d])",seq(1:p)), sprintf("Bias(rho[%d])",seq(1:K)-1),rep(sprintf("RMSE(beta[%d])",seq(1:p)),each = d), sprintf("RMSE(phi[%d])",seq(1:p)), sprintf("RMSE(rho[%d])",seq(1:K)-1) ,'S-norm' ,'F-norm', 'Iter', 'Alpha_Iter','Alpha_time','Alphap_Iter','Alphap_time','Beta_time', 'Non-converge ratio', 'num_error'  )
    
    summary_alpha_result = c(averaged_alpha_bias,alpha_RMSE,averaged_S_norm_Sigma,averaged_f_norm_Sigma)
    names(summary_alpha_result) = c(sprintf("Bias(alpha[%d])",seq(1:K)-1),sprintf("RMSE(alpha[%d])",seq(1:K)-1) ,'S-norm_Sigma' ,'F-norm_Sigma')
    
    
    errors[(total_trial+1):B] = FALSE
    
    sim_result = list('total_trial'=total_trial, 'n' = n, 'p' = p, 'q' = q,'B'=B,'total_replication'=total_replication,'nonconverge'=which(nonconverge),'errors'=which(errors),'full_result'=full_result,'true_R'= R0,'true_beta'=beta, 'true_phi'=phi, 'true_rho' = rho ,'data_Y'= data_Y, 'OLSest_norm'= OLSest_norm,'raw'=result, 'summary'=summary_result, 'alpha_raw' = alpha_result, 'alpha_summary' = summary_alpha_result , 'CI_raw' = OLSp_CIs, 'CI_coverage' = CI_coverage, 'chisq_stat' = chisq_stat, 'chisq_beta_est_coverage' = chisq_beta_est_coverage, 'chisq_rho_est_coverage' = chisq_rho_est_coverage , 'time_OLSp' = time.OLSp)
    
    saveRDS(sim_result,save_path)
    
    count = count + 1
  }
}



############################################################################################
#Format Simulation Results
############################################################################################
library(ggplot2)
library(ggpubr)
library(tidyr)
library(tidyverse)



gaussian_list = vector('list',12)
count = 1
for(p in c(10,25,50)){
  for(n in c(50,100,200,400)){
    result_file_path = here('Simulations',paste0('gaussian_n',n,'_p',p,'_sim_result.rds'))
    sim_result = readRDS(result_file_path)
    gaussian_list[[count]] = sim_result
    count = count+1
  }
}

#MSE
dat_gaussian = expand.grid(c(50,100,200,400),c(10,25,50))
colnames(dat_gaussian) = c('n','p')
dat_gaussian$p = as.factor(dat_gaussian$p)

dat_gaussian$type = 'gaussian'
dat_gaussian$mse_beta = c(sapply(gaussian_list,FUN = function(x){
  B = dim(x$raw)[2];
  K = length(x$true_rho)
  d = length(x$true_beta) / x$p
  p = x$p
  biases = x$raw[1:(p*d),] - matrix(x$true_beta, nrow = length(x$true_beta) ,ncol = B, byrow = FALSE)
  return(sum(apply(biases,2,FUN = function(y){(norm(y,type = '2'))^2}))/(B * p*d) )
  }
))

dat_gaussian$mse_rho = c(sapply(gaussian_list,FUN = function(x){
  B = dim(x$raw)[2];
  K = length(x$true_rho)
  d = length(x$true_beta) / x$p
  p = x$p
  biases = x$raw[((p*d) + p + 2):((p*d) + p + K),] - matrix(x$true_rho[2:K], nrow = K-1 ,ncol = B, byrow = FALSE)
  return(sum(apply(biases,2,FUN = function(y){(norm(y,type = '2'))^2}))/(B * (K-1)) )
}
))

#To format the axis labels with a given number of decimals
fmt_dcimals <- function(decimals=0){
  function(x) sprintf("%0.3f",round(x,decimals))
}

gaussian_mse_beta_plot = ggplot(data=dat_gaussian, aes(x=n, y=mse_beta,col = p, group=p)) +
  geom_line() + geom_point() + 
  ylab(expression(paste('MSE(',beta,')'))) + ggtitle('Gaussian')+ scale_y_continuous(labels =fmt_dcimals(3)) + theme_bw() +theme(plot.title = element_text(hjust = 0.5)) #+ guides(col = "none") 

gaussian_mse_rho_plot = ggplot(data=dat_gaussian, aes(x=n, y=mse_rho,col = p, group=p)) +
  geom_line() + geom_point() + 
  ylab(expression(paste('MSE(',rho,')'))) + ggtitle('Gaussian')+ scale_y_continuous(labels =fmt_dcimals(3)) + theme_bw() +theme(plot.title = element_text(hjust = 0.5)) #+ guides(col = "none") 

ggarrange(gaussian_mse_beta_plot,
          gaussian_mse_rho_plot,ncol = 1, common.legend = T, legend = 'bottom')


#Chi_sq coverage for hatbeta_S with S = {1,2,3,4,5} and hatrho
dat_gaussian$chisq_coverage_beta = sapply(gaussian_list,FUN = function(x){
  x$chisq_beta_est_coverage
})

dat_gaussian$chisq_coverage_rho = sapply(gaussian_list,FUN = function(x){
  x$chisq_rho_est_coverage
})

chisq_coverage_beta_table  = spread(dat_gaussian %>% dplyr::select(n,p, chisq_coverage_beta), key = p, value = chisq_coverage_beta)
chisq_coverage_rho_table = spread(dat_gaussian %>% dplyr::select(n,p, chisq_coverage_rho), key = p, value = chisq_coverage_rho)

colnames(chisq_coverage_beta_table) = colnames(chisq_coverage_rho_table) = c('n','p=10','p=25','p=50')

chisq_coverage_beta_table$n = c('n=50','n=100','n=200','n=400')
chisq_coverage_rho_table$n = c('n=50','n=100','n=200','n=400')

chisq_coverage_beta_table
chisq_coverage_rho_table


#Boxplots summarizing coverage for each of the mean regression coefficients
plot_boxplot_coverage_beta = function(result_list,title){
  B_vec = c(sapply(result_list,function(x){x$total_replication}))
  n_vec = c(sapply(result_list,function(x){dim(x$data_Y)}[3]))
  p_vec = c(sapply(result_list,function(x){dim(x$data_Y)}[1]))
  d = length(result_list[[1]]$true_beta)/p_vec[1]
  
  dat = NULL
  
  for(i in 1:length(B_vec)){
    df = data.frame(n=n_vec[i],p=p_vec[i])
    dat = rbind(dat,df[rep(seq_len(nrow(df)), each = p_vec[i]*d), ])
  }
  
  
  colnames(dat) = c('n','p')
  dat$n = as.factor(dat$n)
  dat$p = as.factor(paste0('p=',dat$p))
  dat$beta_coverage = unlist(c(sapply(result_list,FUN = function(x){l = length(x$CI_coverage); return(x$CI_coverage[1:(l-5)])})))
  
  dat$type = 'OLS'
  
  
  
  
  graph.beta.coverage = ggplot(dat,aes(x=n,y=beta_coverage,fill=type, group = interaction(type,n))) +geom_boxplot() +facet_grid(~p) + guides(fill='none') + ylab(expression(paste('Coverage for elements of ',beta))) +geom_hline(yintercept=0.95, linetype="dashed", color = "red") + ggtitle(title)  + ylim(0.845,0.975) + theme_bw()+ theme(plot.title = element_text(hjust = 0.5))
  
  
  return(graph.beta.coverage)
}

plot_boxplot_coverage_beta(gaussian_list, title = 'Gaussian')

#Coverage for each of the correlation regression parameters
n_string = c('n=50','n=100','n=200','n=400')
gaussian_rho1_coverage = round(matrix(sapply(gaussian_list,function(x){l=length(x$CI_coverage) - 5;return(x$CI_coverage[l+1])}), nrow = 4, byrow = FALSE),3)
gaussian_rho2_coverage = round(matrix(sapply(gaussian_list,function(x){l=length(x$CI_coverage) - 5;return(x$CI_coverage[l+2])}), nrow = 4, byrow = FALSE),3)
gaussian_rho3_coverage = round(matrix(sapply(gaussian_list,function(x){l=length(x$CI_coverage) - 5;return(x$CI_coverage[l+3])}), nrow = 4, byrow = FALSE),3)
gaussian_rho4_coverage = round(matrix(sapply(gaussian_list,function(x){l=length(x$CI_coverage) - 5;return(x$CI_coverage[l+4])}), nrow = 4, byrow = FALSE),3)
gaussian_rho5_coverage = round(matrix(sapply(gaussian_list,function(x){l=length(x$CI_coverage) - 5;return(x$CI_coverage[l+5])}), nrow = 4, byrow = FALSE),3)

rho_element_coverage = rbind(cbind(n_string,gaussian_rho1_coverage),
                             cbind(n_string,gaussian_rho2_coverage),
                             cbind(n_string,gaussian_rho3_coverage),
                             cbind(n_string,gaussian_rho4_coverage),
                             cbind(n_string,gaussian_rho5_coverage))
colnames(rho_element_coverage) = c('n', 'p = 10', 'p = 25', 'p = 50')
rho_element_coverage
