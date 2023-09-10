rm(list=ls())
library(here)
i_am("Simulations/generate_X_and_W_matrices.R")
beetle_env_file = read.csv(here("Applications","beetle_env.csv"))
X_matrix = cbind(1,data.frame(scale(beetle_env_file[,c('pH','Elevation','Management')])))
X_matrix = t(as.matrix(sapply(X_matrix,as.numeric)))

#Generate additional covariate vectors by randomly sampling from available covariate vectors
set.seed(1)
while(dim(X_matrix)[2]<400){
  add_index = sample(87,3,replace = TRUE)
  X_matrix = cbind(X_matrix,c(1,X_matrix[2,add_index[1]],X_matrix[3,add_index[2]],X_matrix[4,add_index[3]]))
}
X = X_matrix

saveRDS(X,file = here("Simulations","X_matrix.rds"))


#
beetle_spe_file = read.csv(here("Applications","beetle_spe.csv"))
beetle_trait_file = read.csv(here("Applications","beetle_trait.csv"))


beetle_spe = data.frame(beetle_spe_file[,2:88],row.names = beetle_spe_file$Code)
beetle_trait = data.frame(beetle_trait_file[,4:8], row.names = beetle_trait_file$CODE)


include.ind = apply(beetle_spe,1,function(x){sum(x>0)}) >= 15
include.ind[c(58,60)] = F
p_ind = sum(include.ind)

quan_trait = beetle_trait[include.ind,'LTL']
quan_trait = scale(quan_trait,center = T,scale = T)

#Generate additional quantitative trait by randomly sampling from available quantitative traits
set.seed(1)
trait_add_index = sample(p_ind,50-p_ind,replace = TRUE)
quan_trait = c(quan_trait,quan_trait[trait_add_index,])


qual_trait = beetle_trait[include.ind,c("CLG","WIN","OVE","BRE")]
qual_trait = sapply(qual_trait,function(x){gsub("\\?","",x)})

#Generate additional qualitative traits by randomly sampling from available qualitative traits
while(dim(qual_trait)[1]<50){
  add_index = sample(p_ind,4,replace = TRUE)
  qual_trait = rbind(qual_trait,qual_trait[add_index[1],1],qual_trait[add_index[2],2],qual_trait[add_index[3],3],qual_trait[add_index[4],4])
}

p = 50; K = 5
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

saveRDS(W,file = here("Simulations","W_matrix.rds"))

