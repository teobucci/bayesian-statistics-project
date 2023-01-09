library(logr)
suppressWarnings(suppressPackageStartupMessages(library(tidyverse)))
suppressWarnings(suppressPackageStartupMessages(library(ACutils)))
suppressWarnings(suppressPackageStartupMessages(library(mvtnorm)))
suppressWarnings(suppressPackageStartupMessages(library(salso)))

# copy and paste from functions.R
Generate_data = function(n,p,z=NULL,W=NULL,seed=123421){
    
    nu = max(n,p) # default choice of nu
    if(n<p)
        stop("p can not be larger than n")
    
    if(is.null(z)){ #default is 3 groups
        if(p<3)
            stop("p must be greater than 3 if z is defaulted")
        z = vector(length = p)
        z[1:ceiling(p/3)] = 1
        z[(ceiling(p/3)+1):(2*ceiling(p/3))] = 2
        z[(2*ceiling(p/3)+1):p] = 3
    }else{
        if(p != length(z))
            stop("the length of z should be p")
    }
    
    if(is.null(W)){ #default is to sample a W from Wishart(nu, I_p)
        D = diag(p)
        set.seed(seed)
        W = MCMCpack::rwish(nu,D)
    }
    
    # Generate V
    V = ComputeV(nu,W,z)
    
    # Generate Omega
    Omega = MCMCpack::rwish(nu,V)
    
    # Compute Sigma (inefficient)
    Sigma = solve(Omega)
    
    # Generate data
    data = matrix(NA,nrow = n, ncol = p)
    data = t(apply(data, 1, FUN = function(x){ mvtnorm::rmvnorm(n = 1, sigma = Sigma)  }))
    
    #return 
    return(list("data" = data,
                "Omega" = Omega,
                "Sigma" = Sigma,
                "V"     = V,
                "z_true"= z
    ))
}


ComputeV = function(nu,W,z){
    p = dim(W)[1] #get dimension
    if(length(z)!=p)
        stop("length of z is different from nrow(W)")
    if(nu<=0)
        stop("nu can not be negative or zero")
    
    V = matrix(0,nrow = p,ncol = p) 
    # fill the upper triangular part of V
    for(i in 1:p){
        for(j in (i+1):p){
            if(j>p)break; #needed because r for loops sucks
            if(z[i]==z[j])
                V[i,j] = W[i,j]/nu #just see the definition
        }
    }
    V = V + t(V) #set the lower trg part
    diag(V) = diag(W)/nu #set the diagonal
    return (V)
}

#W = solve(cov(data))
#nu = max(n,p)
library("FGM")

utility_path = file.path("src/utility_functions.R")
file.exists(utility_path)
source(utility_path)

bulky_path = file.path("src/bulky_functions.R")
file.exists(bulky_path)
source(bulky_path)

data_generation_path = file.path("src/data_generation.R")
file.exists(data_generation_path)
source(data_generation_path)


z_true = c(1,1,1,2,2,3,3,3,4,4) # these are the groups memberships (10 nodes)
counts_true = as.vector(table(z_true))
Nclust_true  = length(counts_true)

p = length(z_true)
n = 100
data_from_prior = Generate_data(n=n,p=p,z=z_true,seed=27091999)
nu = n # why?
data = data_from_prior$data
Omega_true = data_from_prior$Omega # Omega is Kappa (precision matrix TODO change)
Sigma_true = data_from_prior$Sigma
V = data_from_prior$V # cos'è?
U = t(data)%*%data

options = set_options(sigma_prior_0=0.5,
                      sigma_prior_parameters=list("a"=1,"b"=1,"c"=1,"d"=1),
                      theta_prior_0=1,
                      theta_prior_parameters=list("c"=1,"d"=1),
                      rho0=c(5,5),
                      weights_a0=rep(1,p-1),
                      weights_d0=rep(1,p-1),
                      alpha_target=0.40,
                      mu_beta=0.5, # mu of the Beta
                      sig2_beta=0.2, # variance of the Beta
                      d=3,
                      alpha_add=0.5,
                      adaptation_step=0.5,
                      update_sigma_prior=F,
                      update_theta_prior=F,
                      update_weights=T,
                      update_partition=T,
                      update_graph=T,
                      perform_shuffle=T)

nburn = 1
niter = 2500
nburn = 1
thin = 1 # tieni tutti i multipli di thin

log_open(file_name = "log_aiuto2", show_notes=FALSE)
res = Gibbs_sampler(data,niter,nburn,thin,options,seed=270999,print=TRUE)
log_close()