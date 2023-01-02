# library("MCMCpack")

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

UpdatePrecision = function(nu,W,n,U,z){
    p = dim(W)[1]
    if(length(z)!=p)
        stop("length of z is different from nrow(W)")
    if(dim(U)[1]!=p)
        stop("nrow(U) is different from nrow(W)")
    if(nu<=0)
        stop("nu can not be negative or zero")
    if(n<=0)
        stop("nu can not be negative or zero")
    
    dof = n + nu #dof a posterior
    V = ComputeV(nu,W,z) #deterministically compute V
    invScale = solve(V) + U #define inverse scale matrix a posteriori
    
    # MCMCpack::rwish ---> requires dof and inverse scale matrix
    # WATCH OUT,in this package there is an error with the parametrization. They say
    # "The mean of a Wishart random variable with v degrees of freedom and inverse scale matrix S is vS",
    # but this is not true, what they call inverse scale matrix is actually a scale matrix
    
    Omega_new = MCMCpack::rwish(dof,solve(invScale) )  
    return (Omega_new)
}

log_eval_Wishdens = function(Omega,z,nu,W){
    V = ComputeV(nu,W,z)
    Vinv = solve(V)
    p = dim(Omega)[1]
    
    return( 0.5*(nu-p-1)*log(det(Omega)) -
                0.5*sum(diag(Vinv%*%Omega)) - 
                0.5*(nu*p)*log(2) - 
                0.5*nu*log(det(V)) -
                Boom::lmgamma(nu/2, p)  
    )
    
    
}

marginal_Wishartdes = function(z,counts,Nclust,marginal_params){
    if(length(marginal_params)!=3)
        stop("3 parameters are needed: Omega,nu,W ")
    
    Omega = marginal_params[[1]]
    nu    = marginal_params[[2]]
    W     = marginal_params[[3]]
    return( log_eval_Wishdens(Omega,z,nu,W) )
}

#' UpdatePartition - Chinese Restaurant Process
#'
#' This function implements the CRP induced by the Dirichlet process with parameter \code{\alpha}. 
#' This is for the conjugate
#' @param z [vector] the current partition
#' @param counts [vector] number of elements in each cluster
#' @param Nclust [integer] total number of clusters
#' @param alpha [double] DP parameter
#' @param MARGINAL generic function defining the marginal law of a single data. Set NULL for the prior case
#' @param ... parameters needed by \code{MARGINAL} function
#'
#' @export
UpdatePartition = function(z,counts,Nclust,alpha,MARGINAL=NULL, ...){
    
    marginal_params = list(...) #get list with parameters to be used in MARGINAL function
    n_params = length(marginal_params) #get number of parameters to be used in MARGINAL function
    
    if(is.null(MARGINAL)){
        MARGINAL = function(x){1};
    }
    p = length(z)
    # take a Gibbs step at each data point
    for(j in 1:p) {
        # get rid of the jth data point
        c = z[j]
        counts[c] = counts[c] - 1
        # if the jth data point was the only point in a cluster,
        # get rid of that cluster
        if(counts[c]==0) {
            counts[c] = counts[Nclust] #move last cluster in position c
            loc_z = (z==Nclust)
            z[loc_z] = c #remove cluster in position Nclust
            counts = counts[-Nclust]
            Nclust = Nclust - 1
        }
        z[j] = -1  # ensures z[j] doesn't get counted as a cluster
        
        # unnormalized log probabilities for the clusters
        log_weights = rep(NA,Nclust+1)
        # find the unnormalized log probabilities
        # for each existing cluster
        for(c in 1:Nclust) {
            local_z = z
            local_z[j] = c
            log_weights[c] = log(counts[c])  + MARGINAL(local_z,counts,Nclust,marginal_params)
            #log_eval_Wishdens(Omega,local_z,nu,W ) #- log(alpha+p-1)
        }
        # find the unnormalized log probability
        # for the "new" cluster
        local_z = z
        local_z[j] = Nclust+1
        log_weights[Nclust+1] = log(alpha) + MARGINAL(local_z,counts,Nclust,marginal_params)
        #log_eval_Wishdens(Omega,local_z,nu,W ) #- log(alpha+p-1)
        
        # transform unnormalized log probabilities
        # into probabilities
        max_weight = max(log_weights)
        log_weights = log_weights - max_weight
        loc_probs = exp(log_weights)
        loc_probs = loc_probs / sum(loc_probs)
        
        if(any(is.na(loc_probs)))
            stop("NA values in UpdatePartition")
        # sample which cluster this point should
        # belong to
        newz = sample(1:(Nclust+1), 1, replace=TRUE, prob=loc_probs)
        # if necessary, instantiate a new cluster
        if(newz == Nclust + 1) {
            counts = c(counts,0)
            Nclust = Nclust + 1
        }
        z[j] = newz
        # update the cluster counts
        counts[newz] = counts[newz] + 1
    }
    
    return(list("z" = z,"counts" = counts,"Nclust" = Nclust))
    
}

log_FC_alpha = function(alpha,K,p,a_alpha,b_alpha){
    if(alpha <= 0)
        stop("Error in log_FC_alpha, alpha can not be zero")
    #K is the number of clusters
    return( K*log(alpha) + lgamma(alpha) - lgamma(alpha+p) + 
                dgamma(x=alpha, shape = a_alpha, rate = b_alpha, log = T) 
    )
}

Updatealpha_MH = function(alpha_old,K,p,a_alpha,b_alpha, var_alpha_adp, iter, adaptive = T){
    
    log_alpha_new = rnorm(n = 1, mean = log(alpha_old), sd = sqrt(var_alpha_adp))  #draw proposed value in log scale
    alpha_new = exp(log_alpha_new) #compute proposed value in standard scale
    
    # compute acceptance probability in log scale
    ln_acp = log_FC_alpha(alpha_new,K,p,a_alpha,b_alpha) - log_FC_alpha(alpha_old,K,p,a_alpha,b_alpha) +
        alpha_new - alpha_old
    
    if(is.na(ln_acp))
        stop("Error in Updatealpha_MH, ln_acp is NA")
    
    # compute acceptance probability in normal scale
    acp = exp(ln_acp)
    
    # acceptance-rejection step
    u = runif(n=1) 
    if( u < min(1,acp) )
        res = alpha_new
    else
        res = alpha_old
    
    # adaptive variance step
    ww_g = (iter+1)^(-0.7) 
    tau_bar = 0.44
    if(adaptive){
        
        var_alpha_adp_new = var_alpha_adp * exp( ww_g *
                                                     (  exp(min(0,ln_acp)) - tau_bar  )
        )
        
    }else{
        var_alpha_adp_new = var_alpha_adp
    }
    
    
    
    if(var_alpha_adp_new < 1e-10)
        var_alpha_adp_new = 1e-10
    
    if(var_alpha_adp_new > 1e+10)
        var_alpha_adp_new = 1e+10
    
    return(list("alpha"=res,"var_alpha_adp"=var_alpha_adp_new))
    
}

Updatealpha_augmentation = function(alpha_old,K,p,a_alpha,b_alpha){
    
    # Data augmentation sampling
    U = rbeta(n=1,shape1 = alpha_old + 1, shape2 = p) # draw auxiliary variable
    if(U <= 0)
        stop("Error in Updatealpha_augmentation, U can not be zero")
    
    # Exact sampling for alpha
    weight = (K+a_alpha-1)/(p*(b_alpha - log(U))) # compute mixture weight
    
    # choose mixture component
    rnd = runif(n=1) 
    if(rnd < weight){
        shape = K + a_alpha
        rate  = b_alpha - log(U)
        return( rgamma(n=1, shape = shape, rate = rate) )
    }else{
        shape = K + a_alpha - 1
        rate  = b_alpha - log(U)
        return( rgamma(n=1, shape = shape, rate = rate) )
    }
}

set_options = function( nu, W, a_alpha = 1, b_alpha = 1,
                        Omega0, z0, alpha0,
                        var_alpha_adp0 = 1, adaptiveAlpha = T,
                        UpdatePartition = T, UpdateOmega = T, UpdateAlpha = T){
    
    option = list("nu"=nu, "W"=W, "a_alpha" = a_alpha, "b_alpha" = b_alpha,
                  "Omega0"=Omega0, "z0"=z0, "alpha0"=alpha0,
                  "var_alpha_adp0"=var_alpha_adp0, "adaptiveAlpha" = adaptiveAlpha,
                  "UpdatePartition" = UpdatePartition, "UpdateOmega" = UpdateOmega, "UpdateAlpha" = UpdateAlpha)
    return (option)
}

Gibbs_sampler = function(data,niter,nburn,thin,
                         options,
                         seed=1234,print=T)
{
    n = nrow(data)
    p = ncol(data)
    n_total_iter =  nburn + niter*thin
    
    if(length(options$z0)!=p)
        stop("length p0 not coherent with ncol(data)")
    if(nrow(options$W)!=p)
        stop("nrow W not coherent with ncol(data)")  
    if(nrow(options$Omega0)!=p)
        stop("nrow Omega0 not coherent with ncol(data)")
    
    # get initial values
    Omega = options$Omega0
    z     = options$z0
    alpha = options$alpha0
    var_alpha_adp = options$var_alpha_adp0
    # other quantities
    U = t(data)%*%data # compute U (data must have zero mean)
    counts = as.vector(table(z))  # initial data counts at each cluster
    Nclust = length(counts)	  # initial number of clusters
    
    # Define structure to save sampled values
    save_res = vector("list",length = 4)
    names(save_res) = c("Omega","Partition","alpha","var_alpha_adp")
    save_res$Omega = vector("list",length = niter) 
    save_res$Partition = matrix(NA,nrow = niter, ncol = p)
    save_res$alpha = rep(NA,niter)
    save_res$var_alpha_adp = rep(NA,niter)
    
    it_saved = 0 #initialize iteration counter
    pb = txtProgressBar(min = 1, max = n_total_iter, initial = 1, style = 3)  #initialize progress bar
    for(iter in 1:n_total_iter){
        
        # Update precision matrix
        if(options$UpdateOmega)
            Omega = UpdatePrecision(options$nu,options$W,n,U,z)
        
        # Update partition
        if(options$UpdatePartition){
            # cat('\n pre z = ', z, '\n')
            # cat('\n pre counts = ', counts, '\n')
            # cat('\n pre Nclust = ', Nclust, '\n')
            list_update_part = UpdatePartition(z,counts,Nclust,alpha,
                                               MARGINAL = marginal_Wishartdes,
                                               Omega,options$nu,options$W)
            z = list_update_part$z
            counts = list_update_part$counts
            Nclust = list_update_part$Nclust
            # cat('\n post z = ', z, '\n')
            # cat('\n post counts = ', counts, '\n')
            # cat('\n post Nclust = ', Nclust, '\n')
        }
        
        # Update alpha 
        if(options$UpdateAlpha){
            #alpha = Updatealpha_augmentation(alpha,Nclust,p,options$a_alpha,options$b_alpha)
            list_update_alpha = Updatealpha_MH(alpha,Nclust,p,options$a_alpha,options$b_alpha, var_alpha_adp, iter, options$adaptiveAlpha)
            alpha = list_update_alpha$alpha
            var_alpha_adp = list_update_alpha$var_alpha_adp
        }
        
        # save results
        if(iter>nburn && (iter-nburn)%%thin == 0){
            it_saved = it_saved + 1 
            save_res$Omega[[it_saved]] = Omega
            save_res$Partition[it_saved,] = z
            save_res$alpha[it_saved] = alpha
            save_res$var_alpha_adp[it_saved] = var_alpha_adp
        }
        
        if(print){
            setTxtProgressBar(pb, iter)
        }
        
    }
    
    close(pb)
    return(save_res)
}



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
    
    # Genetate Omega
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


Prior_expected_n_clust = function(alpha,n){
    alpha * log(1+n/alpha)
}









