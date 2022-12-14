UpdatePartition = function(z,counts,Nclust,alpha,MARGINAL=NULL, ...){
    
    marginal_params = list(...) #get list with parameters to be used in MARGINAL function
    n_params = length(marginal_params) #get number of parameters to be used in MARGINAL function
    
    if(is.null(MARGINAL)){
        MARGINAL = function(x){1};
    }
    p = length(z)
    
    
    # parameters to be set in the function header
    alpha_add = 0.5
    r = c(0,0,0,0,0) # con p colonne
    
    
    # check if you want to perform an add or delete move
    perform_add = F
    if(all(r==0)){ # must perform add move
        perform_add = T
    }else if(all(r==1)){ # must perform delete move
        perform_add = F
    }else{
        if(runif(n=1)<alpha_add) # uniform draw
            perform_add = T # add
        }else{
            perform_add = F # delete (not necessary this step, kept for clarity)
        }
    }
    
    
    
    

#' Adaptation step to update the weights vector
#' The function takes as an input the current weights and update them as a function of
#' the current iteration number t, the initial adaptation h, 
#' The function works in log scale
#' 
#' @param logweights vector of the logaritm of the current weights
#' @param alfaTarget scalar indicating the target acceptance probability (optimal range around 0.10-0.15)
#' @param alfaADD probability of adding a move (usually 0.5)
#' @param t number of the current iteration
#' @param h initial adaptation (must be >0)
#'
#' @return the vector of updated logweights
#' @export
#'
#' @examples
#' 
logAdaptation = function(logweights,t,h,alfaTarget,alfaADD){ 
  if(h>0){ #checks that the adaptation step is positive
    return (logweights + h*length(logweights)/t *(alfaADD - alfaTarget))
  }
}


#' Partition current data
#' Given y (vector of data) and rho_n (vector of the partition) the function splits the observations and
#' partitions them into the current groups 
#' partitions them into the partition rho
#' @param y - vector of n ordered data
#' @param rho_n - partition written in compact form 
#' e.g. rho=c(2,4,5) means the first group has 2 elements, the second has 4 and the fifth has 5
#'
#' @return a list, where each element contains the corresponding group of y elements.
#' If the dimensions of rho and y are not comparable, return an empty vector
#' @export
#'
#' @examples
partitionData <- function(y,rho_n){ 
  if(sum(rho_n)==length(y)){ #checks that the dimensionality is okay
  
    dataPartition <- list()
    
    for (i in 1:length(rho_n)){
      if(i == 1) 
        {
          dataPartition[[i]] <- y[1:rho_n[i]]
          cumsum=rho_n[i]
        }
      else
       {
         first_index=cumsum + 1
         last_index=rho_n[i] + cumsum
         dataPartition[[i]] <- y[first_index : last_index]
         cumsum=cumsum+rho_n[i]
       }
    }
    return(dataPartition) 
    }
}














    
    
    
    
    
    
    
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
