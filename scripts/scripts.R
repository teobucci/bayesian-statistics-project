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
    
    
    
    

#' Adaptation step to update the weights vector a and d
#' The function takes as an input the current weights and updates them as a function of
#' the current iteration number t, the initial adaptation h, 
#' The function works in log scale
#' 
#' @param logweights vector of the logarithm of the current weights
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
    if(!h>0)
        stop("Adaptation step h must be positive")
    return (logweights + h*length(logweights) / t * (alfaADD - alfaTarget))
}


#' Partition current data
#' Given y (vector of data) and rho_n (vector of the partition) the function splits the observations and
#' partitions them into the current groups 
#' partitions them into the partition rho
#' @param y - vector of n ordered data
#' @param rho_n - partition written in compact form 
#' e.g. rho=c(2,4,5) means the first group has 2 elements, the second has 4 and the third has 5
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



#' Computes the proposal ratio from the file "Samplingstrategy_nonparam" at page 4
#' Follows strictly the paper, just adds the possibility of having alfaADD set by the author
#' NOTE: I have not added the checks for the range of the parameters ranging from 0 and 1. Anyway, they are easy to include in the first if
#'
#' @param rho the partition in compact form (e.. rho=c(1,4,5) means that the first group has 1 element, the second has 4 elements, the last has 5 elements)
#' @param alfaADD fixed probability of choosing ADD or DELETE as a move
#' @param a_weights vector of size n-1 (number of datapoints - 1) containing at element j the weights to consider when ADDING a changepoint between point j and point j+1 (weights are non-normalized probabilities)
#' @param d_weights vector of size n-1 (number of datapoints - 1) containing at element j the weights to consider when DELETING a changepoint between point j and point j+1 (weights are non-normalized probabilities)
#'
#' @return the proposal ratio, which is 
#' @export
#'
#' @examples
proposalRatio=function(rho, alfaADD, a_weights, d_weights){
    
    #alfaADD=0.5
  num_groups = length(rho)
  n_elems = length(a_weights)
  unifsample = runif(n=1)
  choose_add = unifsample<alfaADD
  
  # preliminary check for extreme cases
  if((!choose_add && num_groups==1) || (choose_add && num_groups==n_elems) ){
      # incompatible to delete when only one group is present
      # or to add when every point is a group
      ratio=0
      return(ratio)
  }
  
  # select the element range which we will use to extract the useful indexes
  elems_range <- 1:n_elems # same size of a_weights and d_weights
  
  # indexes of the changepoints
  cp_indexes <- cumsum(rho)
  
  # PER ORA SUPPONGO CHE A_WEIGHTS SIA DELLA STESSA DIMENSIONE N, NON N-1
  # QUELLO MANCANTE SI SETTA A MANO
  
  #rho
  #cumsum(rho)
  #a_weights
  #n_elem=length(a_weights)
  
  
  a_weights_available = a_weights
  a_weights_available[cp_indexes] = 0
  a_weights_available_sum = sum(a_weights_available)
  
  d_weights_available = d_weights
  d_weights_available[-cp_indexes] = 0
  d_weights_available_sum = sum(d_weights_available)
  
  
  if (choose_add){
      draw_weights = a_weights_available
  } else {
      draw_weights = d_weights_available
  }
  
  candidate = sample(1:n_elem, 1, prob=draw_weights)
  
  if (choose_add && num_groups==1){
      # case in which you choose to propose an add move
      # (with just 1 group) that may or may not be accepted
      ratio = alfaADD/1 * a_weights_available_sum / a_weights[candidate]
      return(ratio)
  }
  
  if (!choose_add && num_groups==n_elems){
      # case in which you choose to propose an delete move
      # (with every point being a group) that may or may not be accepted
      ratio = alfaADD/1 * d_weights_available_sum / d_weights[candidate]
      return(ratio)
  }
  
  # only the general cases remain
  if (choose_add){
      ratio = (1-alfaADD)/alfaADD * a_weights_available_sum / a_weights[candidate] * d_weights[candidate] / (d_weights[candidate] + d_weights_available_sum)
      return (ratio)
  } else {
      ratio = alfaADD/ (1-alfaADD) * d_weights_available_sum / d_weights[candidate] * a_weights[candidate] / (a_weights[candidate] + a_weights_available_sum)
      return (ratio)
  }
  
}






likelihoodRatio=function(rho, alfaADD, a_weights, d_weights){
    
    alpha = 1
    beta = 1
    M = length(rho)
    
    # auxiliary function to evaluate the beta
    function rhoB(arg1,arg2,log=T){
        if(log){
            return(lbeta(alpha + arg1, beta + arg2))
        }else{
            return(beta(alpha + arg1, beta + arg2))
        }
        
    }
    
    if("il nodo non è agli estremi")

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
