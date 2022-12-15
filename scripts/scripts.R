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
    
  n_groups=length(rho)
  n_elems=length(a_weights)
  unifsample = runif(n=1)
  
  # unifsample>alfaADD --> delete move
  # unifsample<alfaADD --> add move
  # preliminary check for extreme cases
  if((unifsample>alfaADD && n_groups==1) || (unifsample<alfaADD && n_groups==n_elems) ){
    return(0)
  }
  
  #Now I select the element range which we will use to extract the useful indexes
  elems_range=1:n_elem #same size of a_weights and d_weights, i.e., n-1
  changepoint_indexes=cumsum(rho) 
    
  a_weights_nocp=a_weights[!changepoint_indexes]
  a_weights_sum=sum(a_weights_nocp)
  
  d_weights_cp=d_weights[changepoint_indexes]
  d_weights_sum=sum(d_weights_nocp)
 

  if (unifsample<=alfaADD)
  { #I choose ADD - the possible elements are just the ones which are NOT the changepoints, which are all the indexes except the ones from rho
    possible_indexes= elems_range[!changepoint_indexes]
    weight=a_weights_nocp
    weight_sum=a_weights_sum
    probratio=(1-alfaADD)/1 #to use just for the ratio in the case n_groups==1
  }else
  { #I choose DELETE - the possible elements are just the ones which are the changepoints, which correspond to the values in rho
    possible_indexes=elems_range[changepoint_indexes]
    weight=d_weights_cp
    weigth_sum=sum(d_weights_cp)
    probratio=alfaADD/1   #to use just for the ratio in the case n_groups==n_elems
  }
  #I extract the index of the element sampling n=1 element from a categorical distribution which assigns to "possible elements" the vector of probabilities "weight/weight_sum")
  elem_index=sample(possible_elems, 1, replace = TRUE, prob = weight/weight_sum) #samples with replacement from a categorical
  
  #I deal with the remining case where I have all changepoints or zero changepoints
  if (n_groups==n_elems || n_groups==1){ # I have no choice on delete or merge
    ratio= (weight_sum/weight[elem_index])*(probratio) 
    return (ratio)
  }
  if(unifsample<=alfaADD)
  {ratio=  (1-alfaADD)/alfaADD * a_weight_sum/a_weights[elem_index]*(d_weights[elem_index]/(d_weights[elem_index]+d_weight_sum))} #case ADD
  else
  {ratio= alfaADD/ (1-alfaADD) * d_weight_sum/d_weights[elem_index]*(a_weights[elem_index]/(a_weights[elem_index]+a_weight_sum))} #case DELETE 
  return (ratio)
}



#' splitPartition in the compact form
#' 
#' NOTE: maybe it would be bette to switch the representation, add a changepoint and switch back... 
#' let's think about it
#'
#' @param k index of the the point where to split the group (equivalent to adding a changepoint)
#' @param rho_n partition in compact form e.g., rho=c(2,3,4) meands the first group has 2 elements, the second has three and the third has four
#'
#' @return a list whose first element is the updated partition 
#' and the second is the index of the groups that has changed (do not know if it is necessary, though)
#' @export
#'
#' @examples
splitPartition <- function(k,rho_n){
  n_elems=sum(rho_n)
  n_groups=length(rho_n)
  
  if(n_elems==(length(rho_n)) || k>n_elems-1){ #First check: all groups with 1 element or index out of bound (number of changepoints=n_elems-1)
    output[[1]] = rho_n
    output[[2]] = 0 # returns 0 if the change has not been performed
    return (output) 
  }
  cumsum_rho=cumsum(rho_n)
  found=FALSE
  #For every index i, I store the new i-th group
  for(i in 1:ngroups){
    if(found==FALSE && cumsum[i]==k) #k is already a changepoint, nothing to split - returns the original rho
    {
      output[[1]] = rho_n
      output[[2]] = 0 # returns 0 if the change has not been performed
      return (output) 
    }
    #I update the partition in the general case (either I have already split the group or not, just the indexing changes)
    if(found==FALSE) {
        new_rho[i]=rho_n[i]
    }
    else {
        new_rho[i+1]=rho_n[i]
    }
    if (found==FALSE && cumsum[i]>k){ #in the case I have just passed the element index - e.g. i am in the group to be split
      new_rho[i]=k-cumsum_rho[i-1]  # is the index of the element minus the cumulative number of elements in the previous groups
      new_rho[i+1]=rho[i]-new_rho[i] # is the dimension of the original group minus the elements moved to new_rho[i]
      j=i # I save the index of the group that has changed (the i-th group)- not sure it is necessary, though
      found=TRUE
    }
  }
  output[[1]] = new_rho
  output[[2]] = j
  return(output)
}



#Useful Functions

#' Computes the rising factorial (also called Pochammer symbol), eventually in logarithmic form.
#' 
#' @param x the value for which the rising factorial is to be calculated
#' @param n the power x is to be raised to
#' @param log if true, the rising factorial is calculated in logarithmic form
#'
#' @return the rising factorial of x to the n power
#' @export
#'
#' @examples
pochhammer <- function(x,n,log=FALSE){   
  if(log==TRUE){
    num_vec <- as.numeric()
    num_vec[1] = x
  
    if (n!=0){
      for(i in 1:(n-1)){num_vec[i+1] <- (x+i)} 
    }
  
    if (n==0){num_vec[1] <- 1}
  
    num <- sum(log(num_vec))
  
    return(num)
  }
  
  else{
    num_vec <- as.numeric()
    
    num_vec[1] = x
    
    if (n!=0){
      for(i in 1:(n-1)){
        num_vec[i+1] <- (x+i)
      } 
    }
    
    if (n==0){num_vec[1] <- 1}
    
    num <- prod(num_vec)
    
    return(num)
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
