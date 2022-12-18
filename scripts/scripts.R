UpdatePartition = function(z,counts,Nclust,alpha,MARGINAL=NULL, ...){
    
    
    
    # sampling dalla partizione
    # metropolis hastings per aggiornare la partizione
    
    likelihoodRatioNow = likelihoodRatio
    priorRatioNow = priorRatio
    proposalRatioNow = proposalRatio
    
    alpha_accept <- min(1,likelihoodRatioNow*priorRatioNow*proposalRatioNow)
    if (runif(n=1) < alpha_accept){
        # accept the move
        if(merge){ # accepted move is a merge
            mergePartition
        }
        else{ # accepted move is a split
            splitPartition
        }
        partitionData
    }
    
    
    
    
    
    
    
    
    
    
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
          cumsum_rho=rho_n[i]
        }
      else
       {
         first_index=cumsum_rho + 1
         last_index=rho_n[i] + cumsum_rho
         dataPartition[[i]] <- y[first_index : last_index]
         cumsum_rho=cumsum_rho+rho_n[i]
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



#' splitPartition in the compact form
#' 
#' NOTE: maybe it would be better to switch the representation, add a changepoint and switch back... 
#' let's think about it
#'
#' @param k index of the the point where to split the group (equivalent to adding a changepoint)
#' @param rho_n partition in compact form e.g., rho=c(2,3,4) means the first group has 2 elements, the second has three and the third has four
#'
#' @return a list whose first element is the updated partition 
#' and the second is the index of the groups that has changed (do not know if it is necessary, though)
#' @export
#'
#' @examples
splitPartition <- function(k,rho_n){
  n_elems=sum(rho_n)
  n_groups=length(rho_n)
  output=list()
  new_rho={}

  
  if(n_elems==(length(rho_n)) || k>n_elems-1){ #First check: all groups with 1 element or index out of bound (number of changepoints=n_elems-1)
    output[[1]] = rho_n
    output[[2]] = 0 # returns 0 if the change has not been performed
    return (output) 
  }
  cumsum_rho=cumsum(rho_n)
  found=FALSE
  #For every index i, I store the new i-th group
  for(i in 1:n_groups){
    if(found==FALSE && cumsum_rho[i]==k) #k is already a changepoint, nothing to split - returns the original rho
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
    if (found==FALSE && cumsum_rho[i]>k){ #in the case I have just passed the element index - e.g. i am in the group to be split
      new_rho[i]=k-(i!=1)*cumsum_rho[i-1*(i!=1)]  # is the index of the element minus the cumulative number of elements in the previous groups only if i!=1
      new_rho[i+1]=rho_n[i]-new_rho[i] # is the dimension of the original group minus the elements moved to new_rho[i]
      j=i # I save the index of the group that has changed (the i-th group)- not sure it is necessary, though
      found=TRUE
    }
  }
  output[[1]] = new_rho
  output[[2]] = j
  return(output)
}





#' mergePartition in the compact form
#' 
#' NOTE: maybe it would be better to switch the representation, add a changepoint and switch back... 
#' let's think about it
#' 
#' @param k index of the the point where to split the group (equivalent to adding a changepoint)
#' @param rho_n partition in compact form e.g., rho=c(2,3,4) means the first group has 2 elements, the second has three and the third has four
#'
#' @return a list whose first element is the updated partition 
#' and the second is the index of the groups that has changed (do not know if it is necessary, though)
#' @export
#'
#' @examples
mergePartition <- function(k,rho_n){
  n_elems=sum(rho_n)
  n_groups=length(rho_n)
  output=list()
  new_rho={}
  
  if(( length(rho_n)==1) || k>n_elems-1){ #First check: only 1 group or index out of bound (number of changepoints=n_elems-1)
    output[[1]] = rho_n
    output[[2]] = 0 # returns 0 if the change has not been performed
    return (output) 
  }
  cumsum_rho=cumsum(rho_n)
  found=FALSE
  #For every index i, I store the new i-th group
  for(i in 1:(n_groups-1)){ #ACHTUNG! Must stop at n_groups - 1 
    if(found==FALSE && cumsum_rho[i]!=k){ #k is already a non-changepoint, nothing to merge - returns the original rho
      output[[1]] = rho_n
      output[[2]] = 0 # returns 0 if the change has not been performed
      return (output) 
      }
    if(found==FALSE){
      new_rho[i]=rho_n[i]
      }
    else {
      new_rho[i]=rho_n[i+1]
      }
    if (found==FALSE && cumsum[i]==k){ #in the case I am at the changepoint between the two groups to be merged
      new_rho[i]=rho_n[i]+rho_n[i+1]  # is the index of the element minus the cumulative number of elements in the previous groups
      j=i # I save the index of the group that has changed (the i-th group)- not sure it is necessary, though
      found=TRUE
      }
    }
    output[[1]] = new_rho
    output[[2]] = j
    return(output)
}

#' shufflePartition (QUESTA RESTA MOLTO SIMILE A CORRADIN, NO?)
#' 
#' @param k index of the the point where the shuffling between groups happens
#' @param rho_n partition in compact form e.g., rho=c(2,3,4) means the first group has 2 elements, the second has three and the third has four
#'
#' @return the updated partition
#' @export
#'
#' @examples
shuffle <- function(k,rho_n){
  new_rho={}
  if(length(rho_n)==1){
    new_rho=rho_n #the shuffling obviously can be done only if the number of groups is at least 2
  }
  else{
  i <- sample(1:(k-1),1)
  j <- sample(1:(rho_n[i] + rho_n[i + 1] - 1),1)
  new_rho[i+1] <- (rho_n[i] + rho_n[i+1] - j)
  new_rho[i] <- j 
  }
  return(new_rho)
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
pochhammer <- function(x,n,log=F){
    if(n<0)
        stop("Non si fa")
  if(log){
    num_vec <- as.numeric()
    num_vec[1] = x
  
    if (n!=0){for(i in 1:(n-1)){num_vec[i+1] <- (x+i)}}
  
    if (n==0){num_vec[1] <- 1}
  
    num <- sum(log(num_vec))
  
    return(num)
  }
  
  else{
    num_vec <- as.numeric()
    
    num_vec[1] = x
    
    if (n!=0){for(i in 1:(n-1)){num_vec[i+1] <- (x+i)}}
    
    if (n==0){num_vec[1] <- 1}
    
    num <- prod(num_vec)
    
    return(num)
  }
}


likelihoodRatio=function(rho, alfaADD, a_weights, d_weights){
    
    alpha = 1
    beta = 1
    M = length(rho)
    
    # auxiliary function to evaluate the beta
    rhoB = function(arg1,arg2,log=T){
        if(log){
            return(lbeta(alpha + arg1, beta + arg2))
        }else{
            return(beta(alpha + arg1, beta + arg2))
        }
        
    }
    
    if("il nodo non è agli estremi"){
        
        # ipotizzo la situa in cui aggiungo un cp NEL MEZZO
        C=c(20,20,5)
        C_star=c(20,7,13,5)
        
        M = length(C)
        S = 2 # è quello che è stato splittato
        # due nuovi sono in S e S+1 (2 e 3)
        M, M+1
        
        # bisogna gestire S_lm e S^star_lm
        
        ratio = -(M+1)*rhoB(0,0)
        
        for(l in 1:(S-1)){
            ratio = ratio + rhoB(C_l_star,C_S_star) + rhoB(C_l_star,C_S+1_star) # primo termine numeratore
            ratio = ratio - rhoB(C_l,C_S) # primo termine denominatore
        }
        
        for(m in (S+2):(M+1))
            ratio = ratio + rhoB(C_S_star,C_m_star) + rhoB(C_S+1_star,C_m_star) # secondo termine numeratore
        
        # terzo termine numeratore
        ratio = ratio + rhoB(C_S_star,C_S+1_star) + rhoB(C_S_star,C_S_star) + rhoB(C_S+1_star,C_S+1_star)
        
        for(m in (S+1):M)
            ratio = ratio - rhoB(C_S,C_m) # secondo termine denominatore
        
        # terzo termine denominatore
        ratio = ratio - rhoB(C_S,C_S)
        
        return(exp(ratio))
    }

}



priorRatio = function(theta, sigma, current_rho, proposed_rho){
    
    #current_rho = c(2,4,2,2)
    #proposed_rho = c(2,4,1,1,2)
    M = length(current_rho)
    
    # ORA È IMPLEMENTATO SOLO IL CASO ADD (CIOÈ SPLIT)
    
    current_r = rep(0, sum(current_rho))
    current_r[cumsum(current_rho)] = 1
    
    proposed_r = rep(0, sum(proposed_rho))
    proposed_r[cumsum(proposed_rho)] = 1
    
    tau = which.max(proposed_r-current_r)
    
    indexes = cumsum(current_rho)
    temp = which(indexes < tau)
    S = temp[length(temp)] + 1
    
    n_star_s = proposed_rho[S]
    n_star_s_plus_1 = proposed_rho[S+1]
    n_s = n_star_s + n_star_s_plus_1
    
    
    ratio = - log(M) + (theta + M*sigma) + pochhammer(1-sigma, n_star_s - 1, log=T)
                                         + pochhammer(1-sigma, n_star_s_plus_1 - 1, log=T)
                                         - pochhammer(1-sigma, n_s - 1, log=T)
                                         + lfactorial(n_s)
                                         - lfactorial(n_star_s)
                                         - lfactorial(n_star_s_plus_1)
    return(ratio)

}

    
    






# set_options = function( nu, W, a_alpha = 1, b_alpha = 1,
#                         Kappa0, z0, alpha0,
#                         var_alpha_adp0 = 1, adaptiveAlpha = T,
#                         UpdatePartition = T, UpdateKappa = T, UpdateAlpha = T){
#     
#     option = list("nu"=nu,
#                   "W"=W,
#                   "a_alpha"=a_alpha,
#                   "b_alpha"=b_alpha,
#                   "Kappa0"=Kappa0,
#                   "z0"=z0,
#                   "alpha0"=alpha0,
#                   "var_alpha_adp0"=var_alpha_adp0,
#                   "adaptiveAlpha"=adaptiveAlpha,
#                   "UpdatePartition"=UpdatePartition,
#                   "UpdateKappa"=UpdateKappa,
#                   "UpdateAlpha"=UpdateAlpha)
#     return (option)
# }


set_options = function(){
    
    options = list(
        "sigma0"           = sigma0,
        "theta0"           = theta0,
        "rho0"             = rho0,
        "alphaAdd"         = alphaAdd,
        "weights_a0"       = weights_a0,
        "weights_d0"       = weights_d0,
        "alpha_target"     = alpha_target,
        "a"                = a,
        "b"                = b,
        "update_sigma"     = update_sigma,
        "update_theta"     = update_theta,
        "update_weights"   = update_weights,
        "update_partition" = update_partition
        )
    return(options)
}

Gibbs_sampler = function(data,niter,nburn,thin,
                         options,
                         seed=1234,print=T)
{
    n = nrow(data) # number of observations
    p = ncol(data) # number of nodes
    n_total_iter = nburn + niter*thin # total iterations to be made
    
    # TODO check qui
    # if(length(options$z0)!=p)
    #     stop("length p0 not coherent with ncol(data)")
    # if(nrow(options$W)!=p)
    #     stop("nrow W not coherent with ncol(data)")  
    # if(nrow(options$Kappa0)!=p)
    #     stop("nrow Kappa0 not coherent with ncol(data)")
    
    # dynamic parameters
    sigma     = options$sigma0 # initial parameter of the nonparametric prior
    theta     = options$theta0 # initial parameter of the nonparametric prior
    rho       = options$rho0 # initial partition (eg. c(150,151))
    weights_a = options$weights_a0 # add weights
    weights_d = options$weights_d0 # del weights
    
    # constant parameters
    alphaAdd = options$alphaAdd # probability of choosing add over delete
    alpha_target = options$alpha_target # target alpha for adapting weights
    
    a = options$a # parameter for the likelihood of the graph (Beta(a,b))
    b = options$b # parameter for the likelihood of the graph (Beta(a,b))
    
    # TODO add graph parameters
    # REGAZ DIAMOCI UNA MOSSAAAA
    
    
    #   var_alpha_adp = options$var_alpha_adp0
    #   # other quantities
    #   U = t(data)%*%data # compute U (data must have zero mean)
    #   counts = as.vector(table(z))  # initial data counts at each cluster
    #   Nclust = length(counts)	  # initial number of clusters
    #   
    #   # Define structure to save sampled values
    #   save_res = vector("list",length = 4)
    #   names(save_res) = c("Kappa","Partition","alpha","var_alpha_adp")
    #   save_res$Kappa = vector("list",length = niter) 
    #   save_res$Partition = matrix(NA,nrow = niter, ncol = p)
    #   save_res$alpha = rep(NA,niter)
    #   save_res$var_alpha_adp = rep(NA,niter)
    
    # initialize iteration counter
    it_saved = 0
    
    # initialize progress bar
    pb = txtProgressBar(min=1, max=n_total_iter, initial=1, style=3)
    
    for(iter in 1:n_total_iter){
        
        # pseudocode
        
        if('voglio adattare i pesi'){
            logAdaptation
        }
        
        # Update precision matrix
        if(options$UpdateKappa){
            cat("Updating the precision matrix...")
            # TODO aggiorna precision matrix da bdgraph
            Kappa = UpdatePrecision(options$nu,options$W,n,U,z)
        
        }
            
        # Update partition
        if(options$UpdatePartition){
            cat("Updating the partition...")
            # cat('\n pre z = ', z, '\n')
            # cat('\n pre counts = ', counts, '\n')
            # cat('\n pre Nclust = ', Nclust, '\n')
            list_update_part = UpdatePartition(z,counts,Nclust,alpha,
                                               MARGINAL = marginal_Wishartdes,
                                               Kappa,options$nu,options$W)
            z = list_update_part$z
            counts = list_update_part$counts
            Nclust = list_update_part$Nclust
            # cat('\n post z = ', z, '\n')
            # cat('\n post counts = ', counts, '\n')
            # cat('\n post Nclust = ', Nclust, '\n')
        }
        
        # Update alpha 
        if(options$UpdateAlpha){
            cat("Updating the a and d weights...")
            #alpha = Updatealpha_augmentation(alpha,Nclust,p,options$a_alpha,options$b_alpha)
            #list_update_alpha = Updatealpha_MH(alpha,Nclust,p,options$a_alpha,options$b_alpha, var_alpha_adp, iter, options$adaptiveAlpha)
            #alpha = list_update_alpha$alpha
            #var_alpha_adp = list_update_alpha$var_alpha_adp
        }
        
        # save results
        if(iter>nburn && (iter-nburn)%%thin == 0){
            it_saved = it_saved + 1 
            #save_res$Kappa[[it_saved]] = Kappa
            #save_res$Partition[it_saved,] = z
            #save_res$alpha[it_saved] = alpha
            #save_res$var_alpha_adp[it_saved] = var_alpha_adp
        }
        
        if(print){
            setTxtProgressBar(pb, iter)
        }
        
    }
    
    close(pb)
    return(save_res)
}
