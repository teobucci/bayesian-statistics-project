#' Title
#'
#' @param rho 
#' @param alpha_add 
#' @param a_weights 
#' @param d_weights 
#' @param Theta_groups 
#' @param theta_prior 
#' @param sigma 
#'
#' @return
#' @export
#'
#' @examples
update_partition = function(rho,alpha_add,a_weights,d_weights,Theta_groups,theta_prior,sigma){
    
    unifsample = runif(n=1)
    choose_add = unifsample < alpha_add
    
    # OK
    proposal_list = proposalRatio(rho, alpha_add, a_weights, d_weights, unifsample)
    proposalRatioNow = log(proposal_list$ratio)
    candidate = proposal_list$candidate
    
    # compute proposed partition based on candidate and index of the group
    if(choose_add){
        t = splitPartition(candidate, rho)
    } else {
        t = mergePartition(candidate, rho)
    }
    proposed_rho = t$rho
    THE_GROUP = t$group_index # TODO cambiare questo nome per renderlo coerente con ciò che c'è sotto
    
    # OK qua dentro guarda che forse c'è un ricalcolo inutile dell'indice del gruppo da splittare o mergiare
    priorRatioNow = log_priorRatio(theta_prior, sigma, current_rho, proposed_rho, choose_add)
    
    # Status? 
    likelihoodRatioNow = log_likelihoodRatio(choose_add,Theta_groups)
    
    alpha_accept <- min(1, exp(likelihoodRatioNow + priorRatioNow + proposalRatioNow))
    
    if (runif(n=1) < alpha_accept){ # accept the move

        if(choose_add){ # accepted move is a split
            # TODO splitPartition(candidate,rho)
        }
        else{ # accepted move is a merge
            # TODO mergePartition
        }
        partitionData
    } else { # don't do anything
        
    }
    
}
    
    
    
    
    
    
    
    
    
    
    

#' Adaptation step to update the weights vector a and d
#' The function takes as an input the current weights and updates them as a function of
#' the current iteration number t, the initial adaptation h, 
#' The function works in log scale
#' 
#' @param logweights vector of the logarithm of the current weights
#' @param alpha_target scalar indicating the target acceptance probability (optimal range around 0.10-0.15)
#' @param alpha_add probability of adding a move (usually 0.5)
#' @param t number of the current iteration
#' @param h initial adaptation (must be >0)
#'
#' @return the vector of updated logweights
#' @export
#'
#' @examples
#' 
logAdaptation = function(logweights, t, h, alpha_target, alpha_add){ 
    if(!h>0)
        stop("Adaptation step h must be positive")
    return (logweights + h*length(logweights) / t * (alpha_add - alpha_target))
}


#' Partition current data
#' Given y (vector of data) and rho (vector of the partition) the function splits the observations and
#' partitions them into the current groups 
#' partitions them into the partition rho
#' @param y - vector of n ordered data
#' @param rho - partition written in compact form 
#' e.g. rho=c(2,4,5) means the first group has 2 elements, the second has 4 and the third has 5
#'
#' @return a list, where each element contains the corresponding group of y elements.
#' If the dimensions of rho and y are not comparable, return an empty vector
#' @export
#'
#' @examples
partitionData <- function(y,rho){ 
  if(sum(rho)==length(y)){ #checks that the dimensionality is okay
  
    dataPartition <- list()
    
    for (i in 1:length(rho)){
      if(i == 1) 
        {
          dataPartition[[i]] <- y[1:rho[i]]
          cumsum_rho=rho[i]
        }
      else
       {
         first_index=cumsum_rho + 1
         last_index=rho[i] + cumsum_rho
         dataPartition[[i]] <- y[first_index : last_index]
         cumsum_rho=cumsum_rho+rho[i]
       }
    }
    return(dataPartition) 
    }
}



#Create Theta_groups
#TODO Add checks
#' Title
#'
#' @param rho 
#' @param G 
#'
#' @return
#' @export
#'
#' @examples
create_Theta=function(rho, G){
    n_groups=length(rho)
    Theta= as.numeric(n_groups*n_groups)
    z=rho_to_z(rho)
    p = sum(rho)
    for(i in 1:p){
        for(j in 1:(i-1)){
            if (G[ i*p + j] == 1){
                ++Theta[z[i] * n_groups + z[j]];
                ++Theta[z[j] * n_groups + z[i]];
            }
        }
    }
    return(Theta)
}







#Rho to changepoint function
#add check
#' Title
#'
#' @param rho 
#'
#' @return
#' @export
#'
#' @examples
rho_to_r=function(rho){
    cumsum_rho=cumsum(rho)
    total_n=sum(rho)
    r<-numeric(total_n)
    for(i in cumsum_rho){
        r[i]=1
    }
    return(r)
}

#Rho to z function
#Add check
#' Title
#'
#' @param rho 
#'
#' @return
#' @export
#'
#' @examples
rho_to_z=function(rho){
    total_n=sum(rho)
    n_groups<-length(rho)
    z<-{}
    for(i in 1:n_groups){
        z<-c(z, rep(i,rho[i]))
    }
    return(z)
}


#' Computes the proposal ratio from the file "Samplingstrategy_nonparam" at page 4
#' Follows strictly the paper, just adds the possibility of having alpha_add set by the author
#' NOTE: I have not added the checks for the range of the parameters ranging from 0 and 1. Anyway, they are easy to include in the first if
#'
#' @param rho the partition in compact form (e.. rho=c(1,4,5) means that the first group has 1 element, the second has 4 elements, the last has 5 elements)
#' @param alpha_add fixed probability of choosing ADD or DELETE as a move
#' @param a_weights vector of size n-1 (number of datapoints - 1) containing at element j the weights to consider when ADDING a changepoint between point j and point j+1 (weights are non-normalized probabilities)
#' @param d_weights vector of size n-1 (number of datapoints - 1) containing at element j the weights to consider when DELETING a changepoint between point j and point j+1 (weights are non-normalized probabilities)
#'
#' @return the proposal ratio, which is 
#' @export
#'
#' @examples
proposalRatio=function(rho, alpha_add, a_weights, d_weights, unifsample){
    
    #alpha_add=0.5
    num_groups = length(rho)
    choose_add = unifsample < alpha_add
    n_elems = length(a_weights)
    
    # preliminary check for extreme cases
    if((!choose_add && num_groups==1) || (choose_add && num_groups==n_elems) ){
        # incompatible to delete when only one group is present
        # or to add when every point is a group
        return(list("ratio"=0,"candidate"=-1)) # fictitious candidate
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
        ratio = alpha_add/1 * a_weights_available_sum / a_weights[candidate]
        return(list("ratio"=ratio,"candidate"=candidate))
    }
    
    if (!choose_add && num_groups==n_elems){
        # case in which you choose to propose an delete move
        # (with every point being a group) that may or may not be accepted
        ratio = alpha_add/1 * d_weights_available_sum / d_weights[candidate]
        return(list("ratio"=ratio,"candidate"=candidate))
    }
    
    
    # only the general cases remain
    if (choose_add){
        ratio = (1-alpha_add)/alpha_add * 
            a_weights_available_sum / a_weights[candidate] * 
            d_weights[candidate] / (d_weights[candidate] + d_weights_available_sum)
    } else {
        ratio = alpha_add/ (1-alpha_add) * 
            d_weights_available_sum / d_weights[candidate] *
            a_weights[candidate] / (a_weights[candidate] + a_weights_available_sum)
    }
    
    return(list("ratio"=ratio,"candidate"=candidate))
}



#' splitPartition in the compact form
#' 
#' NOTE: maybe it would be better to switch the representation, add a changepoint and switch back... 
#' let's think about it
#'
#' @param k index of the the point where to split the group (equivalent to adding a changepoint)
#' @param rho partition in compact form e.g., rho=c(2,3,4) means the first group has 2 elements, the second has three and the third has four
#'
#' @return a list whose first element is the updated partition 
#' and the second is the index of the groups that has changed (do not know if it is necessary, though)
#' @export
#'
#' @examples
splitPartition <- function(candidate_index, rho) {
    n_elems = sum(rho)
    n_groups = length(rho)
    new_rho = rep(NA,n_groups+1)
    
    # number of changepoints = n_elems - 1
    # case: all groups with 1 element or index out of bound
    if (n_elems == (length(rho)) || candidate_index > n_elems - 1) {
        return(list('rho' = rho, 'group_index' = -1))
    }
    
    cumsum_rho = cumsum(rho)
    found = F
    
    for (i in 1:n_groups) {
        if (!found && cumsum_rho[i] == candidate_index) {
            # candidate_index is already a changepoint, return the original rho
            return(list('rho' = rho, 'group_index' = -1))
        }
        # update the partition in the general case
        # (either I have already split the group or not, just the indexing changes)
        if (!found) {
            new_rho[i] = rho[i]
        } else {
            new_rho[i + 1] = rho[i]
        }
        if (!found && cumsum_rho[i] > candidate_index) {
            # just passed the element index - I am in the group to be split
            
            # index of the element minus the cumulative
            # umber of elements in the previous groups only if i!=1
            new_rho[i] = candidate_index - (i != 1) * cumsum_rho[i - 1 * (i != 1)]
            
            # dimension of the original group minus the elements moved to new_rho[i]
            new_rho[i + 1] = rho[i] - new_rho[i]
            
            # save the index of the group that has changed
            j = i
            
            found = T
        }
    }
    return(list('rho' = new_rho, 'group_index' = j))
}

# version that preallocates z    
rho_to_z = function(rho){
    cumsum_rho = cumsum(rho)
    z = numeric(sum(rho))
    n_groups = length(rho)
    for(i in 1:n_groups){
        if(i == 1)
            z[1:cumsum_rho[1]] = i
        else
            z[(cumsum_rho[i-1]+1):cumsum_rho[i]] = i
    }
}


#' mergePartition in the compact form
#' 
#' NOTE: maybe it would be better to switch the representation, add a changepoint and switch back... 
#' let's think about it
#' 
#' @param k index of the the point where to split the group (equivalent to adding a changepoint)
#' @param rho partition in compact form e.g., rho=c(2,3,4) means the first group has 2 elements, the second has three and the third has four
#'
#' @return a list whose first element is the updated partition 
#' and the second is the index of the groups that has changed (do not know if it is necessary, though)
#' @export
#'
#' @examples
mergePartition <- function(candidate_index, rho) {
    n_elems = sum(rho)
    n_groups = length(rho)
    new_rho = rep(NA, n_groups - 1)
    
    # number of changepoints = n_elems - 1
    # case: only 1 group or index out of bound
    if ((length(rho) == 1) || candidate_index > n_elems - 1) {
        return(list('rho' = rho, 'group_index' = -1))
    }
    
    cumsum_rho = cumsum(rho)
    found = F
    
    for (i in 1:(n_groups - 1)) {
        if (!found && cumsum_rho[i] != candidate_index) {
            # candidate_index is already a changepoint, return the original rho
            return(list('rho' = rho, 'group_index' = -1))
        }
        # update the partition in the general case
        # (either I have already merged the group or not, just the indexing changes)
        if (!found) {
            new_rho[i] = rho[i]
        } else {
            new_rho[i] = rho[i + 1]
        }
        if (!found && cumsum[i] == candidate_index) {
            # I am at the changepoint between the two groups to be merged
            
            # index of the element minus the cumulative
            # number of elements in the previous groups
            new_rho[i] = rho[i] + rho[i + 1]
            
            # save the index of the group that has changed
            j = i
            
            found = T
        }
    }
    return(list('rho' = new_rho, 'group_index' = j))
}

#' shufflePartition (QUESTA RESTA MOLTO SIMILE A CORRADIN, NO?)
#' 
#' @param rho partition in compact form e.g., rho=c(2,3,4) means the first group has 2 elements, the second has three and the third has four
#'
#' @return the updated partition
#' @export
#'
#' @examples
shuffle <- function(rho){ # vedi Corradin p.16

    k = length(rho)
    
    if(k < 2){ # shuffling can be done only if the number of groups is at least 2
        return(rho)
    }
    
    new_rho = rho
        
    j <- sample(1:(k - 1),1)
    l <- sample(1:(rho[j]+rho[j+1] - 1),1)
    
    new_rho[j+1] <- rho[j+1] + rho[j] - l
    new_rho[j] <- l
    
    
    # compute alpha_shuffle
    
    prior_ratio = lpochhammer(1 - sigma, l)
                + lpochhammer(1 - sigma, rho[j] + rho[j+1] - l)
                - lpochhammer(1 - sigma, rho[j] - 1)
                - lpochhammer(1 - sigma, rho[j+1] - 1)
                + lfactorial(rho[j])
                + lfactorial(rho[j+1])
                - lfactorial(l)
                - lfactorial(rho[j] + rho[j+1] - l)
    
    likelihood_ratio = 
        # TODO bisogna scrivere S e avere i due casi tenendo conto
        # di quanti elementi vanno da una parte a un'altra
    
    
    
    alpha_shuffle = min(1, likelihood_ratio * prior_ratio)
    
    if(runif(n=1) < alpha_shuffle){
        # accept the shuffle
        return(new_rho)
    } else {
        # reject the shuffle
        return(rho)
    }
    
}

#Useful Functions

#' Computes the rising factorial (also called Pochhammer symbol), eventually in logarithmic form.
#' 
#' @param x the value for which the rising factorial is to be calculated
#' @param n the power x is to be raised to
#' @param log if true, the rising factorial is calculated in logarithmic form
#'
#' @return the rising factorial of x to the n power
#' @export
#'
#' @examples
lpochhammer <- function(x, n, log = T) {
    if (n < 0)
        stop("The pochhammer operator doesn't allow n < 0")
    
    num_vec <- numeric(n)
    
    if (n != 0) {
        num_vec[1] = x
        for (i in 1:(n - 1)) {
            num_vec[i + 1] = (x + i)
        }
    }
    
    if (n == 0) {
        num_vec[1] = 1
    }
    
    if (log)
        return(sum(log(num_vec)))
    else
        return(prod(num_vec))
}

# ---------------------------------------------
# TEST
# ---------------------------------------------
    n_gruppi = 3
    eStoyAqui = c(
        5,1,0,
        1,6,3,
        0,3,1
    )
    eStoyAqui = matrix(eStoyAqui,nrow=n_gruppi,ncol=n_gruppi,byrow=T)
    eStoyAqui
    
    i = 3
    j = 3

    if(i!=j){
        S_ij = eStoyAqui[i,j] # number of edges between cluster i and cluster j
        S_star_ij = rho[i] * rho[j] - S_ij
    } else {
        if(rho[i] != 1){
            S_ij = eStoyAqui[i,j]
            S_star_ij = rho[i] * (rho[i]-1) / 2 - S_ij
        } else { # no internal edges
            S_ij = 0
            S_star_ij = 0
        }
    }
    
# ---------------------------------------------

log_likelihoodRatio = function(rho, alpha_add, a_weights, d_weights,Theta_groups){

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
    
    if("il nodo non è agli estremi"){ # questo check non serve (semicit. Corradin)
        
        # ipotizzo la situa in cui aggiungo un cp NEL MEZZO
        C = c(20,20,5)
        C_star = c(20,7,13,5)
        
        M = length(C)
        S = 2 # è quello che è stato splittato
        # due nuovi sono in S e S+1 (2 e 3)
        M, M+1
        
        # bisogna gestire S_lm e S^star_lm
        
        ratio = -(M+1)*rhoB(0,0)
        
        for(l in 1:(S-1)){
            # first numerator term
            ratio = ratio + rhoB(C_l_star,C_S_star) + rhoB(C_l_star,C_S+1_star)
            # first denominator term
            ratio = ratio - rhoB(C_l,C_S)
        }
        
        for(m in (S+2):(M+1)){
            # second numerator term
            ratio = ratio + rhoB(C_S_star,C_m_star) + rhoB(C_S+1_star,C_m_star)
        }
        
        # third numerator term
        ratio = ratio + rhoB(C_S_star,C_S+1_star) + rhoB(C_S_star,C_S_star) + rhoB(C_S+1_star,C_S+1_star)
        
        for(m in (S+1):M){
            # second denominator term
            ratio = ratio - rhoB(C_S,C_m)
        }
        
        # third denominator term
        ratio = ratio - rhoB(C_S,C_S)
        
        return(exp(ratio))
    }

}



log_priorRatio = function(theta_prior, sigma, current_rho, proposed_rho, choose_add){
    
    if(!choose_add){ # delete/merge case
        # swap rhos 'cause we're lazy
        temp = current_rho
        current_rho = proposed_rho
        proposed_rho = temp
    }
    
    M = length(current_rho)
    
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
    
    ratio = - log(M) + (theta_prior + M*sigma) + lpochhammer(1-sigma, n_star_s - 1)
                                         + lpochhammer(1-sigma, n_star_s_plus_1 - 1)
                                         - lpochhammer(1-sigma, n_s - 1)
                                         + lfactorial(n_s)
                                         - lfactorial(n_star_s)
                                         - lfactorial(n_star_s_plus_1)
    
    if(!choose_add){ # in the delete/merge case we have to invert everything
        ratio = -ratio
    }
    
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


set_options = function(sigma0,
                       theta_prior0,
                       rho0,
                       weights_a0,
                       weights_d0,
                       alpha_target,
                       a,
                       b,
                       alpha_add=0.5,
                       update_sigma=T,
                       update_theta_prior=T,
                       update_weights=T,
                       update_partition=T,
                       update_graph=T,
                       perform_shuffle=T){
    
    options = list(
        "sigma0"             = sigma0,
        "theta_prior0"       = theta_prior0,
        "rho0"               = rho0,
        "weights_a0"         = weights_a0,
        "weights_d0"         = weights_d0,
        "alpha_target"       = alpha_target,
        "a"                  = a,
        "b"                  = b,
        "alpha_add"          = alpha_add,
        "update_sigma"       = update_sigma,
        "update_theta_prior" = update_theta_prior,
        "update_weights"     = update_weights,
        "update_partition"   = update_partition,
        "update_graph"       = update_graph,
        "perform_shuffle"    = perform_shuffle,
        )
    return(options)
}

Gibbs_sampler = function(data, niter, nburn, thin,
                         options,
                         seed=1234, print=T)
{
    n = nrow(data) # number of observations
    p = ncol(data) # number of nodes
    n_total_iter = nburn + niter*thin # total iterations to be made
    
    # dynamic parameters
    sigma            = options$sigma0 # initial parameter of the nonparametric prior
    theta_prior      = options$theta_prior0 # initial parameter of the nonparametric prior
    rho              = options$rho0 # initial partition (eg. c(150,151))
    weights_a        = options$weights_a0 # add weights
    weights_d        = options$weights_d0 # del weights
           
    # constant parameters
    alpha_add = options$alpha_add # probability of choosing add over delete
    alpha_target = options$alpha_target # target alpha for adapting weights
    
    a = options$a # parameter for the likelihood of the graph (Beta(a,b))
    b = options$b # parameter for the likelihood of the graph (Beta(a,b))

    
    if(sum(rho) != p)
        stop("The partition rho must sum to p")
    # TODO check qui dei parametri del grafo
    # if(nrow(options$W)!=p)
    #     stop("nrow W not coherent with ncol(data)")  
    # if(nrow(options$Kappa0)!=p)
    #     stop("nrow Kappa0 not coherent with ncol(data)")
    
    
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
    
    #Start the simulation
    for(iter in 1:n_total_iter){
        
        # Update graph
        if(options$update_graph){
            #we ran a single iteration of BDgraph with iter=1 and burning=0
            #TODO think about extracting fixes parameters such as S, n, p which
            # at the moment are computed for every bdgraph iteration
             output = bdgraph( data, rho, n, method = "ggm", algorithm = "bdmcmc", iter=1,
                               burnin = 0, not.cont = NULL, g.prior = 0.5, df.prior = 3,
                               CCG_D = NULL, g.start = "empty", jump = NULL, save = TRUE, print = 1000,
                               cores = NULL, threshold = 1e-8 )
             
             #Extracting the matrix with edges between groups
             Theta_groups = output$last_theta 
            # Kappa = UpdatePrecision(options$nu,options$W,n,U,z)
        
        }
        
        if(options$update_partition){
            # TODO
            update_partition(rho,alpha_add,a_weights,d_weights,Theta_groups,theta_prior,sigma)
        }
        
        
        if(options$update_weights){
            weights_a = exp(logAdaptation(
                log(weights_a), iter, 1/p, alpha_target, alpha_add))
            
            weights_d = exp(logAdaptation(
                log(weights_d), iter, 1/p, alpha_target, 1-alpha_add))
        }
        
        if(options$perform_shuffle){
            rho = shuffle(rho)
        }
        
        if(options$update_sigma){
            # TODO mettere la full conditional di sigma qui
        }
        
        if(options$update_theta_prior){
            # TODO mettere la full conditional di theta prior qui
        }
        
        # save results only on thin iterations
        if(iter>nburn && (iter - nburn)%%thin == 0) {
            # TODO
            it_saved = it_saved + 1
        }
        
        if(print){
            setTxtProgressBar(pb, iter)
        }
        
        
        
            
        # Update partition
        # if(options$UpdatePartition){
        #     cat("Updating the partition...")
        #     # cat('\n pre z = ', z, '\n')
        #     # cat('\n pre counts = ', counts, '\n')
        #     # cat('\n pre Nclust = ', Nclust, '\n')
        #     list_update_part = UpdatePartition(z,counts,Nclust,alpha,
        #                                        MARGINAL = marginal_Wishartdes,
        #                                        Kappa,options$nu,options$W)
        #     z = list_update_part$z
        #     counts = list_update_part$counts
        #     Nclust = list_update_part$Nclust
        #     # cat('\n post z = ', z, '\n')
        #     # cat('\n post counts = ', counts, '\n')
        #     # cat('\n post Nclust = ', Nclust, '\n')
        # }
        
        # Update alpha 
        #if(options$UpdateAlpha){
        #    cat("Updating the a and d weights...")
            #alpha = Updatealpha_augmentation(alpha,Nclust,p,options$a_alpha,options$b_alpha)
            #list_update_alpha = Updatealpha_MH(alpha,Nclust,p,options$a_alpha,options$b_alpha, var_alpha_adp, iter, options$adaptiveAlpha)
            #alpha = list_update_alpha$alpha
            #var_alpha_adp = list_update_alpha$var_alpha_adp
        #}
        
        # save results
        # if(iter>nburn && (iter-nburn)%%thin == 0){
        #     it_saved = it_saved + 1 
        #     #save_res$Kappa[[it_saved]] = Kappa
        #     #save_res$Partition[it_saved,] = z
        #     #save_res$alpha[it_saved] = alpha
        #     #save_res$var_alpha_adp[it_saved] = var_alpha_adp
        # }
        
        # if(print){
        #     setTxtProgressBar(pb, iter)
        # }
        
    }
    
    close(pb)
    return(save_res)
}
