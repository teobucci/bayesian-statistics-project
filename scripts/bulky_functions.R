#' Main function that updates the partition
#'
#' @param rho The partition in compact form (e.. rho=c(1,4,5) means that the first group has 1 element, the second has 4 elements, the last has 5 elements).
#' @param alpha_add Fixed probability of choosing an add move or delete move.
#' @param a_weights Vector of size (number of nodes - 1) containing at element j the weights to consider when ADDING a changepoint between point j and point j+1 (weights are non-normalized probabilities).
#' @param d_weights Vector of size (number of nodes - 1) containing at element j the weights to consider when DELETING a changepoint between point j and point j+1 (weights are non-normalized probabilities).
#' @param Theta_groups 
#' @param theta_prior 
#' @param sigma 
#'
#' @return
#' @export
#'
#' @examples
update_partition = function(rho,
                            alpha_add,
                            a_weights,
                            d_weights,
                            Theta_groups,
                            theta_prior,
                            sigma) {
    unifsample = runif(n = 1)
    choose_add = unifsample < alpha_add
    
    # OK
    proposal_list = proposal_ratio(rho, alpha_add, a_weights, d_weights, unifsample)
    log_proposal_ratioNow = log(proposal_list$ratio)
    candidate = proposal_list$candidate
    
    # compute proposed partition based on candidate and index of the group
    if (choose_add) {
        list_output_modify_partition = split_partition(candidate, rho)
    } else {
        list_output_modify_partition = merge_partition(candidate, rho)
    }
    proposed_rho = list_output_modify_partition$rho
    
    # TODO cambiare questo nome per renderlo coerente con ciÃ² che c'Ã¨ sotto
    THE_GROUP = list_output_modify_partition$group_index
    
    # OK
    # c'e' un ricalcolo inutile dell'indice del gruppo da splittare o mergiare
    log_priorRatioNow = log_priorRatio(theta_prior, sigma, current_rho, proposed_rho, choose_add)
    
    # OK
    # visto che il calcolo inutile dell'indice c'era sopra, l'ho messo pure qui
    log_likelihood_ratioNow = log_likelihood_ratio(choose_add, Theta_groups)
    
    alpha_accept <- min(1, exp(log_likelihood_ratioNow +
                               log_priorRatioNow +
                               log_proposal_ratioNow))
    
    if (runif(n = 1) < alpha_accept) {
        # accept the move
        
        if (choose_add) {
            # accepted move is a split
            # TODO split_partition(candidate,rho)
        }
        else{
            # accepted move is a merge
            # TODO merge_partition
        }
        partition_data
    } else {
        # don't do anything
        
    }
    
}



#' Adaptation step to update the weights vector a and d
#' The function takes as an input the current weights and updates them as a function of
#' the current iteration number t, the initial adaptation h, 
#' The function works in log scale
#' 
#' @param logweights vector of the logarithm of the current weights
#' @param alpha_target scalar indicating the target acceptance probability (optimal range around 0.10-0.15)
#' @param t number of the current iteration
#' @param h initial adaptation (must be >0)
#' @inheritParams update_partition
#'
#' @return Vector of updated logweights.
#' @export
#'
#' @examples
#' 
log_weights_adaptation = function(logweights, t, h, alpha_target, alpha_add) {
    if (!h > 0)
        stop("Adaptation step h must be positive")
    return (logweights + h * length(logweights) / t * (alpha_add - alpha_target))
}



#' Partition current data
#' Given y (vector of data) and rho (vector of the partition) the function splits the observations and
#' partitions them into the current groups 
#' partitions them into the partition rho
#' @param y - Vector of n ordered data
#' @inheritParams update_partition
#'
#' @return List where each element contains the corresponding group of y elements. If the dimensions of rho and y are not comparable, return an empty vector
#' @export
#'
#' @examples
partition_data <- function(y, rho) {
    if (sum(rho) != length(y))
        stop("The partition is not coherent with the data")
    
    dataPartition <- list()
    
    # number of groups
    M = length(rho)
    
    for (i in 1:M) {
        if (i == 1) {
            dataPartition[[i]] <- y[1:rho[i]]
            cumsum_rho = rho[i]
        } else {
            first_index = cumsum_rho + 1
            last_index = rho[i] + cumsum_rho
            dataPartition[[i]] <- y[first_index:last_index]
            cumsum_rho = cumsum_rho + rho[i]
        }
    }
    return(dataPartition)
}














#' Proposal Ratio
#'

#' @param choose_add Boolean to tell if we are performing an add move or not.
#'
#' @return The proposal ratio (not in log).
#' @export
#'
#' @examples
proposal_ratio = function(rho,
                         alpha_add,
                         a_weights,
                         d_weights,
                         choose_add) {
    # number of groups
    M = length(rho)
    
    n_elems = length(a_weights)
    
    # preliminary check for extreme cases
    if ((!choose_add & M == 1) | (choose_add & M == n_elems)) {
        # incompatible to delete when only one group is present
        # or to add when every point is a group
        return(list("ratio" = 0, "candidate" = -1)) # fictitious candidate
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
    
    # not all points can be selected for an add move
    # assign probability zero to those who cannot be
    a_weights_available = a_weights
    a_weights_available[cp_indexes] = 0
    a_weights_available_sum = sum(a_weights_available)
    
    # not all points can be selected for a delete move
    # assign probability zero to those who cannot be
    d_weights_available = d_weights
    d_weights_available[-cp_indexes] = 0
    d_weights_available_sum = sum(d_weights_available)
    
    if (choose_add) {
        draw_weights = a_weights_available
    } else {
        draw_weights = d_weights_available
    }
    
    # my candidate cannot be the first node by definition
    candidate = sample(2:(n_elem+1), 1, prob = draw_weights)
    
    if (choose_add & M == 1) {
        # case in which you choose to propose an add move
        # (with just 1 group) that may or may not be accepted
        ratio = (alpha_add / 1) * (a_weights_available_sum / a_weights[candidate])
        return(list("ratio" = ratio, "candidate" = candidate))
    }
    
    if (!choose_add & M == n_elems) {
        # case in which you choose to propose an delete move
        # (with every point being a group) that may or may not be accepted
        ratio = (alpha_add / 1) * (d_weights_available_sum / d_weights[candidate])
        return(list("ratio" = ratio, "candidate" = candidate))
    }
    
    
    # only the general cases remain
    if (choose_add) {
        ratio = (1 - alpha_add) / alpha_add *
            a_weights_available_sum / a_weights[candidate] *
            d_weights[candidate] / (d_weights[candidate] + d_weights_available_sum)
    } else {
        ratio = alpha_add / (1 - alpha_add) *
            d_weights_available_sum / d_weights[candidate] *
            a_weights[candidate] / (a_weights[candidate] + a_weights_available_sum)
    }
    
    return(list("ratio" = ratio, "candidate" = candidate))
}



#' Split Partition
#'
#' @param candidate_index Index of the the point where to split the group (equivalent to adding a changepoint).
#' @inheritParams update_partition
#' @return A list whose first element is the updated partition and the second is the index of the groups that has changed.
#' @export
#'
#' @examples
split_partition <- function(candidate_index, rho) {
    n_elems = sum(rho)
    # number of groups
    M = length(rho)
    new_rho = rep(NA, M + 1)
    
    # number of changepoints = n_elems - 1
    # case: all groups with 1 element or index out of bound
    if (n_elems == M | candidate_index > n_elems - 1) {
        return(list('rho' = rho, 'group_index' = -1))
    }
    
    cumsum_rho = cumsum(rho)
    found = F
    
    for (i in 1:M) {
        
        if (!found & cumsum_rho[i] == candidate_index) {
            # candidate_index is already a changepoint, return the original rho
            return(list('rho' = rho, 'group_index' = -1))
        }
        
        # update the partition in the general case
        # (either I have already split the group or not, just the index changes)
        if (!found) {
            new_rho[i] = rho[i]
        } else {
            new_rho[i + 1] = rho[i]
        }
        
        if (!found & cumsum_rho[i] > candidate_index) {
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



#' Merge Partition
#' 
#' @param candidate_index Index of the the point where to split the group (equivalent to adding a changepoint).
#' @inheritParams update_partition
#' @return A list whose first element is the updated partition and the second is the index of the groups that has changed.
#' @export
#'
#' @examples
merge_partition <- function(candidate_index, rho) {
    n_elems = sum(rho)
    # number of groups
    M = length(rho)
    new_rho = rep(NA, M - 1)
    
    # number of changepoints = n_elems - 1
    # case: only 1 group or index out of bound
    if (M == 1 | candidate_index > n_elems - 1) {
        return(list('rho' = rho, 'group_index' = -1))
    }
    
    cumsum_rho = cumsum(rho)
    found = F
    
    for (i in 1:(M - 1)) {
        if (!found & cumsum_rho[i] != candidate_index) {
            # candidate_index is already a changepoint, return the original rho
            return(list('rho' = rho, 'group_index' = -1))
        }
        
        # update the partition in the general case
        # (either I have already merged the group or not, just the index changes)
        if (!found) {
            new_rho[i] = rho[i]
        } else {
            new_rho[i] = rho[i + 1]
        }
        
        if (!found & cumsum[i] == candidate_index) {
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



#' Shuffle Partition
#' 
#' @param G Adjacency matrix of the Graph.
#' @inheritParams update_partition
#' @return The updated partition after a shuffle move, if accepted.
#' @export
#'
#' @examples
shuffle_partition <- function(rho, G) {
    # vedi Corradin p.16
    
    current_rho = rho
    
    # number of groups
    M = length(current_rho)
    
    if (M < 2) {
        # shuffling can be done only if the number of groups is at least 2
        return(current_rho)
    }
    
    # build new proposed rho
    proposed_rho = current_rho

    # K is the index of the shuffle group with K+1
    K <- sample(1:(M - 1), 1)
    # sample how many elements to keep in the K-th group
    l <- sample(1:(current_rho[K] + current_rho[K + 1] - 1), 1)
    # move the elements
    proposed_rho[K + 1] <- current_rho[K + 1] + current_rho[K] - l
    proposed_rho[K] <- l
    
    # compute log_prior_ratio
    log_prior_ratio = lpochhammer(1 - sigma, l)
                    + lpochhammer(1 - sigma, current_rho[j] + current_rho[j + 1] - l)
                    - lpochhammer(1 - sigma, current_rho[j] - 1)
                    - lpochhammer(1 - sigma, current_rho[j + 1] - 1)
                    + lfactorial(current_rho[j])
                    + lfactorial(current_rho[j + 1])
                    - lfactorial(l)
                    - lfactorial(current_rho[j] + current_rho[j + 1] - l)
    
    # to compute log_likelihood_ratio I need S

    S_current = get_S_from_G_rho(G, current_rho)
    S_star_current = get_S_star_from_S_and_rho(S, current_rho)
    
    S_proposed = get_S_from_G_rho(G, proposed_rho)
    S_star_proposed = get_S_star_from_S_and_rho(S, proposed_rho)
    
    # wrap the general rhoB into this version where I don't have to specify
    # alpha and beta (doing this for adaptiveness)
    rhoB = function(group1, group2, S, S_star) {
        return(rhoB_general(
            group1,
            group2,
            S,
            S_star,
            log = T,
            alpha = alpha,
            beta = beta
        ))
    }
    
    # here there's not the ratio with the Beta coefficients
    # compute log_likelihood_ratio
    log_likelihood_ratio = 0
    
    for (l in 1:(K - 1)) {
        # first numerator term
        log_likelihood_ratio = log_likelihood_ratio + rhoB(l, K, S_proposed, S_star_proposed)
        log_likelihood_ratio = log_likelihood_ratio + rhoB(l, K + 1, S_proposed, S_star_proposed)
        # first denominator term
        log_likelihood_ratio = log_likelihood_ratio + rhoB(l, K, S, S_star)
        log_likelihood_ratio = log_likelihood_ratio + rhoB(l, K + 1, S, S_star)
    }
    
    for (l in (K + 2):M) {
        # second numerator term
        log_likelihood_ratio = log_likelihood_ratio + rhoB(l, K, S_proposed, S_star_proposed)
        log_likelihood_ratio = log_likelihood_ratio + rhoB(l, K + 1, S_proposed, S_star_proposed)
        # second denominator term
        log_likelihood_ratio = log_likelihood_ratio + rhoB(l, K, S, S_star)
        log_likelihood_ratio = log_likelihood_ratio + rhoB(l, K + 1, S, S_star)
    }
    
    # third numerator term
    log_likelihood_ratio = log_likelihood_ratio + rhoB(K, K + 1, S_proposed, S_star_proposed)
    log_likelihood_ratio = log_likelihood_ratio + rhoB(K, K, S_proposed, S_star_proposed)
    log_likelihood_ratio = log_likelihood_ratio + rhoB(K + 1, K + 1, S_proposed, S_star_proposed)
    # third denominator term
    log_likelihood_ratio = log_likelihood_ratio + rhoB(K, K + 1, S, S_star)
    log_likelihood_ratio = log_likelihood_ratio + rhoB(K, K, S, S_star)
    log_likelihood_ratio = log_likelihood_ratio + rhoB(K + 1, K + 1, S, S_star)

    # compute alpha_shuffle
    alpha_shuffle = min(1, exp(log_likelihood_ratio + log_prior_ratio))
    
    if (runif(n = 1) < alpha_shuffle) {
        # accept the shuffle
        return(proposed_rho)
    } else {
        # reject the shuffle
        return(current_rho)
    }
    
}


#' Get Non-Edges from S and partition
#'
#' @param S MxM matrix, where M is number of groups, containing the sum of edges.
#' @inheritParams update_partition
#'
#' @return
#' @export
#'
#' @examples
get_S_star_from_S_and_rho = function(S, rho){
    # number of groups
    M = length(rho)
    
    # initialize S_star matrix
    S_star = matrix(numeric(M * M), nrow = M, byrow = T)
    
    # loop through the groups
    for (l in 1:M) {
        for (m in 1:l) {
            if (l == m){
                S_star[l,m] = rho[l] * (rho[l] - 1) / 2 - S[l,m]
            } else {
                S_star[l,m] = rho[l] * rho[m] - S[l,m]
                S_star[m,l] = S_star[l,m]
            }
        }
    }
    
    return(S_star)
}

# auxiliary function to evaluate the beta
rhoB_general = function(group1,
                group2,
                S,
                S_star,
                log = T,
                alpha,
                beta) {
    if (log) {
        return(lbeta(alpha + S[group1, group2], beta + S_star[group1, group2]))
    } else{
        return(beta(alpha + S[group1, group2], beta + S_star[group1, group2]))
    }
}

#' Log Likelihood Ratio
#'
#' @inheritParams update_partition
#' @inheritParams shuffle_partition
#' @param current_rho Current partition.
#' @param proposed_rho Proposed partition.
#' @param alpha Parameter of the Beta of the Graph.
#' @param beta Parameter of the Beta of the Graph.
#'
#' @return The likelihood ratio in log.
#' @export
#'
#' @examples
log_likelihood_ratio = function(alpha_add,
                               a_weights,
                               d_weights,
                               G,
                               current_rho,
                               proposed_rho,
                               choose_add,
                               alpha = 1,
                               beta = 1) {
    # differentiate delete/merge case
    if (!choose_add) {
        # swap rhos 'cause we're lazy
        temp = current_rho
        current_rho = proposed_rho
        proposed_rho = temp
    }
    
    # number of groups
    M = length(current_rho)
    
    # wrap the general rhoB into this version where I don't have to specify
    # alpha and beta (doing this for adaptiveness)
    rhoB = function(group1, group2, S, S_star) {
        return(rhoB_general(
            group1,
            group2,
            S,
            S_star,
            log = T,
            alpha = alpha,
            beta = beta
        ))
    }
        
    
    S_current = get_S_from_G_rho(G, current_rho)
    S_star_current = get_S_star_from_S_and_rho(S, current_rho)
    
    S_proposed = get_S_from_G_rho(G, proposed_rho)
    S_star_proposed = get_S_star_from_S_and_rho(S, proposed_rho)
    
    K = get_index_changed_group(current_rho, proposed_rho)
    
    log_ratio = -(M + 1) * rhoB(0, 0)
    # TODO mettere questo in tutti quelli sotto per
    # consentire l'aggiornamento di alpha e beta
    
    for (l in 1:(K - 1)) {
        # first numerator term
        log_ratio = log_ratio + rhoB(l, K, S_proposed, S_star_proposed)
        log_ratio = log_ratio + rhoB(l, K + 1, S_proposed, S_star_proposed)
        # first denominator term
        log_ratio = log_ratio - rhoB(l, K, S_current, S_star_current)
    }
    
    for (m in (K + 2):(M + 1)) {
        # second numerator term
        log_ratio = log_ratio + rhoB(K, m, S_proposed, S_star_proposed) +
                                rhoB(K + 1, m, S_proposed, S_star_proposed)
    }
    
    # third numerator term
    log_ratio = log_ratio + rhoB(K, K + 1, S_proposed, S_star_proposed) +
                            rhoB(K, K, S_proposed, S_star_proposed) +
                            rhoB(K + 1, K + 1, S_proposed, S_star_proposed)
    
    for (m in (K + 1):M) {
        # second denominator term
        log_ratio = log_ratio - rhoB(K, m, S_current, S_star_current)
    }
    
    # third denominator term
    log_ratio = log_ratio - rhoB(K, K, S_current, S_star_current)
    
    if (!choose_add) {
        # in the delete/merge case we have to invert everything
        log_ratio = -log_ratio
    }
    
    return(log_ratio)
    
}


get_index_changed_group = function(current_rho,proposed_rho){
    # indexes of the changepoints in the current partition
    cp_idxs_current = cumsum(current_rho)
    
    # move from the rho representation to r representation
    # i.e. from c(2,3) to c(0,1,0,0,1)
    # both for the current and the proposed partition
    
    current_r = rep(0, sum(current_rho))
    current_r[cp_idxs_current] = 1
    
    proposed_r = rep(0, sum(proposed_rho))
    proposed_r[cumsum(proposed_rho)] = 1
    
    # now by making the difference I can extract the index
    # of the new changepoint (or the one deleted)
    # there's no need for an absolute value because in the delete move
    # they have been fictitiously swapped to avoid repeating code
    tau = which.max(proposed_r - current_r)
    
    # take all the indexes of the cp smaller than tau
    temp = which(cp_idxs_current < tau)
    # take the last element of this list to have the
    # index of the group affected
    K = temp[length(temp)] + 1
    
    return(K)
}



log_priorRatio = function(theta_prior,
                          sigma,
                          current_rho,
                          proposed_rho,
                          choose_add)
{
    # differentiate delete/merge case
    if (!choose_add) {
        # swap rhos 'cause we're lazy
        temp = current_rho
        current_rho = proposed_rho
        proposed_rho = temp
    }
    
    # number of groups in the current partition
    M = length(current_rho)
    
    K = get_index_changed_group(current_rho,proposed_rho)
    
    # move from the current quantities to the math formulas ones
    
    # cardinality of the group in the proposed partition
    n_star_s = proposed_rho[K]
    # same but the next one
    n_star_s_plus_1 = proposed_rho[K + 1]
    # in the add move, the cardinality in the current partition is the sum
    n_s = n_star_s + n_star_s_plus_1
    
    # compute the prior ratio
    log_ratio = - log(M) + log(theta_prior + M * sigma)
                + lpochhammer(1 - sigma, n_star_s - 1)
                + lpochhammer(1 - sigma, n_star_s_plus_1 - 1)
                - lpochhammer(1 - sigma, n_s - 1)
                + lfactorial(n_s)
                - lfactorial(n_star_s)
                - lfactorial(n_star_s_plus_1)
    
    if (!choose_add) {
        # in the delete/merge case we have to invert everything
        log_ratio = -log_ratio
    }
    
    return(log_ratio)
}


#' Full-conditional for theta as in Martinez and Mena
#'
#' @param c first parameter of the shifted gamma prior 
#' @param d second parameter of the shifted gamma prior
#' @param candidate 
#' @param k 
#'
#' @return
#' @export
#'
#' @examples
#' #TODO  check that all the parameters make sense
full_conditional_theta <- function(c, d, candidate, k, n){
    weights <- rep(0,(k+1))
    z <- rbeta(1,candidate + 2, n)
    f <- rexp (1,candidate + 1)
    for(j in 1:k){ 
        # compute theta
        weight_j=compute_weights_theta(c, d, n, sigma, k, j-1, f, z)
        gamma_weights[j] <- weight_j
    }
    #Normalizing the weights
    gamma_weights = gamma_weights/sum(gamma_weights)
    
    #  I choose a random sample
    u = runif(1)
    
    #TODO understand the meaning of this step
    component <- min(which(cumsum(gamma_weights) > u))
    
    #BEWARE! There might be an error here on sigma, but it may depend on how it is passed
    theta <- shifted_gamma(prior_c+ (component-1), prior_d + f -log(z), sigma) #!!!shouldn't it be -sigma??
    
    return(theta)
}



#' Full-conditional for sigma as in Martinez and Mena
#'
#' The formula is on page 13 - ACHTUNG! They did everything in log and 
#' returned the logged result, but I am returning the unlogged version 
#' at the moment
#'
#' @param sigma ? not sure whether it is the previous
#' @param theta other parameter for the prior
#' @param k I think the changepoint index??? To be checked
#' @param rho current partition
#' @param a first parameter of the distribution
#' @param b second parameter of the distribution
#' @param c third parameter of the distribution
#' @param d fourth parameter of the distribution
#'
#' @return
#' @export
#'
#' @examples
full_conditional_sigma <- function(sigma,theta,k,rho,a,b,c,d){
    
    #First product term
    log_prod_1 = numeric(0)
    for(i in 1:(k-1)){
        log_prod_1 = log_prod_1 + log(theta + i*sigma)
    }
    
    #Second product term
    log_prod_2 <- numeric(0)
    for(i in 1:k){
        log_prod_2 = log_prod_2 + log_pochhammer((1-sigma),(rho[i]-1))
    }
    
    #Final output
    output <- (a-1) * log(sigma) + (b-1) * log(1-sigma) + 
        (c-1)*log(theta + sigma) + log(exp(-d*sigma)) + 
        log_prod_1 + log_prod_2
    
    return(exp(output))
}






set_options = function(sigma0,
                       theta_prior0,
                       rho0,
                       weights_a0,
                       weights_d0,
                       alpha_target,
                       mu_beta=0.5, # mu of the Beta
                       sig2_beta=0.2, # variance of the Beta
                       d=3,
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
    "mu_beta"            = mu_beta,
    "sig2_beta"          = sig2_beta,
    "d"                  = d,
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

#' Estimate alpha and beta from a Beta function given mean and variance
#' See https://stats.stackexchange.com/a/12239 for details.
#' 
#' @param mu Mean of the Beta.
#' @param var Variance of the Beta.
#'
#' @return List with alpha and beta of the corresponding Beta distribution.
#' @export
#'
#' @examples
estimate_Beta_params <- function(mu, var) {
    alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
    beta <- alpha * (1 / mu - 1)
    return(list(alpha = alpha, beta = beta))
}


Gibbs_sampler = function(data, niter, nburn, thin,
                         options,
                         seed=1234, print=T)
{
  n = nrow(data) # number of observations
  p = ncol(data) # number of nodes
  n_total_iter = nburn + niter*thin # total iterations to be made
  
  # dynamic parameters
  sigma            = options$sigma0 # initial parameter of the PY prior
  theta_prior      = options$theta_prior0 # initial parameter of the PY prior
  rho              = options$rho0 # initial partition (eg. c(150,151))
  weights_a        = options$weights_a0 # add weights
  weights_d        = options$weights_d0 # del weights
  
  # constant parameters
  alpha_add = options$alpha_add # probability of choosing add over delete
  alpha_target = options$alpha_target # target alpha for adapting weights
  
  
 
  d = options$d # parameter for the Wishart
  
  
  if(sum(rho) != p)
    stop("The partition rho must sum to p")
  if(d<3)
    stop("The Wishart's d must be greater or equal than 3")
  if(!(mu_beta > 0 & mu_beta < 1))
      stop("The mean of the Beta must be between 0 and 1")
  if(!(sig2_beta > 0 & sig2_beta < 0.25))
      stop("The mean of the Beta must be between 0 and 1")
  
  beta_params = estimate_Beta_params(mu_beta,sig2_beta)
  a = beta_params$alpha
  b = beta_params$beta 
  
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
  
  # initialize the sum of the weights of the graphs
  total_weights = 0
  
  #initialize the sum of all graphs
  total_graphs = matrix(0,p,p)
  
  # Start the simulation
  for(iter in 1:n_total_iter){
    
    # Update graph
    if(options$update_graph){
      # we ran a single iteration of BDgraph with iter=1 and burning=0
      #TODO think about extracting fixes parameters such as S, n, p which
      # at the moment are computed for every bdgraph iteration
      output = bdgraph( data, rho, n, method = "ggm", algorithm = "bdmcmc", iter=1,
                        burnin = 0, not.cont = NULL, g.prior = 0.5, df.prior = d,
                        CCG_D = NULL, g.start = "empty", jump = NULL, save = TRUE, print = 1000,
                        cores = NULL, threshold = 1e-8 )
      
      # Extracting the matrix with edges between groups
      Theta_groups = output$last_theta
      # Extracting the  precision matrix ?? Quale delle due ? effettivamente l'ultima precision matrix?
      # K_hat = output$K_hat
      last_K = output$last_K
      # Kappa = UpdatePrecision(options$nu,options$W,n,U,z)
      # Updating total_weights
      total_weights = total_weights + output$all_weights
      # Updating total_graphs taking into consideration the weights
      total_graphs = total_graphs + output$last_graph * output$all_weights
      
      
    }
    
    if(options$update_partition){
      # TODO
      update_partition(rho,alpha_add,a_weights,d_weights,Theta_groups,theta_prior,sigma)
    }
    
    
    if(options$update_weights){
      weights_a = exp(log_weights_adaptation(
        log(weights_a), iter, 1/p, alpha_target, alpha_add))
      
      weights_d = exp(log_weights_adaptation(
        log(weights_d), iter, 1/p, alpha_target, 1-alpha_add))
    }
    
    if(options$perform_shuffle){
      rho = shuffle_partition(rho)
    }
    
    if(options$update_sigma){
      # TODO mettere la full conditional di sigma qui
    }
    
    if(options$update_theta_prior){
      # TODO mettere la full conditional di theta prior qui
    }
    
    # save results only on thin iterations
    # (i.e. only save multiples of thin)
    if(iter > nburn & (iter - nburn)%%thin == 0) {
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
    # if(iter>nburn & (iter-nburn)%%thin == 0){
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
  graph_final = total_graphs/total_weights
  save_res$graph = graph_final
  
  close(pb)
  return(save_res)
}
