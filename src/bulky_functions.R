library(logr)

#' Main function that updates the partition
#'
#' @param rho The partition in compact form (e.g. rho=c(1,4,5) means that the first group has 1 element, the second has 4 elements and the last has 5 elements).
#' @param alpha_add Fixed probability of choosing an add move or delete move.
#' @param weights_a Vector of size (number of nodes - 1) containing at element j the weights to consider when ADDING a changepoint between point j and point j+1 (weights are non-normalized probabilities).
#' @param weights_d Vector of size (number of nodes - 1) containing at element j the weights to consider when DELETING a changepoint between point j and point j+1 (weights are non-normalized probabilities).
#' @param theta_prior Prior parameter as in Martinez and Mena (2014)
#' @param sigma_prior Prior parameter as in Martinez and Mena (2014)
#' @param G Adjacency matrix of the graph
#' @param beta_params Parameters of the Beta
#'
#' @return a new partition in the compact form
#' @export
#'
#' @examples
update_partition = function(rho_current,
                            alpha_add,
                            weights_a,
                            weights_d,
                            theta_prior,
                            sigma_prior,
                            G,
                            beta_params) {
    unifsample = runif(n = 1)
    choose_add = unifsample < alpha_add
    
    # number of groups
    M = length(rho_current)
    p = sum(rho_current)
    
    # force opposite choice if merging/splitting is not feasible
    if ((!choose_add && M == 1) || (choose_add && M == p)) {
        choose_add = !choose_add
    }

    if (choose_add){
        #log_print("Chosen move: ADD/SPLIT", console = FALSE)
    }else{
        #log_print("Chosen move: DELETE/MERGE", console = FALSE)
    }
    
    proposal_list = proposal_ratio(rho_current, alpha_add, weights_a, weights_d, choose_add)
    log_proposal_ratioNow = log(proposal_list$ratio)
    candidate = proposal_list$candidate
    
    # compute proposed partition based on candidate and index of the group
    if (choose_add) {
        list_output_modify_partition = split_partition(candidate, rho_current)
    } else {
        list_output_modify_partition = merge_partition(candidate, rho_current)
    }
    rho_proposed        = list_output_modify_partition$new_rho
    changed_group_index = list_output_modify_partition$changed_group_index
    
    #log_print("rho_current", console = FALSE)
    #log_print(rho_current, console = FALSE)
    #log_print("rho_proposed", console = FALSE)
    #log_print(rho_proposed, console = FALSE)

    log_prior_ratioNow = log_prior_ratio(
        theta_prior,
        sigma_prior,
        rho_current,
        rho_proposed,
        choose_add,
        changed_group_index
    )
    
    log_likelihood_ratioNow = log_likelihood_ratio(
        alpha_add,
        weights_a,
        weights_d,
        G,
        rho_current,
        rho_proposed,
        choose_add,
        beta_params$alpha,
        beta_params$beta,
        changed_group_index
    )
    
    alpha_accept <- min(1, exp(log_likelihood_ratioNow +
                               log_prior_ratioNow +
                               log_proposal_ratioNow))

    if (runif(n = 1) < alpha_accept) {
        accepted = TRUE
        #log_print("Move ACCEPTED", console = FALSE)
        rho_updated = rho_proposed
    } else {
        accepted = FALSE
        #log_print("Move REJECTED", console = FALSE)
        rho_updated = rho_current
    }
    return(
        list(
            "rho_updated" = rho_updated,
            "accepted" = as.numeric(accepted),
            "choose_add" = choose_add,
            "candidate" = candidate
        )
    )
    
}



#' Adaptation step to update the weights vectors a and d
#'
#' The function takes as an input the current weights and updates them as a function of
#' the current iteration number t, the initial adaptation h
#' The function works in log scale
#'
#' @param logweights Vector of the logarithm of the current weights
#' @param alpha_target Scalar indicating the target acceptance probability (optimal range empirically observed around 0.10-0.15)
#' @param t Number of the current iteration
#' @param h Initial adaptation (must be >0)
#' @inheritParams update_partition
#'
#' @return Vector of updated logweights
#' @export
#'
#' @examples
#'
update_weight = function(weights, index, h, t_over_p, alpha_add, alpha_target) {
    if (!h > 0)
        stop("Adaptation step h must be positive")
    # select the weight that has to be updated
    weight = weights[index]
    # update according to Benson
    weight = exp(log(weight) + h / t_over_p * (alpha_add - alpha_target))
    # put it back
    weights[index] = weight
    # return the vector of weights
    return(weights)
}



#' Partition current data
#' Given y (vector of data) and rho (vector of the partition), the function splits the observations and
#' partitions them into the current groups
#' partitions them into the partition rho
#' @param y Vector of n ordered data
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
                          weights_a,
                          weights_d,
                          choose_add) {
    # number of groups
    M = length(rho)
    
    n_elems = length(weights_a)
    
    # indexes of the changepoints
    cp_indexes <- get_group_indexes(rho)
    
    # exclude the last one because it's technically always 1 (a changepoint)
    cp_indexes <- cp_indexes[-length(cp_indexes)]
    
    # not all points can be selected for an add move
    # assign probability zero to those who cannot be
    weights_a_available = weights_a
    weights_a_available[cp_indexes] = 0
    weights_a_available_sum = sum(weights_a_available)
    
    # not all points can be selected for a delete move
    # assign probability zero to those who cannot be
    weights_d_available = weights_d
    weights_d_available[-cp_indexes] = 0
    weights_d_available_sum = sum(weights_d_available)
    
    if (choose_add) {
        draw_weights = weights_a_available
    } else {
        draw_weights = weights_d_available
    }
    
    # draw the candidate among the first 1:(p-1)
    candidate = sample(1:n_elems, 1, prob = draw_weights)
    
    if (choose_add && M == 1) {
        # case in which you choose to propose an add move
        # (with just 1 group) that may or may not be accepted
        ratio = (alpha_add / 1) * (weights_a_available_sum / weights_a[candidate])
        return(list("ratio" = ratio, "candidate" = candidate))
    }
    
    if (!choose_add && M == n_elems) {
        # case in which you choose to propose an delete move
        # (with every point being a group) that may or may not be accepted
        ratio = (alpha_add / 1) * (weights_d_available_sum / weights_d[candidate])
        return(list("ratio" = ratio, "candidate" = candidate))
    }
    
    # only the general cases remain
    if (choose_add) {
        ratio = (1 - alpha_add) / alpha_add *
            weights_a_available_sum / weights_a[candidate] *
            weights_d[candidate] / (weights_d[candidate] + weights_d_available_sum)
    } else {
        ratio = alpha_add / (1 - alpha_add) *
            weights_d_available_sum / weights_d[candidate] *
            weights_a[candidate] / (weights_a[candidate] + weights_a_available_sum)
    }
    
    return(list("ratio" = ratio, "candidate" = candidate))
}



#' Split Partition
#'
#' @param candidate_index Index of the the point where to split the group (equivalent to adding a changepoint).
#' @inheritParams update_partition
#' @return A list whose first element is the updated partition and the second is the index of the group that has changed.
#' @export
#'
#' @examples
split_partition <- function(candidate_index, rho) {

    # number of groups
    M = length(rho)
    new_rho = rep(NA, M + 1)
    
    group_indexes = get_group_indexes(rho)
    found = FALSE
    
    for (i in 1:M) {
        
        # update the partition in the general case
        # (either I have already split the group or not, just the index changes)
        if (!found) {
            new_rho[i] = rho[i]
        } else {
            new_rho[i + 1] = rho[i]
        }
        
        if (!found && group_indexes[i] > candidate_index) {
            # just passed the element index - I am in the group to be split
            
            # index of the element minus the cumulative
            # number of elements in the previous groups only if i!=1
            new_rho[i] = candidate_index - (i != 1) * group_indexes[i - 1 * (i != 1)]
            
            # dimension of the original group minus the elements moved to new_rho[i]
            new_rho[i + 1] = rho[i] - new_rho[i]
            
            # save the index of the group that has changed
            j = i
            
            found = TRUE
        }
    }
    return(list("new_rho" = new_rho, "changed_group_index" = j))
}



#' Merge Partition
#'
#' @param candidate_index Index of the the point where to split the group (equivalent to adding a changepoint).
#' @inheritParams update_partition
#' @return A list whose first element is the updated partition and the second is the index of the group that has changed.
#' @export
#'
#' @examples
merge_partition <- function(candidate_index, rho) {

    # number of groups
    M = length(rho)
    new_rho = rep(NA, M - 1)
    
    group_indexes = get_group_indexes(rho)
    found = FALSE
    
    for (i in 1:(M - 1)) {
        
        # update the partition in the general case
        # either I have already merged the group or not, just the index changes
        if (!found) {
            new_rho[i] = rho[i]
        } else {
            new_rho[i] = rho[i + 1]
        }
        
        if (!found && group_indexes[i] == candidate_index) {
            # I am at the changepoint between the two groups to be merged
            
            # index of the element minus the cumulative
            # number of elements in the previous groups
            new_rho[i] = rho[i] + rho[i + 1]
            
            # save the index of the group that has changed
            j = i
            
            found = TRUE
        }
    }
    return(list("new_rho" = new_rho, "changed_group_index" = j))
}



#' Shuffle Partition
#'
#' @param G Adjacency matrix of the Graph.
#' @inheritParams update_partition
#' @return The updated partition after a shuffle move, if accepted.
#' @export
#'
#' @examples
shuffle_partition <- function(rho_current, G, sigma_prior, alpha, beta) {
    # vedi Corradin p.16 step (ii)
    
    # number of groups
    M = length(rho_current)
    
    # shuffling can be done only if the number of groups is at least 2
    if (M < 2) return(rho_current)
    
    # build new proposed rho
    rho_proposed = rho_current

    # going to shuffle group K with group K+1
    K <- sample(1:(M - 1), 1)
    
    if (rho_current[K] == 1 && rho_current[K+1] == 1){
        #log_print("SHUFFLE: cannot shuffle anything without reproposing the same partition", console = FALSE)
        return(rho_current)
    }
    
    # sample how many elements to keep in the K-th group
    sample_vector <- 1:(rho_current[K] + rho_current[K + 1] - 1)
    # avoid proposing the same partition (shuffling nothing)
    sample_vector <- sample_vector[-c(rho_current[K])]
    l <- sample(sample_vector, 1)

    # move the elements
    rho_proposed[K + 1] <- rho_current[K + 1] + rho_current[K] - l
    rho_proposed[K] <- l
    
    # compute log_prior_ratio
    log_prior_ratio = lpochhammer(1 - sigma_prior, rho_proposed[K])
                    + lpochhammer(1 - sigma_prior, rho_proposed[K + 1])
                    - lpochhammer(1 - sigma_prior, rho_current[K] - 1)
                    - lpochhammer(1 - sigma_prior, rho_current[K + 1] - 1)
                    + lfactorial(rho_current[K])
                    + lfactorial(rho_current[K + 1])
                    - lfactorial(rho_proposed[K])
                    - lfactorial(rho_proposed[K + 1])
    
    # to compute log_likelihood_ratio I need S

    S_current = get_S_from_G_rho(G, rho_current)
    S_star_current = get_S_star_from_S_and_rho(S_current, rho_current)
    
    S_proposed = get_S_from_G_rho(G, rho_proposed)
    S_star_proposed = get_S_star_from_S_and_rho(S_proposed, rho_proposed)
    
    # wrap the general fB into this version where I don't have to specify
    # alpha and beta (doing this for adaptiveness)
    fB = function(group1, group2, S, S_star) {
        return(fB_general(
            group1,
            group2,
            S,
            S_star,
            alpha = alpha,
            beta = beta,
            log = TRUE
        ))
    }
    
    # here there's not the ratio with the Beta coefficients
    # compute log_likelihood_ratio
    log_likelihood_ratio = 0
    
    for (l in 1:(K - 1)) {
        if (l > (K - 1)) break; # needed because R for-loops suck
        # first numerator term
        log_likelihood_ratio = log_likelihood_ratio + fB(l, K, S_proposed, S_star_proposed)
        log_likelihood_ratio = log_likelihood_ratio + fB(l, K + 1, S_proposed, S_star_proposed)
        # first denominator term
        log_likelihood_ratio = log_likelihood_ratio + fB(l, K, S_current, S_star_current)
        log_likelihood_ratio = log_likelihood_ratio + fB(l, K + 1, S_current, S_star_current)
    }
    
    for (l in (K + 2):M) {
        if (l > M) break; # needed because R for-loops suck
        # second numerator term
        log_likelihood_ratio = log_likelihood_ratio + fB(l, K, S_proposed, S_star_proposed)
        log_likelihood_ratio = log_likelihood_ratio + fB(l, K + 1, S_proposed, S_star_proposed)
        # second denominator term
        log_likelihood_ratio = log_likelihood_ratio + fB(l, K, S_current, S_star_current)
        log_likelihood_ratio = log_likelihood_ratio + fB(l, K + 1, S_current, S_star_current)
    }
    
    # third numerator term
    log_likelihood_ratio = log_likelihood_ratio + fB(K, K + 1, S_proposed, S_star_proposed)
    log_likelihood_ratio = log_likelihood_ratio + fB(K, K, S_proposed, S_star_proposed)
    log_likelihood_ratio = log_likelihood_ratio + fB(K + 1, K + 1, S_proposed, S_star_proposed)
    # third denominator term
    log_likelihood_ratio = log_likelihood_ratio + fB(K, K + 1, S_current, S_star_current)
    log_likelihood_ratio = log_likelihood_ratio + fB(K, K, S_current, S_star_current)
    log_likelihood_ratio = log_likelihood_ratio + fB(K + 1, K + 1, S_current, S_star_current)

    # compute alpha_shuffle
    alpha_shuffle = min(1, exp(log_likelihood_ratio + log_prior_ratio))
    
    #log_print("SHUFFLE proposal", console = FALSE)
    #log_print("rho_current", console = FALSE)
    #log_print(rho_current, console = FALSE)
    #log_print("rho_proposed", console = FALSE)
    #log_print(rho_proposed, console = FALSE)

    if (runif(n = 1) < alpha_shuffle) {
        # accept the shuffle
        #log_print("Shuffle ACCEPTED", console = FALSE)
        return(rho_proposed)
    } else {
        # reject the shuffle
        #log_print("Shuffle REJECTED", console = FALSE)
        return(rho_current)
    }
}






# Build the entire S (sum of edges between clusters) from scratch
# idea: extracting all submatrices needed from G and sum all the elements
# (which are all ones). Do this only for the triangular part, then make it
# symmetric. You need to know just G and the partition rho
get_S_from_G_rho = function(G, rho) {
    
    if (!all(t(G) == G))
        stop("G is not symmetric")
    
    # number of groups
    M = length(rho)
    
    # initialize S matrix
    S = matrix(numeric(M * M), nrow = M, byrow = TRUE)
    
    # indexes of the right bounds of the partition
    bounds = cumsum(rho)
    
    # loop through the groups
    for (l in 1:M) {
        # extract the submatrix and sum all the elements
        for (m in 1:l) {
            start_row = ifelse(m != 1, bounds[m - 1] + 1, 0)
            end_row = bounds[m]
            start_col = ifelse(l != 1, bounds[l - 1] + 1, 0)
            end_col = bounds[l]
            S[l, m] = sum(G[start_row:end_row, start_col:end_col])
            
            if(l == m){
                # the inside connections are now counted twice, correct for it
                S[l, m] = S[l, m] / 2
            } else {
                # otherwise, write to symmetric part of the matrix as well
                S[m, l] = S[l, m]
            }
        }
    }
    return(S)
}

# Build S (sum of edges between clusters) from knowledge of the previous S,
# the G matrix, the previous rho and the new proposed rho.
# Handles all the cases: add, delete, shuffle, i.e. "a new group is added",
# "a group is deleted", "same number of groups, but two groups exchange some
# elements".
get_S_from_G_rho_oldrho_oldS = function(G,rho,oldrho,oldS){
    
    if (!all(t(G) == G))
        stop("G is not symmetric")
    
    # number of groups in new and old rho
    M    = length(   rho)
    oldM = length(oldrho)
    
    # indexes of the right bounds of the partition
    bounds = get_group_indexes(rho)
    
    # groups that needs to be updated with the new rho information
    groups_to_be_refilled = {}
    
    if(M > oldM){ # case Add
        
        # initialize S matrix
        S = matrix(numeric(M * M), nrow = M, byrow = TRUE)
        
        # find the group that has changed
        K = min(which(rho != c(oldrho,NA)))
        
        if(K > 1 && K < oldM){ # in the standard case perform all four
            
            # upper left block
            S[1:(K-1),1:(K-1)] = oldS[1:(K-1),1:(K-1)]
            
            # lower left block
            S[(K+1+1):M,1:(K-1)] = oldS[(K+1):oldM,1:(K-1)]
            
            # upper right block
            S[1:(K-1),(K+1+1):M] = oldS[1:(K-1),(K+1):oldM]
            
            # lower right block
            S[(K+1+1):M,(K+1+1):M] = oldS[(K+1):oldM,(K+1):oldM]
            
        } else if(K == 1){ # bring to the new S only the lower right block
            
            # lower right block
            S[(K+1+1):M,(K+1+1):M] = oldS[(K+1):oldM,(K+1):oldM]
            
        } else if(K == oldM){ # bring to the new S only the upper left block
            
            # upper left block
            S[1:(K-1),1:(K-1)] = oldS[1:(K-1),1:(K-1)]
            
        }
        
        groups_to_be_refilled = c(K,K+1)
        
    } else if(M < oldM) { # case Delete
        
        # initialize S matrix
        S = matrix(numeric(M * M), nrow = M, byrow = TRUE)
        
        # find the group that has changed
        K = min(which(oldrho != c(rho,NA)))
        
        if(K > 1 && K+1 < oldM){ # in the standard case perform all four
            
            # upper left block
            S[1:(K-1),1:(K-1)] = oldS[1:(K-1),1:(K-1)]
            
            # lower left block
            S[(K+1):M,1:(K-1)] = oldS[(K+1+1):oldM,1:(K-1)]
            
            # upper right block
            S[1:(K-1),(K+1):M] = oldS[1:(K-1),(K+1+1):oldM]
            
            # lower right block
            S[(K+1):M,(K+1):M] = oldS[(K+1+1):oldM,(K+1+1):oldM]
            
        } else if(K == 1){ # bring to the new S only the lower right block
            
            # lower right block
            S[(K+1):M,(K+1):M] = oldS[(K+1+1):oldM,(K+1+1):oldM]
            
        } else if(K+1 == oldM){ # bring to the new S only the upper left block
            
            # upper left block
            S[1:(K-1),1:(K-1)] = oldS[1:(K-1),1:(K-1)]
            
        }
        
        groups_to_be_refilled = c(K)
        
    } else { # case Shuffle
        
        S = oldS
        
        # find the group that has changed
        K = min(which(oldrho != rho))
        
        # set to zero the columns and rows of the shuffled groups
        S[K:(K+1),1:M] = 0 # row K to K+1
        S[1:M,K:(K+1)] = 0 # column K to K+1
        
        groups_to_be_refilled = c(K,K+1)
        
    }
    
    # loop through the groups
    for (l in groups_to_be_refilled) {
        # extract the submatrix and sum all the elements
        for (m in 1:M) {
            start_row = ifelse(m != 1, bounds[m - 1] + 1, 0)
            end_row = bounds[m]
            start_col = ifelse(l != 1, bounds[l - 1] + 1, 0)
            end_col = bounds[l]
            
            S[l, m] = sum(G[start_row:end_row, start_col:end_col])
            if(l == m){
                # the inside connections are now counted twice, correct for it
                S[l, m] = S[l, m] / 2
            } else {
                # otherwise, write to symmetric part of the matrix as well
                S[m, l] = S[l, m]
            }
        }
    }
    
    return(S)
}


#' Get Non-Edges from S and partition
#'
#' @param S MxM matrix, where M is number of groups, containing the sum of edges.
#' @inheritParams update_partition
#'
#' @return a positive scalar indicating the number of non-edges between two groups
#' computed as the possible number of edges between two groups
#' (depending on group cardinality) minus the effective number of edges
#' (depending on the edges actually present in the current Graph)
#' @export
#'
#' @examples
get_S_star_from_S_and_rho = function(S, rho){
    # number of groups
    M = length(rho)
    
    # initialize S_star matrix
    S_star = matrix(numeric(M * M), nrow = M, byrow = TRUE)
    
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


# auxiliary function to evaluate the beta function for the likelihood ratio
fB_general = function(group1,
                      group2,
                      S,
                      S_star,
                      alpha,
                      beta,
                      log = TRUE) {
    if (log) {
        return(lbeta(alpha + S[group1, group2], beta + S_star[group1, group2]))
    } else{
        return(beta(alpha + S[group1, group2], beta + S_star[group1, group2]))
    }
}

# auxiliary function to evaluate the beta function for the likelihood ratio
fB_zero = function(alpha,
                   beta,
                   log = TRUE) {
    if (log) {
        return(lbeta(alpha, beta))
    } else{
        return(beta(alpha, beta))
    }
}


#' Log Likelihood Ratio
#'
#' @inheritParams update_partition
#' @inheritParams shuffle_partition
#' @param rho_current Current partition.
#' @param rho_proposed Proposed partition.
#' @param alpha Parameter of the Beta of the Graph.
#' @param beta Parameter of the Beta of the Graph.
#'
#' @return The likelihood ratio in log.
#' @export
#'
#' @examples
log_likelihood_ratio = function(alpha_add,
                                weights_a,
                                weights_d,
                                G,
                                rho_current,
                                rho_proposed,
                                choose_add,
                                alpha,
                                beta,
                                changed_group_index) {
    # differentiate delete/merge case
    if (!choose_add) {
        # swap rhos 'cause we're lazy
        temp = rho_current
        rho_current = rho_proposed
        rho_proposed = temp
    }
    
    # number of groups
    M = length(rho_current)
    
    # wrap the general fB into this version where I don't have to specify
    # alpha and beta (doing this for adaptiveness)
    fB = function(group1, group2, S, S_star) {
        return(fB_general(
            group1,
            group2,
            S,
            S_star,
            alpha = alpha,
            beta = beta,
            log = TRUE
        ))
    }
    
    S_current = get_S_from_G_rho(G, rho_current)
    S_star_current = get_S_star_from_S_and_rho(S_current, rho_current)
    
    S_proposed = get_S_from_G_rho(G, rho_proposed)
    S_star_proposed = get_S_star_from_S_and_rho(S_proposed, rho_proposed)
    
    K = changed_group_index
    #K = get_index_changed_group(rho_current, rho_proposed)
    
    # TODO
    # anziché tenere questo termine, metterlo in ognuno di quelli sotto.
    # Inefficiente, ma serve per consentire l'aggiornamento di alpha e beta
    log_ratio = -(M + 1) * fB_zero(alpha, beta)
    
    for (l in 1:(K - 1)) {
        if (l > (K - 1)) break; # needed because R for-loops suck
        # first numerator term
        log_ratio = log_ratio + fB(l, K, S_proposed, S_star_proposed)
        log_ratio = log_ratio + fB(l, K + 1, S_proposed, S_star_proposed)
        # first denominator term
        log_ratio = log_ratio - fB(l, K, S_current, S_star_current)
    }
    
    for (m in (K + 2):(M + 1)) {
        if (m > (M + 1)) break; # needed because R for-loops suck
        # second numerator term
        log_ratio = log_ratio + fB(K, m, S_proposed, S_star_proposed) +
                                fB(K + 1, m, S_proposed, S_star_proposed)
    }
    
    # third numerator term
    log_ratio = log_ratio + fB(K, K + 1, S_proposed, S_star_proposed) +
                            fB(K, K, S_proposed, S_star_proposed) +
                            fB(K + 1, K + 1, S_proposed, S_star_proposed)
    
    for (m in (K + 1):M) {
        if (m > M) break; # needed because R for-loops suck
        # second denominator term
        log_ratio = log_ratio - fB(K, m, S_current, S_star_current)
    }
    
    # third denominator term
    log_ratio = log_ratio - fB(K, K, S_current, S_star_current)
    
    # in the delete/merge case we have to invert everything
    if (!choose_add) {
        log_ratio = -log_ratio
    }
    
    return(log_ratio)
    
}


#' Get index of the changed group between two partitions
#'
#' @param rho_current Current partition in the form of group cardinalities.
#' @param rho_proposed Proposed partition in the form of group cardinalities.
#'
#' @return Index of the group that has been affected by a change (split/merge/shuffle) from the current patition to the proposed one.
#' @export
#'
#' @examples
get_index_changed_group = function(rho_current, rho_proposed) {
    
    # indexes of the changepoints in the current partition
    cp_idxs_current = get_group_indexes(rho_current)
    
    # move from the rho representation to r representation
    # i.e. from c(2,3) to c(0,1,0,0,1)
    # both for the current and the proposed partition
     
    current_r = rho_to_r(rho_current)
    proposed_r = rho_to_r(rho_proposed)
    
    # now by making the difference I can extract the index
    # of the new changepoint (or the one deleted)
    # there's no need for an absolute value because in the delete move
    # they have been fictitiously swapped to avoid repeating code
    tau = which.max(abs(proposed_r - current_r))
    
    # take all the indexes of the cp smaller than tau
    temp = which(cp_idxs_current < tau)
    if(length(temp) == 0){
        temp = c(0)
    }
    
    # take the last element of this list to have the
    # index of the group affected
    K = temp[length(temp)] + 1
    
    return(K)
}


# TODO Add documentation
log_prior_ratio = function(theta_prior,
                          sigma_prior,
                          rho_current,
                          rho_proposed,
                          choose_add,
                          changed_group_index)
{
    # differentiate delete/merge case
    if (!choose_add) {
        # swap rhos 'cause we're lazy
        temp = rho_current
        rho_current = rho_proposed
        rho_proposed = temp
    }
    
    # number of groups in the current partition
    M = length(rho_current)
    
    K = changed_group_index
    #K = get_index_changed_group(rho_current,rho_proposed)
    
    # compute the prior ratio
    log_ratio = - log(M) + log(theta_prior + M * sigma_prior)
                + lpochhammer(1 - sigma_prior, rho_proposed[K] - 1)
                + lpochhammer(1 - sigma_prior, rho_proposed[K + 1] - 1)
                - lpochhammer(1 - sigma_prior, rho_proposed[K] + rho_proposed[K + 1] - 1)
                + lfactorial(rho_proposed[K] + rho_proposed[K + 1])
                - lfactorial(rho_proposed[K])
                - lfactorial(rho_proposed[K + 1])
    
    
    # in the delete/merge case we have to invert everything
    if (!choose_add) {
        log_ratio = -log_ratio
    }
    
    return(log_ratio)
}




#' Compute weights of the full conditional of theta
#' 
#' For further details see Proposition 1 Martinez and Mena (2014).
#'
#' @param p number of nodes
#' @param sigma_prior other parameter used to compute the prior ratio
#' @param k number of groups
#' @param j index of the iteration for which we are computing the weight
#' @param f value drawn from Exp(theta+1)
#' @param z value drawn from Be(theta+2,n)
#' @param c prior first parameter of the shifted gamma
#' @param d prior second parameter of the shifted gamma
#'
#' @return
#' @export
#'
#' @examples
compute_weights_theta <- function(c, d, p, sigma_prior, k, j, f, z) {
    abs_stir <- abs_stirling_number_1st(k, j)
    num <- (
        (p - sigma_prior) * (p + 1 - sigma_prior) * abs_stirling_number_1st(k - 1, j) +
            (2 * p + 1 - 2 * sigma_prior) * sigma_prior * abs_stirling_number_1st(k - 1, j - 1) +
            (sigma_prior ^ 2) * abs_stirling_number_1st(k - 1, j - 2)
    ) * gamma(c + j)

    denom <- (sigma_prior * (d + f - log(z)))^j
    return(num / denom)
}



#' Full-conditional for theta
#' 
#' For further details see Proposition 1 Martinez and Mena (2014).
#' 
#' TODO COMPLETE DOCUMENTATION
#'
#' @param c First parameter of the shifted gamma prior
#' @param d Second parameter of the shifted gamma prior
#' @param candidate proposed value for theta
#' @param k number of groups
#'
#' @return scalar value for theta at the current iteration
#' @export
#'
#' @examples
full_conditional_theta <- function(prior_c, prior_d, candidate, k, p, sigma_prior){
    weights_gamma <- rep(0,k+2)
    z = rbeta(1,candidate + 2, p)
    f = rexp (1,candidate + 1)
    
    for (j in 0:(k+1)){
        weight_j = compute_weights_theta(prior_c, prior_d, p, sigma_prior, k, j, f, z)
        weights_gamma[j] = weight_j
    }

    # normalizing the weights
    if(sum(weights_gamma) != 0) {
        weights_gamma = weights_gamma / sum(weights_gamma)
    }
    
    component = min(which(cumsum(weights_gamma) > runif(n = 1)))
    
    theta = shifted_gamma(prior_c + (component - 1), prior_d + f - log(z), -sigma_prior)
    return(theta)
}



#' Full-conditional for sigma
#' 
#' For further details see Section 4 Martinez and Mena (2014).
#'
#' @param sigma value of sigma to be updated
#' @param theta value of theta
#' @param rho Current partition
#' @param a First parameter of the distribution
#' @param b Second parameter of the distribution
#' @param c Third parameter of the distribution
#' @param d Fourth parameter of the distribution
#'
#' @return the value for sigma for the current iteration
#' @export
#'
#' @examples
full_conditional_sigma <- function(sigma, theta, rho, a, b, c, d){
    
    # number of groups
    M <- length(rho)

    # first product term
    log_prod_1 <- 0

    for (i in 1:(M - 1)) {
        if (i > (M - 1)) break; # needed because R for-loops suck
        log_prod_1 <- log_prod_1 + log(theta + i * sigma)
    }

    # second product term
    log_prod_2 <- 0
    for (i in 1:M) {
        log_prod_2 <- log_prod_2 + lpochhammer((1 - sigma), (rho[i] - 1))
    }

    # final output
    output <- (a - 1) * log(sigma) + (b - 1) * log(1 - sigma) +
        (c - 1) * log(theta + sigma) + log(exp(-d * sigma)) +
        log_prod_1 + log_prod_2

    return(output)
}

#' Set options for the Gibbs sampler
#'
#' @param sigma_prior_0 initial parameter of the prior (5) from Martinez and Mena (2014).
#' @param sigma_prior_parameters list of 4 parameters for updating sigma_prior, e.g. list("a"=1,"b"=1,"c"=1,"d"=1). For further details see Section 4 in Martinez and Mena (2014).
#' @param theta_prior_0 initial parameter of the prior (5) from Martinez and Mena (2014).
#' @param theta_prior_parameters list of 2 parameters for updating theta_prior, e.g. list("c"=1,"d"=1). For further details see Proposition 1 in Martinez and Mena (2014).
#' @param rho0 initial partition (e.g. c(150,151))
#' @param weights_a0 weights for choosing the candidate node in an add/split move.
#' @param weights_d0 weights for choosing the candidate node in a delete/merge move.
#' @param alpha_target target acceptance rate of the split and merge Metropolis-Hastings.
#' @param beta_mu expected value for the Beta prior of the graph.
#' @param beta_sig2 variance for the Beta prior of the graph. Must be between 0 and 0.25.
#' @param d parameter of the G-Wishart. Default is 3.
#' @param alpha_add probability of choosing an add/split move over a delete/merge move. Default is 0.5.
#' @param adaptation_step adaptation step for tweaking how much the weights are updated each time. Fur further details see Section 4.2 in Benson and Friel (2018).
#' @param update_sigma_prior boolean for choosing whether to update sigma_prior or not. Default is TRUE.
#' @param update_theta_prior boolean for choosing whether to update theta_prior or not. Default is TRUE.
#' @param update_weights boolean for choosing whether to update weights or not. Default is TRUE.
#' @param update_partition boolean for choosing whether to update the partition or not. Default is TRUE.
#' @param update_graph boolean for choosing whether to update the graph or not. Default is TRUE.
#' @param perform_shuffle boolean for choosing whether to perform shuffle or not. Default is TRUE.
#'
#' @return a list with all options correctly set for working with the Gibbs sampler.
#' @export
#'
#' @examples
set_options = function(sigma_prior_0,
                       sigma_prior_parameters,
                       theta_prior_0,
                       theta_prior_parameters,
                       rho0,
                       weights_a0,
                       weights_d0,
                       alpha_target,
                       beta_mu,
                       beta_sig2,
                       d=3,
                       alpha_add=0.5,
                       adaptation_step,
                       update_sigma_prior=TRUE,
                       update_theta_prior=TRUE,
                       update_weights=TRUE,
                       update_partition=TRUE,
                       update_graph=TRUE,
                       perform_shuffle=TRUE
                       ) {
  
    options = list(
        "sigma_prior_0"          = sigma_prior_0,
        "sigma_prior_parameters" = sigma_prior_parameters,
        "theta_prior_0"          = theta_prior_0,
        "theta_prior_parameters" = theta_prior_parameters,
        "rho0"                   = rho0,
        "weights_a0"             = weights_a0,
        "weights_d0"             = weights_d0,
        "alpha_target"           = alpha_target,
        "beta_mu"                = beta_mu,
        "beta_sig2"              = beta_sig2,
        "d"                      = d,
        "alpha_add"              = alpha_add,
        "adaptation_step"        = adaptation_step,
        "update_sigma_prior"     = update_sigma_prior,
        "update_theta_prior"     = update_theta_prior,
        "update_weights"         = update_weights,
        "update_partition"       = update_partition,
        "update_graph"           = update_graph,
        "perform_shuffle"        = perform_shuffle
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
    if(!(var > 0 && var < mu*(1-mu)))
        stop("The variance of the Beta must be between 0 and beta_mu*(1-beta_mu)")

    alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
    beta <- alpha * (1 / mu - 1)
    return(list(alpha = alpha, beta = beta))
}





#' Gibbs sampler
#'
#' @param data an n x p matrix of data.
#' @param niter desired number of effective iterations.
#' @param nburn number of iterations to be burned.
#' @param thin keep all multiples of thin.
#' @param options all the parameters necessary to run the Gibbs sampler.
#' @param seed seed for reproducibility.
#' @param print print the progressbar. Default to TRUE.
#'
#' @return
#' @export
#'
#' @examples
Gibbs_sampler = function(data,
                         niter,
                         nburn,
                         thin,
                         options,
                         seed=1234,
                         print=TRUE)
{
    n = nrow(data) # number of observations
    p = ncol(data) # number of nodes
    n_total_iter = nburn + niter * thin # total iterations to be made
    
    # dynamic parameters
    sigma_prior            = options$sigma_prior_0
    theta_prior            = options$theta_prior_0
    rho                    = options$rho0
    weights_a              = options$weights_a0
    weights_d              = options$weights_d0
    adaptation_step        = options$adaptation_step
    sigma_prior_parameters = options$sigma_prior_parameters
    theta_prior_parameters = options$theta_prior_parameters
    
    # constant parameters
    alpha_add    = options$alpha_add
    alpha_target = options$alpha_target
   
    # parameter for the Wishart
    d = options$d
    
    # parameters for the Beta
    beta_mu   = options$beta_mu
    beta_sig2 = options$beta_sig2
    

    # checks
    if(sum(rho) != p)
        stop("The partition rho must sum to the number of variables p")
    if(d < 3)
        stop("The Wishart's d parameter must be greater or equal than 3")
    if(!(beta_mu > 0 && beta_mu < 1))
        stop("The mean of the Beta must be between 0 and 1")
    if(!(beta_sig2 > 0 && beta_sig2 < beta_mu*(1-beta_mu)))
        stop("The variance of the Beta must be between 0 and beta_mu*(1-beta_mu)")
    if(length(weights_a) != (p-1) || length(weights_d) != (p-1))
        stop("The number of elements in the weights vectors must be equal to p-1")
    if(!(adaptation_step > 0))
        stop("The adapation step h must be positive")
    if(!(alpha_add > 0 && alpha_add < 1))
        stop("The probability of choosing an add move alpha_add must be between 0 and 1")
    if(!(alpha_target > 0 && alpha_target < 1))
        stop("The target acceptance rate of the Metropolis-Hastings alpha_target must be between 0 and 1")

    t_over_p = n_total_iter / p
  
    beta_params = estimate_Beta_params(beta_mu, beta_sig2)
    
    # define structure to save sampled values
    save_res = list(
        G = vector("list", length = niter),
        K = vector("list", length = niter),
        rho = vector("list", length = niter),
        accepted = vector("numeric", length = niter),
        S = vector("list", length = niter),
        sigma = vector("numeric", length = niter),
        theta = vector("numeric", length = niter)
        )
    
    # initialize iteration counter
    it_saved = 0
    
    # initialize progress bar
    if(print){
        pb = txtProgressBar(min=1, max=n_total_iter, initial=1, style=3)
    }
    
    # initialize the sum of the weights of the graphs
    total_weights = 0
    total_K = matrix(0,p,p)
    
    # initialize the sum of all graphs
    total_graphs = matrix(0,p,p)
    g.start = "empty"

    # save start time for measuring execution time
    start_time = Sys.time()
    
    # start the simulation
    for(iter in 1:n_total_iter){
        #log_print("Iter:", console = FALSE)
        #log_print(iter, console = FALSE)

        # update graph
        if (options$update_graph){
            
            # we run a single iteration of BDgraph with iter = 1 and burnin = 0
            output = bdgraph(
                data,
                rho,
                n,
                method = "ggm",
                algorithm = "bdmcmc",
                iter = 1,
                burnin = 0,
                not.cont = NULL,
                g.prior = 0.5,
                df.prior = d,
                CCG_D = NULL,
                g.start = g.start,
                jump = NULL,
                save = TRUE,
                print = 1000,
                cores = NULL,
                threshold = 1e-8
            )


            # extract adjacency matrix G
            last_G = output$last_graph
            # update for the next iteration
            g.start = last_G
            # extract precision matrix K
            last_K = output$last_K

            if(niter > nburn){ # only if niter > nburn right? TODO what about burnin? Siamo 'sgor?
                # update total_weights
                total_weights = total_weights + output$all_weights
                # update total_graphs taking into consideration the weights
                total_graphs = total_graphs + output$last_graph * output$all_weights
                total_K = total_K + output$last_K * output$all_weights
            }
        }
        
        if(options$update_partition){
            list_output_update_partition = update_partition(rho,
                                                            alpha_add,
                                                            weights_a,
                                                            weights_d,
                                                            theta_prior,
                                                            sigma_prior,
                                                            last_G,
                                                            beta_params)
            
            rho = list_output_update_partition$rho_updated

            # it makes sense to perform the adaptive step only if we're updating the partition
            # update the single weight at the point only if the move has been accepted
            if(options$update_weights && list_output_update_partition$accepted){
                if(list_output_update_partition$choose_add){
                    weights_a = update_weight(
                        weights_a,
                        list_output_update_partition$candidate,
                        adaptation_step,
                        t_over_p,
                        alpha_add,
                        alpha_target
                    )
                } else {
                    weights_d = update_weight(
                        weights_d,
                        list_output_update_partition$candidate,
                        adaptation_step,
                        t_over_p,
                        1 - alpha_add,
                        alpha_target
                    )
                }
            }
        }
        
        if(options$perform_shuffle){
            rho = shuffle_partition(rho, last_G, sigma_prior, beta_params$alpha, beta_params$beta)
        }
        
        if(options$update_sigma_prior){
            candidate <- runif(1,max(0,-theta_prior),1)
            alpha_MH <- full_conditional_sigma(candidate,
                                               theta_prior,
                                               rho,
                                               sigma_prior_parameters$a,
                                               sigma_prior_parameters$b,
                                               sigma_prior_parameters$c,
                                               sigma_prior_parameters$d) -
                        full_conditional_sigma(sigma_prior,
                                               theta_prior,
                                               rho,
                                               sigma_prior_parameters$a,
                                               sigma_prior_parameters$b,
                                               sigma_prior_parameters$c,
                                               sigma_prior_parameters$d)
            
            if(log(runif(1)) <= min(alpha_MH,log(1))){
                sigma_prior = candidate
            }else{
                sigma_prior = sigma_prior
            }
        }

        if(options$update_theta_prior) {
            theta_prior = full_conditional_theta(
                theta_prior_parameters$c,
                theta_prior_parameters$d,
                theta_prior,
                length(rho),
                p,
                sigma_prior
            )
        }
        
        if (options$update_graph){
            last_S = get_S_from_G_rho(last_G,rho)
        }
        
        # save results only on thin iterations
        # (i.e. only save multiples of thin)
        if(iter > nburn && (iter - nburn) %% thin == 0) {
            it_saved = it_saved + 1
            if(options$update_graph){
                # cumulative precision matrix K and probability of inclusion links
                save_res$K[[it_saved]] = total_K / total_weights
                save_res$G[[it_saved]] = total_graphs / total_weights
            }
            if(options$update_partition){
                save_res$rho[[it_saved]] = rho
                save_res$accepted[[it_saved]] = list_output_update_partition$accepted
            }
            if(options$update_sigma_prior){
                save_res$sigma[[it_saved]] = sigma_prior
            }
            if(options$update_theta_prior){
                save_res$theta[[it_saved]] = theta_prior
            }
            save_res$S[[it_saved]] = last_S
        }
        
        if(print){
            setTxtProgressBar(pb, iter)
        }
        #log_print("last_G:", console = FALSE)
        #log_print(last_G, console = FALSE)
        #log_print("last_S:", console = FALSE)
        #log_print(last_S, console = FALSE)
        #log_print("---------------------------------------------------------------------", console = FALSE)
    }

    save_res$execution_time = Sys.time() - start_time
    
    #output = list( sample_graphs = sample_graphs, graph_weights = graph_weights, K_hat = "empty",
    #                   all_graphs = niter-nburn, all_weights = all_weights, last_graph = last_G,
    #                   last_K = last_K, last_Theta = last_Theta )
    #
    
    if(print){close(pb)}
    return(save_res)
}

