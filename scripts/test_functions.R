#
#Some tests, useful for debugging
#

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










# make a matrix symmetric
makeSymm <- function(m, from_lower = T) {
    if (from_lower)
        m[upper.tri(m)] = t(m)[upper.tri(m)] # from lower
    else
        m[lower.tri(m)] = t(m)[lower.tri(m)] # from upper
    return(m)
}


# Build the entire S (sum of edges between clusters) from scratch
# idea: extracting all submatrices needed from G and sum all the elements
# (which are all ones). Do this only for the triangular part, then make it
# symmetric. You need to know just G and the partition rho
get_S_from_G_rho = function(G, rho) {
    
    # check on G
    if (!all(t(G) == G))
        stop("G is not symmetric")
    
    # number of groups
    M = length(rho)
    
    # initialize S matrix
    S = matrix(numeric(M * M), nrow = M, byrow = T)
    
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
    
    # alternative correction for the diagonal
    # diagonal = col(S) == row(S)
    # S[diagonal] = S[diagonal] / 2
    
    # alternative for making the matrix symmetric
    # S = makeSymm(S)
    
    return(S)
}

# Build S (sum of edges between clusters) from knowledge of the previous S,
# the G matrix, the previous rho and the new proposed rho.
# Handles all the cases: add, delete, shuffle, i.e. "a new group is added",
# "a group is deleted", "same number of groups, but two groups exchange some
# elements".
get_S_from_G_rho_oldrho_oldS = function(G,rho,oldrho,oldS,debug=F){
    
    # check on G
    if (!all(t(G) == G))
        stop("G is not symmetric")
    
    # number of groups in new and old rho
    M    = length(   rho)
    oldM = length(oldrho)
    
    # indexes of the right bounds of the partition
    bounds = cumsum(rho)
    
    # groups that needs to be updated with the new rho information
    groups_to_be_refilled = {}
    
    if(M > oldM){ # case Add
        
        # initialize S matrix
        S = matrix(numeric(M * M), nrow = M, byrow = T)
        
        # find the group that has changed
        K = min(which(rho != c(oldrho,NA)))
        
        if(K > 1 & K < oldM){ # in the standard case perform all four
            
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
        S = matrix(numeric(M * M), nrow = M, byrow = T)
        
        # find the group that has changed
        K = min(which(oldrho != c(rho,NA)))
        
        if(K > 1 & K+1 < oldM){ # in the standard case perform all four
            
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

# test the upper functions
G = matrix(c(
    0,0,1,0,1,
    0,0,0,1,0,
    1,0,0,1,0,
    0,1,1,0,1,
    1,0,0,1,0
),nrow=5,byrow = T)

oldrho = c(1,3,1)
rho = c(1,2,1,1)
oldS = get_S_from_G_rho(G,oldrho)
S_from_scratch = get_S_from_G_rho(G,rho)
S_from_previous_S = get_S_from_G_rho_oldrho_oldS(G,rho,oldrho,oldS,debug=T)

if(all(S_from_scratch == S_from_previous_S)){
    cat("Ce l'hai fatta Jonny")
}




