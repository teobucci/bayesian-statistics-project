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









# mi serve una funzione che prenda la G del grafo attuale
# e una partizione, e mi restituisca S (ovvero Theta, che con rho attuale già ho)
# gli devo dare la nuova partizione

# per esempio se G è fatta cosi' (5 nodi)
G = matrix(c(
    0,0,1,0,1,
    0,0,0,1,0,
    1,0,0,1,0,
    0,1,1,0,1,
    1,0,0,1,0
),nrow=5,byrow = T)
# metti che la vecchia rho è fatta cosi'
rho = c(
    2,3
)
# devo trovare S
# devo partire in realtà dalla conoscenza di Theta (cioè S)
# perché altrimenti è un mezzo casino ricostruirla tutta
# devo prendere Theta e cambiare 2 colonne e 2 righe (1 concettualmente
# indipendente, l'altra segue) relative ai due nuovi grupi S e S+1

# ora penso al caso ADD, se sono nel caso DELETE prima piango e poi ci penso

# in questo caso Theta potrebbe essere fatta cosi'
c(
    0,3,
    3,2
)

# cavoli però nel caso delete

# [piange]

# per chiarezza di idee prima ne faccio una generica che prende un generico
# G e rho, e mi da Theta/S
# scomodo ricostruire tutto ma intanto entro nell'idea



# make a matrix symmetric
makeSymm <- function(m, from_lower = T) {
    if (from_lower)
        m[upper.tri(m)] = t(m)[upper.tri(m)] # from lower
    else
        m[lower.tri(m)] = t(m)[lower.tri(m)] # from upper
    return(m)
}

get_S_from_G_rho = function(G, rho) {
    
    # check on G
    if (!all(t(G) == G))
        stop("G is not symmetric")
    
    M = length(rho) # number of groups
    S = matrix(numeric(M * M), nrow = M, byrow = T) # initialize S matrix
    bounds = cumsum(rho) # index of the right bounds of the partition
    
    # loop through the groups
    for (l in 1:M) {
        # cat("analysing group", i, "\n")
        # extract the submatrix and sum all the elements
        # the inside connections are now counted twice
        for (m in 1:l) {
            start_row = ifelse(m != 1, bounds[m - 1] + 1, 0)
            end_row = bounds[m]
            start_col = ifelse(l != 1, bounds[l - 1] + 1, 0)
            end_col = bounds[l]
            S[l, m] = sum(G[start_row:end_row, start_col:end_col])
            # cat("connection between", l, "and", m, "=", S[l, m], "\n")
        }
    }
    
    # correct the diagonal
    diagonal = col(S) == row(S)
    S[diagonal] = S[diagonal] / 2
    
    # make symmetric
    S = makeSymm(S)
    
    return(S)
    
    # idealmente una volta che conosco gli indici che sono cambiati
    # devo solo cambiare l'outer loop con gli indici che mi interessano
    # cioè il cavolo perché la matrice S ha un numero diverso di elementi
    # quindi va presa quella vecchia e messa al suo posto lasciando spazio
    # vuoto per la roba nuova
    # che schifo
    # TODO fare questo per evitare di bruciare risorse, già che non ce ne sono troppe
}

# faccio un tentativo di scrivere sta roba senza ricalcolare tutto
# ma partendo da G, vecchia rho e sua S, nuova rho

get_S_from_G_rho_oldrho_oldS = function(G,rho,oldrho,oldS,debug=F){
    # check on G
    if (!all(t(G) == G))
        stop("G is not symmetric")
    
    M = length(rho) # number of groups in new rho
    oldM = length(oldrho) # number of groups in old rho
    S = matrix(numeric(M * M), nrow = M, byrow = T) # initialize S matrix
    bounds = cumsum(rho) # index of the right bounds of the partition
    
    groups_to_be_refilled = {}
    
    if(M > oldM){
        if(debug)
            cat("Entering case ADD\n")
        # caso di add
        
        # find the group that has changed
        # notazione infelice perche' lo chiamiamo anche lui S, va uniformato
        K = min(which(rho != c(oldrho,NA)))
        
        if(debug)
            cat("Group affected K =",K,"\n")
        
        if(K > 1 && K < oldM){ # in the standard case perform all four
            S[1:(K-1),1:(K-1)] = oldS[1:(K-1),1:(K-1)] # upper left block
            S[(K+1+1):M,1:(K-1)] = oldS[(K+1):oldM,1:(K-1)] # lower left block
            S[1:(K-1),(K+1+1):M] = oldS[1:(K-1),(K+1):oldM] # upper right block
            S[(K+1+1):M,(K+1+1):M] = oldS[(K+1):oldM,(K+1):oldM] # lower right block
        } else if(K == 1){ # bring to the new S only the lower right block
            S[(K+1+1):M,(K+1+1):M] = oldS[(K+1):oldM,(K+1):oldM] # lower right block
        } else if(K == oldM){ # bring to the new S only the upper left block
            S[1:(K-1),1:(K-1)] = oldS[1:(K-1),1:(K-1)] # upper left block
        }
        
        
        groups_to_be_refilled = c(K,K+1)
        
    } else if(M < oldM) {
        if(debug)
            cat("Entering case DELETE\n")
        # caso di delete
        
        # find the group that has changed
        # notazione infelice perche' lo chiamiamo anche lui S, va uniformato
        K = min(which(oldrho != c(rho,NA)))
        
        if(debug)
            cat("Group affected K =",K,"\n")
        
        if(K > 1 && K+1 < oldM){ # in the standard case perform all four
            S[1:(K-1),1:(K-1)] = oldS[1:(K-1),1:(K-1)] # upper left block
            S[(K+1):M,1:(K-1)] = oldS[(K+1+1):oldM,1:(K-1)] # lower left block
            S[1:(K-1),(K+1):M] = oldS[1:(K-1),(K+1+1):oldM] # upper right block
            S[(K+1):M,(K+1):M] = oldS[(K+1+1):oldM,(K+1+1):oldM] # lower right block
        } else if(K == 1){ # bring to the new S only the lower right block
            S[(K+1):M,(K+1):M] = oldS[(K+1+1):oldM,(K+1+1):oldM] # lower right block
        } else if(K+1 == oldM){ # bring to the new S only the upper left block
            S[1:(K-1),1:(K-1)] = oldS[1:(K-1),1:(K-1)] # upper left block
        }
        
        groups_to_be_refilled = c(K)
        
    } else {
        if(debug)
            cat("Entering case SHUFFLE\n")
        # caso di shuffle quando sono uguali il numero di partizioni
        
        S = oldS
        
        # find the group that has changed
        # notazione infelice perche' lo chiamiamo anche lui S, va uniformato
        K = min(which(oldrho != rho))
        
        if(debug)
            cat("Group affected K =",K,"\n")
        
        S[K:(K+1),1:M] = 0 # row K to K+1
        S[1:M,K:(K+1)] = 0 # column K to K+1
        
        groups_to_be_refilled = c(K,K+1)
        
    }

    # sistemare i loop qui sotto solo dove serve
    cat("Matrix S before filling\n")
    print(S)
    
    # loop through the groups
    for (l in groups_to_be_refilled) {
        # cat("analysing group", i, "\n")
        # extract the submatrix and sum all the elements
        # the inside connections are now counted twice
        for (m in 1:M) {
            start_row = ifelse(m != 1, bounds[m - 1] + 1, 0)
            end_row = bounds[m]
            start_col = ifelse(l != 1, bounds[l - 1] + 1, 0)
            end_col = bounds[l]
            S[l, m] = sum(G[start_row:end_row, start_col:end_col])
            S[m, l] = S[l, m]
            if(l == m)
                S[l, m] = S[l, m] / 2 # correct the diagonal
            if(debug)
                cat("connection between", l, "and", m, "=", S[l, m], "\n")
        }
    }
    
    
    
    # make symmetric
    # S = makeSymm(S,from_lower=F)
    # ho problemi a capire se devo fare dalla lower o dalla upper, risolvo inserendo sopra e amen
    
    return(S)
    
}


G = matrix(c(
    0,0,1,0,1,
    0,0,0,1,0,
    1,0,0,1,0,
    0,1,1,0,1,
    1,0,0,1,0
),nrow=5,byrow = T)

oldrho = c(1,3,1)
rho = c(1,2,1,1)
oldS = get_S_from_G_rho(G=G,rho=oldrho)
S_method_correct = get_S_from_G_rho(G=G,rho=rho)
S_new_method = get_S_from_G_rho_oldrho_oldS(G,rho,oldrho,oldS,debug=T)

if(all(S_new_method == S_method_correct)){
    cat("Ce l'hai fatta Jonny")
}




