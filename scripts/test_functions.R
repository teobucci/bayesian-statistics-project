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


G = matrix(c(
    0,0,1,0,1,
    0,0,0,1,0,
    1,0,0,1,0,
    0,1,1,0,1,
    1,0,0,1,0
),nrow=5,byrow = T)

rho = c(1, 3, 1)

get_S_from_G_rho(G=G,rho=rho)
