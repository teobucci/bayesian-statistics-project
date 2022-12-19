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
