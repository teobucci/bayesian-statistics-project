#' Generate_BlockDiagonal
#' This function generate data by constructing a block diagonal matrix. 
#' The idea is taken from Sun et at. (2014) -  Adaptive Variable Clustering in Gaussian Graphical Models
#' We add an edge between two variables with probability p_in if they are in the same cluster, otherwise with
#' probability p_out. The corresponding value of the precision matrix is set equal to elem_out. To ensure positive
#' definiteness, we increase its eigenvalues summing the absolute value of the minimum eigenvalue plus min_eigenval_correction
#' Given the precision matrix, the covariance matrix, the data and the graph are generate accordingly. 
#' Note, the precision matrix is not drawn from a G-Wishart distribution
#' @param n [scalar] sample size
#' @param z_true [vector] cluster allocation vector
#' @param p_in [scalar] probability of having an edge between variables in the same cluster
#' @param p_out [scalar] probability of having an edge between variables in the different cluster
#' @param elem_out [scalar] value of elements of the precision matrix outside the diagonal
#' @param min_eigenval_correction [scalar] value to be added to the minimum eigenvalue
#' @param seed [integer] random seed
#'
#' @export
Generate_BlockDiagonal = function(n,
                                  z_true,
                                  p_in = 1,
                                  p_out = 0,
                                  elem_out = 5,
                                  min_eigenval_correction = 3,
                                  seed = 25254) {
    
  set.seed(seed)
  p = length(z_true) #set p
  # Generate precision and covariance matrices
  Omega = matrix(0,nrow = p,ncol = p)
  for(i in 1:(p-1)){
    for(j in (i+1):p){
      
      if(z_true[i]==z_true[j]){
        # draw with prob p_in
        if(runif(n=1)<p_in)
          Omega[i,j] = elem_out
      }else{
        # draw with prob p_out
        if(runif(n=1)<p_out)
          Omega[i,j] = elem_out
      }
      
    }
  }
  Omega = Omega + t(Omega) # define lower diagonal part
  min_eig = eigen(Omega)$values[p] # compute eigenvalues, may be negative or 0
  Omega_true = Omega + diag(x=abs(min_eig) + min_eigenval_correction,nrow = p)  # ensure positive definiteness
  Sigma_true = solve(Omega_true) # compute covariance matrix
  
  # Generate data
  data = matrix(NA,nrow = n, ncol = p)
  data = t(apply(data, 1, FUN = function(x){ mvtnorm::rmvnorm(n = 1, sigma = Sigma_true)  }))
  U = t(data)%*%data
  
  # Compute underlying graph
  G = Omega_true
  G[abs(G)>0] = 1
  
  return(list("data"=data,"U"=U,"Prec"=Omega_true,"Cov"=Sigma_true,"Graph"=G))
}


#' Generate_Block
#' This function generalizes Generate_BlockDiagonal to have more flexibility.
#' First, a matrix Theta of size N_clust x N_clust is generated. Elements of Theta define if a block between two clusters must be
#' added or not. Given Theta, elements in each block are drawn with probability p_inside_block while elements outside are drawn with
#' probability p_outside_block.
#' @inheritParams Generate_BlockDiagonal
#' @param p_block_diag [scalar] probability of having a block along the diagonal
#' @param p_block_extra_diag [scalar] probability of having a block outside the diagonal
#' @param p_inside_block [scalar] probability of having an edge between variables in the same cluster
#' @param p_outside_block [scalar] probability of having an edge between variables in the different cluster
#' @export
Generate_Block = function(n,
                          z_true,
                          p_block_diag = 1,
                          p_block_extra_diag = 0.5,
                          p_inside_block = 1,
                          p_outside_block = 0,
                          elem_out = 5,
                          min_eigenval_correction = 3,
                          seed = 25254) {
    
  set.seed(seed)
  p = length(z_true) #set p
  counts_true = as.vector(table(z_true))
  Nclust_true  = length(counts_true)
  
  # Generate which block should be connected
  Theta = matrix(0,Nclust_true,Nclust_true)
  for(u in 1:Nclust_true){
    for(v in u:Nclust_true){
      
      if(u==v){
        # draw with prob p_block_diag
        if(runif(n=1)<p_block_diag)
          Theta[u,v] = 1
      }else{
        # draw with prob p_block_extra_diag
        if(runif(n=1)<p_block_extra_diag)
          Theta[u,v] = 1
      }
      
    }
  }
  # Generate precision and covariance matrices
  Omega = matrix(0,nrow = p,ncol = p)
  for(i in 1:(p-1)){
    for(j in (i+1):p){
      
      if( Theta[ z_true[i],z_true[j] ]==1 ){
        # draw with prob p_inside_block
        if(runif(n=1)<p_inside_block)
          Omega[i,j] = elem_out
      }else{
        # draw with prob p_outside_block
        if(runif(n=1)<p_outside_block)
          Omega[i,j] = elem_out
      }
      
    }
  }
  Omega = Omega + t(Omega) # define lower diagonal part
  min_eig = eigen(Omega)$values[p] # compute eigenvalues, may be negative or 0
  Omega_true = Omega + diag(x=abs(min_eig) + min_eigenval_correction,nrow = p)  # ensure positive definiteness
  Sigma_true = solve(Omega_true) # compute covariance matrix
  
  # Generate data
  data = matrix(NA,nrow = n, ncol = p)
  data = t(apply(data, 1, FUN = function(x){ mvtnorm::rmvnorm(n = 1, sigma = Sigma_true)  }))
  U = t(data)%*%data
  
  # Compute underlying graph
  G = Omega_true
  G[abs(G)>0] = 1
  
  return(list("data"=data,"U"=U,"Prec"=Omega_true,"Cov"=Sigma_true,"Graph"=G,"Theta"=Theta))
}

# ===========================================
# =           EXAMPLE OF PLOTTING           =
# ===========================================

# # Plot Prec
# ACutils::ACheatmap(
#     sim$Prec,
#     use_x11_device = F,
#     horizontal = F,
#     main = "Precision",
#     center_value = NULL,
#     col.upper = "black",
#     col.center = "grey50",
#     col.lower = "white"
# )

# # Plot Cov
# ACutils::ACheatmap(
#     sim$Cov,
#     use_x11_device = F,
#     horizontal = F,
#     main = "Covariance",
#     center_value = NULL,
#     col.upper = "black",
#     col.center = "grey50",
#     col.lower = "white"
# )

# # Plot empirical estimate
# ACutils::ACheatmap(
#     solve(cov(sim$data)),
#     use_x11_device = F,
#     horizontal = F,
#     main = "Empirical estimate",
#     center_value = NULL,
#     col.upper = "black",
#     col.center = "grey50",
#     col.lower = "white"
# )

# # Plot graph
# ACutils::ACheatmap(
#     sim$G,
#     use_x11_device = F,
#     horizontal = F,
#     main = "Graph",
#     center_value = NULL,
#     col.upper = "black",
#     col.center = "grey50",
#     col.lower = "white"
# )

# # Plot Theta
# ACutils::ACheatmap(
#     sim$Theta,
#     use_x11_device = F,
#     horizontal = F,
#     main = "Theta",
#     center_value = NULL,
#     col.upper = "black",
#     col.center = "grey50",
#     col.lower = "white"
# )
