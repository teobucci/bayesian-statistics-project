#
#Functions that are useful to call from other files.
#

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




