##' get group indexes
#'
#' @param rho partition in che compact version ef c(1,3,3)
#'
#' @return a vector whose components are
#' the increasing indexes of the elements where each group starts
#' e.g. if rho=c(1,3,3), the output will be c(1,4,6),
#' meaning that 1 is the index of the first element of the first group,
#' 4 is the index of the first element of the second group, 
#' 6 is the index of the first element of the third group 
#' 
#' 
#' @export
#'
#' @examples
get_group_indexes = function(rho){
    return(cumsum(rho))
}




#' Rho to changepoint function
#' #TODO Add checks
#'
#' @param rho 
#'
#' @return
#' @export
#'
#' @examples
rho_to_r = function(rho){
    group_indexes = get_group_indexes(rho)
    total_n = sum(rho)
    r <- numeric(total_n)
    for(i in group_indexes){
        r[i]=1
    }
    return(r)
}



#' Rho to z function
#' #TODO Add checks
#'
#' @param rho 
#'
#' @return
#' @export
#'
#' @examples
rho_to_z = function(rho){
    total_n = sum(rho)
    n_groups <- length(rho)
    z <- {}
    for(i in 1:n_groups){
        z <- c(z, rep(i,rho[i]))
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
    #TODO finishing touch: scriverla in C
    if (n < 0)
        stop("The pochhammer operator doesn't allow n < 0")
    if (x < 0)
        stop("The pochhammer operator doesn't allow x < 0")
    if (x == 0){
        if(log)
            return(NA)
        else
            return(0)
    }

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



# Utils for the theta parameter---------

#' Absolute value of Stirling number of the first kind (adapted)
#' Computes an adapted version of the Stirling number of the first kind
#'
#' The Stirling number represents the number of ways that we can arrange 
#' k objects around indistinguishable circles of length j
#' 
#' NOTE: Requires gmp package!
#' 
#' @param k First parameter - indicates the overall number of objects
#' @param j Second parameter - indicates the length of the circles (see above)
#'
#' @return a positive scalar indicating the adapted version of the Stirling number of the first kind
#' (i.e. the "unacceptable" values are turned to zeroes)
#' @export
#'
#' @examples
abs_stirling_number_1st <- function(k,j){
    if(k<0){
        stop("Only positive integers allowed for k!")
    }
    if(j <= 0 || j > k){
        abs_stir_num_first_kind=0
    }
    else {
        abs_stir_num_first_kind=(as.numeric(abs(Stirling1(k,j))))
    }
    return(abs_stir_num_first_kind)
}




#' Compute weights of the full conditional of theta 
#'
#' as described on Martinez and Mena (page 14)
#'
#' @param n 
#' @param sigma other parameter used to compute the prior ration
#' @param k 
#' @param j     index of the iteration for which we are computing the weight
#' @param f     value drawn from Exp(θ+1)
#' @param z     value drawn from Be(θ+2,n)
#' @param c     prior first parameter of the shifted gamma
#' @param d     prior second parameter of the shifted gamma
#'
#' @return
#' @export
#'
#' @examples
compute_weights_theta=function(c, d, n, sigma, k, j, f, z){
    abs_stir= abs_stir_num_first_kind(k,j)
    num = (
        (n-sigma)*(n+1-sigma)*abs_stirling_number_1st(k-1,j) + 
            (2*n + 1 - 2*sigma)*sigma*abs_stirling_number_1st(k-1,j-1) +
            (sigma^2)*abs_stirling_number_1st(k-1,j-2)
    ) * gamma(c+j)
    
    denom = (sigma * (d + f - log(z) ) )^j
    return (num/denom)
}


#' Shifted gamma function
#'
#' Computes a shifted gamma, given the parameters and the shift 
#'
#' @param u 
#' @param o 
#' @param sigma_shift 
#'
#' @return
#' @export
#'
#' @examples
shifted_gamma <- function(u,o,sigma_shift){
    rgamma(1,u,o) + sigma_shift
}


