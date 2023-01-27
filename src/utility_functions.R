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
#' @export
#'
#' @examples
get_group_indexes = function(rho){
    return(cumsum(rho))
}

# =============================================================================
# =           FUNCIONS FOR MOVING BETWEEN PARTITION REPRESENTATIONS           =
# =============================================================================

rho_to_r = function(rho){
    group_indexes = get_group_indexes(rho)
    group_indexes = group_indexes[1:(length(group_indexes)-1)]
    r <- numeric(sum(rho)-1)
    r[group_indexes] = 1
    return(r)
}

z_to_rho = function(z){
    return(as.vector(table(z)))
}

z_to_r = function(z){
    return(diff(z))
}

rho_to_z = function(rho){
    z = numeric()
    
    for(i in seq_along(rho)){
        z = c(z, rep(i, rho[i]))
    }
    return(z)
}


# TODO
r_to_rho = function(r){
    return()
}

r_to_z = function(r){
    return()
}



#' Compute the rising factorial (also called Pochhammer symbol)
#'
#' Peturns the value of Pochhammer's symbol calculated as
#' \deqn{(x)_n = x (x+1) \cdots (x+n-1)}
#'
#'
#' @param x numeric value for the argument of the symbol
#' @param n integer value for the number of terms in the symbol
#' @param log boolean value, if TRUE the rising factorial is returned in log. Default is TRUE.
#'
#' @return the rising factorial of x with n terms.
#' @export
#'
#' @examples
lpochhammer <- function(x, n, log = TRUE) {
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




#' Absolute value of Stirling number of the first kind (adapted)
#' Computes an adapted version of the Stirling number of the first kind
#'
#' The Stirling number represents the number of ways that we can arrange
#' k objects around indistinguishable circles of length j
#'
#' NOTE: Requires gmp package.
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
    if (j == 0 && k == 0) {
        return(1)
    }
    
    if (k < 0) {
        stop("In computing the Stirling number of the first kind, k must be greater or equal than 0.")
    }
    if (j <= 0 || j > k) {
        abs_stir_num_first_kind = 0
    }
    else{
        abs_stir_num_first_kind = (as.numeric(abs(gmp::Stirling1(k, j))))
    }
    return(abs_stir_num_first_kind)
}


#' Shifted gamma function
#'
#' Computes a shifted gamma, given the parameters and the shift.
#' Z ~ shiftedGamma(alpha,beta,mu) is and only if Z-mu ~ Gamma(alpha,beta)
#'
#' @param alpha first parameter of the Gamma
#' @param beta second parameter of the Gamma
#' @param mu shift parameter
#'
#' @return
#' @export
#'
#' @examples
shifted_gamma <- function(alpha, beta, mu, n = 1) {
    rgamma(n, alpha, beta) + mu
}



#' Bayesian FDR Analysis
#'
#' \loadmathjax Given the plinks matrix, this utility computes the False Discovery Rate Index, forcing the false discovery rates to be less than \code{min_rate.}
#' @param plinks matrix containing the posterior inclusion probability for each link. It has to be upper triangular. Its dimension depends on the type of graph it represents.
#' It is indeed possible to pass a \mjseqn{p \times p} matrix, or a \mjseqn{n\_groups \times n\_groups}.
#' @param tol sequence of tolerances to be tested trying to select a graph truncating \code{plinks} at that value.
#' @param min_rate fix false discoveries to remain under this selected threshold.
#' @param diag boolean, if the diagonal of \code{plinks} has to be included in the computations. Set \code{FALSE} if the graph is in complete form, set \code{TRUE} for block graphs.
#'
#' @return a list of two elements: best_threshold, the best value of tol according to this analysis.
#' best_truncated_graph, the proposed posterior graph according to the analysis.
#' @export
BFDR_selection = function (plinks, tol = seq(0.1, 1, by = 0.025), min_rate = 0.05, diag = FALSE)
{
    if(dim(plinks)[1] != dim(plinks)[2])
        stop("plinks matrix must to be squared")
    p = dim(plinks)[1]
    plinks_vet = plinks[upper.tri(plinks, diag = diag)]
    if (any(tol > max(plinks_vet)))
        tol <- tol[-which(tol > max(plinks_vet))]
    if (is.null(tol))
        stop("No feasible tolerances")
    FDR = rep(0, length(tol))
    for (i in 1:length(tol)) {
        tolerance <- tol[i]
        above_tr = plinks_vet[plinks_vet >= tolerance]
        FDR[i] = sum(1 - above_tr)/length(above_tr)
    }
    if(FDR[1] < min_rate) {
        best_soglia_fdr = tol[1]
    }
    else for (i in 2:length(FDR)) {
        if (FDR[i] < min_rate)
            (break)()
    }
    best_soglia_fdr = tol[i]
    best_graph_fdr = matrix(0, p, p)
    best_graph_fdr[plinks >= best_soglia_fdr] = 1
    result = list(
        "best_treshold"=best_soglia_fdr,
        "best_truncated_graph"=best_graph_fdr
    )
    return(result)
}
