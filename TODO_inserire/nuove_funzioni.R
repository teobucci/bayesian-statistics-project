# KLdistance for comparing estimated matrix -------------------------------
library(ACutils)

# Kullback-Leibler for zero-mean Gaussian data
?ACutils::KL_dist

# BFDR for selecting estimated graph --------------------------------------

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
BFDR_selection = function (plinks, tol = seq(0.1, 1, by = 0.025), min_rate = 0.05, diag = F)
{
  if(dim(plinks)[1] != dim(plinks)[2])
    stop("plinksinks matrix has to be squared")
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
  result = list()
  result[[1]] = best_soglia_fdr
  result[[2]] = best_graph_fdr
  names(result) = c("best_treshold", "best_truncated_graph")
  return(result)
}


# Clustering --------------------------------------------------------------


# Get labels for each iterations for each data point
#GDFMM$Partition is a (n_iter x n_data) matrix
part_matrix <- GDFMM$Partition
## ---> this is the object where you saved all the partitions


# Compute similarity matrix
sim_matrix <- psm(part_matrix)

# Heatmap, to understand the uncertainty in the similarity matrix
heatmap(sim_matrix)

# Compute final partition by minimazing the Binder loss or the VI loss functions
binder = minbinder(sim_matrix)
VI     = minVI(sim_matrix)

# Compare partitions just in terms of number of clusters and cluster cardinalities
table(binder$cl)
table(VI$cl)
table(real_partition)

# Compute Rand Index
arandi(binder$cl,real_partition, adjust = F)
arandi(VI$cl,real_partition, adjust = T)


