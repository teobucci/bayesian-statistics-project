#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]

// HOW TO USE IN R:
// X: matrix, each row is an iteration, each column an obs
//
// M <- psm(X)
// dists <- VI_LB(clean_partition(X), psm_mat = M)
// X[which.min(dists),]

//[[Rcpp::export]]
arma::mat psm(arma::mat M){
  // initialize results
  arma::mat result(M.n_cols, M.n_cols, arma::fill::zeros);
  
  for(arma::uword i = 0; i < M.n_cols; i++){
    for(arma::uword j = 0; j <= i; j++){
      result(i,j) = arma::accu(M.col(i) == M.col(j));
      result(j,i) = result(i,j);
    }
    Rcpp::checkUserInterrupt();
  }
  return(result / M.n_rows);
}

//[[Rcpp::export]]
arma::vec VI_LB(arma::mat C_mat, arma::mat psm_mat){
  
  arma::vec result(C_mat.n_rows);
  double f = 0.0;
  int n = psm_mat.n_cols;
  arma::vec tvec(n);
  
  for(arma::uword j = 0; j < C_mat.n_rows; j++){
    f = 0.0;
    for(arma::uword i = 0; i < n; i++){
      tvec = psm_mat.col(i);
      f += (log2(arma::accu(C_mat.row(j) == C_mat(j,i))) +
        log2(arma::accu(tvec)) -
        2 * log2(arma::accu(tvec.elem(arma::find(C_mat.row(j).t() == C_mat(j,i))))))/n;
    }
    result(j) = f;
    Rcpp::checkUserInterrupt();
  }
  return(result);
}