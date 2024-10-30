/* ---------------------------------------------------
 Laplace approximation in C++ for lpsmc
 Copyright Oswaldo Gressani. All rights reserved.
 ------------------------------------------------------*/

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
// [[Rcpp::export]]

List Rcpp_Laplace2(NumericVector lat0, double v, int K,
                  Function Dloglik, Function D2loglik, Function Qv){

  int iter = 0;
  double tol = 1e-3;
  double dist = 3;
  int dimlat = lat0.length();
  double thetaconstr = 1;
  NumericMatrix Q = Qv(v);
  NumericVector gloglik;
  NumericMatrix g2loglik;
  List checkPD;
  LogicalVector condcheck;
  arma::mat Precmat(dimlat,dimlat);
  arma::mat Covmat(dimlat,dimlat);
  arma::vec latnew(dimlat);
  CharacterVector PDnote;

  while(dist > tol){
    gloglik = Dloglik(lat0);
    g2loglik = D2loglik(lat0);
    Precmat = as<arma::mat>(Q)-as<arma::mat>(g2loglik);
    Covmat = arma::inv(Precmat);
    latnew = Covmat * (as<arma::vec>(gloglik)-as<arma::mat>(g2loglik) * as<arma::vec>(lat0));
    NumericVector latnew2 = wrap(latnew);
    dist = sqrt(sum(pow((latnew2-lat0),2)));
    lat0 = latnew2;
    iter = iter + 1;
  }

  // Compute conditional latent vector and Covariance matrix

  double Sigma11 = Covmat(K-1,K-1);

  arma::mat Sigma12 = Covmat.row(K-1);
  Sigma12.shed_col(K-1);
  arma::mat Sigma21 = Sigma12;

  arma::mat Sigma22 = Covmat;
  Sigma22.shed_col(K-1);
  Sigma22.shed_row(K-1);

  arma::vec lat0minK = as<arma::vec>(lat0);
  lat0minK.shed_row(K-1);
  arma::mat latstarc = lat0minK + (1/Sigma11) * Sigma21.t() *
    (thetaconstr - lat0(K-1));

  arma::mat Covstarc = Sigma22 - (1/Sigma11) * (Sigma21.t() * Sigma12);
  NumericVector latstarc2 = wrap(latstarc);
  NumericVector latstarcc(dimlat);

  for(int i = 0; i < K; i++){
    latstarcc[i] = latstarc2[i];
  }
  latstarcc[K-1] = thetaconstr;
  for(int i = K; i < dimlat; i++){
    latstarcc[i] = latstarc2[i-1];
  }

  arma::cx_double logdetCovstarc = arma::log_det(Covstarc);

  return Rcpp::List::create(Named("latstar") = latstarcc,
                            Named("Covstar") = Covstarc,
                            Named("logdetCovstarc") = logdetCovstarc,
                            Named("iterations") = iter,
                            Named("distance") = dist);
}
