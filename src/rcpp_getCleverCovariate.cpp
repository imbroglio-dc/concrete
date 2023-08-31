// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
//
// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"
// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]
//
// via the exports attribute we tell Rcpp to make this function
// available from R
//
// [[Rcpp::export]]
arma::dmat getCleverCovariate(arma::dvec GStar, 
                              arma::dmat NuisanceWeight, 
                              arma::dmat hFS, 
                              int LeqJ) {
    for (arma::uword i = 0; i < NuisanceWeight.n_cols; ++i) {
        NuisanceWeight.col(i) = GStar(i) * NuisanceWeight.col(i);
    }
    return NuisanceWeight % (LeqJ - hFS);
}

// [[Rcpp::export]]
arma::dmat getHazLS(arma::dvec T_Tilde,
                    arma::dvec EvalTimes,
                    arma::dmat HazL) {
    for (arma::uword i = 0; i < HazL.n_cols; ++i) {
        HazL.col(i) = (EvalTimes <= T_Tilde(i)) % HazL.col(i);
    }
    return HazL;
}

// Needs Debugging
// // [[Rcpp::export]]
// arma::cube updateHazardsCpp(arma::dvec J,
//                             arma::dvec TargetTimes,
//                             arma::dvec L, 
//                             arma::cube Hazards, 
//                             arma::dmat TotalSurv,
//                             arma::dvec EvalTimes,
//                             arma::dvec GStar, 
//                             arma::dmat NuisanceWeight,  
//                             arma::dvec PnEIC, 
//                             double StepSize) {
//     arma::cube UpdatedHaz(size(Hazards), arma::fill::zeros);
//     for (arma::uword l = 0; l < Hazards.n_slices; l++) {
//         for (arma::uword j = 0; j < J.size(); j++) {
//             int LeqJ = L(l) == J(j);
//             arma::dmat Fjt = Hazards.slice( arma::as_scalar(arma::find(L == J(j))) ) % TotalSurv;
//             Fjt = arma::cumsum(Fjt, 0);
//             for (arma::uword k = 0; k < TargetTimes.size(); k++) {
//                 arma::dmat hFS(Fjt.n_rows, Fjt.n_cols, arma::fill::zeros);
//                 arma::uword RowK = arma::as_scalar(arma::find(EvalTimes == TargetTimes(k)));
//                 hFS.rows(0, RowK) = arma::repmat(arma::dmat(Fjt.row(RowK)), RowK+1, 1);
//                 hFS.rows(0, RowK) = (hFS.rows(0, RowK) - Fjt.rows(0, RowK)) % TotalSurv.rows(0, RowK);
//                 arma::dmat hG(hFS.n_rows, hFS.n_cols, arma::fill::zeros);
//                 arma::vec TLeqTk(hG.n_rows, arma::fill::zeros);
//                 TLeqTk.head(RowK+1).ones();
//                 for (arma::uword i = 0; i < NuisanceWeight.n_cols; i++) {
//                     hG.col(i) = GStar(i) * TLeqTk % NuisanceWeight.col(i);
//                 }
//                 UpdatedHaz.slice(l) = UpdatedHaz.slice(l) + hG % (LeqJ - hFS) * PnEIC(j*TargetTimes.size() + k);
//             }
//         }
//     }
//     for (arma::uword m=0; m < UpdatedHaz.n_slices; m++) {
//         UpdatedHaz.slice(m) = Hazards.slice(m) % arma::exp(UpdatedHaz.slice(m) * StepSize);
//     }
//     return UpdatedHaz;
// }
