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
// [[Rcpp::export]]
List updateHazardsCpp(List Hazards,
                      const arma::mat& TotalSurv,       // T x N
                      const arma::vec& GStar,           // N
                      const arma::mat& NuisanceWeight,  // T x N
                      const arma::vec& EvalTimes,       // T
                      const CharacterVector& TargetEvent, // J names
                      const arma::vec& TargetTime,      // K
                      const DataFrame& PnEIC,           // columns: Time, Event, PnEIC
                      double OneStepEps,
                      double NormPnEIC) {

    int L = Hazards.size();
    int T = TotalSurv.n_rows;
    int N = TotalSurv.n_cols;

    List NewHazards(L);
    NewHazards.names() = Hazards.names();
    CharacterVector hazardNames = Hazards.names();

    // Extract PnEIC columns
    NumericVector pnTime = PnEIC["Time"];
    CharacterVector pnEvent = PnEIC["Event"];
    NumericVector pnVals = PnEIC["PnEIC"];
    int nPn = pnTime.size();

    for (int l = 0; l < L; l++) {
        arma::mat haz_al = as<arma::mat>(Hazards[l]); // T x N
        arma::mat update_l(T, N, arma::fill::zeros);
        std::string lname = as<std::string>(hazardNames[l]);

        for (int j = 0; j < TargetEvent.size(); j++) {
            std::string jname = as<std::string>(TargetEvent[j]);

            // Find hazard by name
            int idx = -1;
            for (int m = 0; m < hazardNames.size(); m++) {
                if (hazardNames[m] == jname) { idx = m; break; }
            }
            if (idx == -1) stop("TargetEvent not found in Hazards");

            arma::mat haz_j = as<arma::mat>(Hazards[idx]);

            // Compute cumulative hazards * survival (Fjt)
            arma::mat Fjt(T, N);
            for (int n = 0; n < N; n++) {
                Fjt.col(n) = arma::cumsum(haz_j.col(n) % TotalSurv.col(n));
            }

            // Pre-scale NuisanceWeight once per event
            arma::mat scaledNW = NuisanceWeight;
            scaledNW.each_row() %= GStar.t();

            // Loop over target times (serial)
            for (int k = 0; k < TargetTime.n_elem; k++) {
                double tau = TargetTime[k];

                // Find corresponding PnEIC value(s) for this (Time, Event)
                double pn = 0.0;
                for (int r = 0; r < nPn; r++) {
                    if (pnTime[r] == tau && std::string(pnEvent[r]) == jname) {
                        pn = pnVals[r];
                        break;
                    }
                }

                // If no PnEIC found for this (tau,j), stop
                if (pn == 0.0) ) stop("TargetEvent not found in Hazards");

                // Find index of tau in EvalTimes
                arma::uvec tau_idx = arma::find(EvalTimes == tau, 1);
                if (tau_idx.n_elem == 0) continue;
                int t_tau = tau_idx(0);

                // F_tau is row of cumulative hazards at tau
                arma::rowvec F_tau = Fjt.row(t_tau);

                // Compute hFS efficiently using broadcasting
                arma::mat hFS = Fjt;   
                hFS.each_row() = F_tau - hFS.each_row();  // broadcast subtraction
                hFS /= TotalSurv;                          // elementwise division
                hFS.rows(arma::find(EvalTimes > tau)).zeros(); // zero rows after tau

                // clever covariate calculation
                arma::mat ClevCov = scaledNW % ( (lname == jname ? 1.0 : 0.0) - hFS );

                // Update hazard increment (in-place multiplication)
                update_l += ClevCov * pn;
            }
        }

        // Apply one-step TMLE update to hazard
        arma::mat newhaz_al = haz_al % arma::exp(update_l * (OneStepEps / NormPnEIC));
        NewHazards[l] = newhaz_al;
    }

    return NewHazards;
}
