// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
//
// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"
using namespace Rcpp;
// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
// [[Rcpp::depends(RcppArmadillo)]]
//
// via the exports attribute we tell Rcpp to make this function
// available from R

// ============================================================================
// computeIC
//
// Computes the efficient influence curve for all (TargetEvent j, TargetTime tau)
// combinations in a single pass, without looping over hazard causes L.
//
// Math: for each (j, tau) and subject i,
//
//   IC[j,tau](i) = pm(i) - integral(i) + F_j(tau,i) - mean_i F_j(tau,i)
//
// where:
//   pm(i) = W[T~_i, i] * (1(Delta_i == j) - hFS[T~_i, i]) * 1(Delta_i != 0, T~_i <= tau)
//   integral(i) = sum_{t <= min(tau, T~_i)} W[t,i] * (h_j[t,i] - hFS[t,i] * hTotal[t,i])
//   hFS[t,i] = (F_j(tau,i) - F_j(t,i)) / S(t,i)   for t <= tau
//   W[t,i] = GStar_i * NuisanceWeight[t,i]
//   hTotal[t,i] = sum_l h_l[t,i]
//
// Key: the original code loops over all L cause hazards per (j,tau) pair.
// Here we collapse that loop using hTotal, achieving an L-fold speedup.
//
// Arguments:
//   Hazards        : named R list of L matrices, each T x N (cause-specific hazards)
//   TotalSurv      : T x N matrix, S(t|i)
//   NuisanceWeight : T x N matrix, 1 / (g(A|W) * S_C(t|A,W))
//   GStar          : N vector, g*(A|W) for each subject
//   T_tilde        : N vector, observed event/censoring times
//   Delta          : N integer vector, observed event type (0 = censored)
//   EvalTimes      : T vector, evaluation time grid (starts with 0, contains all T_tilde values)
//   TargetEvent    : J integer vector, event types to target
//   TargetTime     : K double vector, time points to target
//
// Returns a list with:
//   IC    : N x (J*K) matrix of IC values (column order: all K taus for j[1], then j[2], ...)
//   Time  : J*K vector of corresponding tau values
//   Event : J*K vector of corresponding j values
// [[Rcpp::export]]
List computeIC(List Hazards,
               arma::mat TotalSurv,
               arma::mat NuisanceWeight,
               arma::vec GStar,
               arma::vec T_tilde,
               arma::ivec Delta,
               arma::vec EvalTimes,
               arma::ivec TargetEvent,
               arma::vec TargetTime) {

    int T  = TotalSurv.n_rows;
    int N  = TotalSurv.n_cols;
    int J  = TargetEvent.n_elem;
    int K  = TargetTime.n_elem;
    int L  = Hazards.size();

    CharacterVector hazNames = Hazards.names();

    // ------------------------------------------------------------------
    // Precompute once: W, hTotal, AtRisk, tidx, W_obs
    // ------------------------------------------------------------------

    // W[t,i] = GStar[i] * NuisanceWeight[t,i]
    arma::mat W = NuisanceWeight;
    for (int i = 0; i < N; i++) W.col(i) *= GStar(i);

    // hTotal[t,i] = sum_l h_l[t,i]
    arma::mat hTotal(T, N, arma::fill::zeros);
    for (int l = 0; l < L; l++) hTotal += as<arma::mat>(Hazards[l]);

    // AtRisk[t,i] = 1(EvalTimes[t] <= T_tilde[i])  (stored as double for % ops)
    arma::mat AtRisk(T, N);
    for (int i = 0; i < N; i++)
        for (int t = 0; t < T; t++)
            AtRisk(t, i) = (EvalTimes(t) <= T_tilde(i)) ? 1.0 : 0.0;

    // tidx[i] = row index of T_tilde[i] in EvalTimes (exact match; guaranteed by construction)
    arma::uvec tidx(N);
    for (int i = 0; i < N; i++) {
        arma::uvec m = arma::find(EvalTimes == T_tilde(i), 1);
        tidx(i) = m.is_empty() ? (arma::uword)0 : m(0);
    }

    // W_obs[i] = W[tidx[i], i]  — weight at observed time, used for point mass
    arma::vec W_obs(N);
    for (int i = 0; i < N; i++) W_obs(i) = W(tidx(i), i);

    // ------------------------------------------------------------------
    // Main loops: J target events x K target times
    // ------------------------------------------------------------------

    int nComb = J * K;
    arma::mat IC_out(N, nComb);   // output: IC values, N x (J*K)
    arma::vec out_Time(nComb);
    arma::ivec out_Event(nComb);

    int comb = 0;

    for (int jj = 0; jj < J; jj++) {
        int j = TargetEvent(jj);

        // Find hazard matrix for event j
        std::string jname = std::to_string(j);
        int j_idx = -1;
        for (int l = 0; l < L; l++) {
            if (as<std::string>(hazNames[l]) == jname) { j_idx = l; break; }
        }
        if (j_idx < 0) stop("TargetEvent " + jname + " not found in Hazards");
        arma::mat haz_j = as<arma::mat>(Hazards[j_idx]);

        // F_j[t,i] = cumulative incidence = cumsum_t( h_j[t,i] * S[t,i] )
        arma::mat Fj = arma::cumsum(haz_j % TotalSurv, 0);  // cumsum along rows (dim 0)

        for (int kk = 0; kk < K; kk++) {
            double tau = TargetTime(kk);

            // Row index of tau in EvalTimes
            arma::uvec tau_m = arma::find(EvalTimes == tau, 1);
            if (tau_m.is_empty()) stop("TargetTime not found in EvalTimes");
            int t_tau = (int)tau_m(0);   // 0-indexed; rows 0..t_tau cover t <= tau

            // F_j_tau[i] = Fj[t_tau, i]
            arma::rowvec Fj_tau = Fj.row(t_tau);   // 1 x N

            // hFS[t,i] = (Fj_tau[i] - Fj[t,i]) / S[t,i]  for t = 0..t_tau
            // Shape: (t_tau+1) x N
            arma::mat Fj_sub = Fj.rows(0, t_tau);          // (t_tau+1) x N: F_j(t,i)
            arma::mat S_sub  = TotalSurv.rows(0, t_tau);   // (t_tau+1) x N: S(t,i)
            arma::mat hFS(t_tau + 1, N);
            hFS.each_row() = Fj_tau;     // broadcast: each row = F_j(tau, :)
            hFS -= Fj_sub;               // Fj_tau[i] - Fj[t,i]
            hFS /= S_sub;                // divide by S(t,i)

            // Integral: sum_{t <= min(tau, T~_i)} W[t,i] * (h_j[t,i] - hFS[t,i] * hTotal[t,i])
            arma::mat integrand = W.rows(0, t_tau) %
                (haz_j.rows(0, t_tau) - hFS % hTotal.rows(0, t_tau));
            // Mask by at-risk indicator (t <= T~_i) then sum over time
            arma::rowvec integral = arma::sum(integrand % AtRisk.rows(0, t_tau), 0);  // 1 x N

            // Point mass: only non-censored subjects with T~_i <= tau
            arma::rowvec pm(N, arma::fill::zeros);
            for (int i = 0; i < N; i++) {
                if (Delta(i) != 0 && T_tilde(i) <= tau) {
                    int ti = (int)tidx(i);   // row index of T~_i; guaranteed ti <= t_tau
                    double hfs_at_obs = hFS(ti, i);
                    double lj_indicator = (Delta(i) == j) ? 1.0 : 0.0;
                    pm(i) = W_obs(i) * (lj_indicator - hfs_at_obs);
                }
            }

            double mean_Fj_tau = arma::mean(Fj_tau);
            arma::rowvec IC = pm - integral + Fj_tau - mean_Fj_tau;

            IC_out.col(comb) = IC.t();
            out_Time(comb)   = tau;
            out_Event(comb)  = j;
            comb++;
        }
    }

    return List::create(Named("IC")    = IC_out,
                        Named("Time")  = out_Time,
                        Named("Event") = out_Event);
}

// ============================================================================
// updateHazardsCpp
//
// Performs one TMLE update step: for each cause-specific hazard h_l, computes
//   h_l_new[t,i] = h_l[t,i] * exp(StepSize * update_l[t,i])
//
// where:
//   update_l[t,i] = W[t,i] * (psum_l[t] - B[t,i])
//
// and B is shared across all l:
//   B[t,i] = (Fjsum[t,i] - Fpsum[t,i]) / S[t,i]
//
//   Fjsum[t,i] = sum_{j,k: TargetTime[k] >= EvalTimes[t]} F_j(tau_k, i) * PnEIC[j,k]
//   Fpsum[t,i] = sum_j F_j[t,i] * psum_j[t]
//   psum_j[t]  = sum_{k: tau_k >= EvalTimes[t]} PnEIC[j,k]   (reverse cumsum of PnEIC)
//   psum_l[t]  = psum_j[t] if l in TargetEvent (j index), else 0
//
// This collapses the original J*K inner loop per hazard into a single O(T*N*J)
// precomputation of B, then O(T*N) per hazard l.
//
// Arguments:
//   Hazards       : named R list of L matrices, each T x N
//   TotalSurv     : T x N
//   GStar         : N x 1 matrix (or vector; scales NuisanceWeight columns)
//   NuisanceWeight: T x N
//   EvalTimes     : T vector
//   TargetEvent   : character vector (J names matching Hazards names, e.g. "1", "2")
//   TargetTime    : K vector
//   PnEIC         : data.frame with columns Time, Event (character), PnEIC
//   OneStepEps    : step size numerator
//   NormPnEIC     : step size denominator (StepSize = OneStepEps / NormPnEIC)
// [[Rcpp::export]]
List updateHazardsCpp(List Hazards,
                      const arma::mat& TotalSurv,
                      const arma::mat& GStar,
                      const arma::mat& NuisanceWeight,
                      const arma::vec& EvalTimes,
                      const CharacterVector& TargetEvent,
                      const arma::vec& TargetTime,
                      const DataFrame& PnEIC,
                      double OneStepEps,
                      double NormPnEIC) {

    int T = TotalSurv.n_rows;
    int N = TotalSurv.n_cols;
    int J = TargetEvent.size();
    int K = TargetTime.n_elem;
    int L = Hazards.size();
    CharacterVector hazNames = Hazards.names();
    double StepSize = OneStepEps / NormPnEIC;

    // Extract PnEIC lookup table
    NumericVector pnTime  = PnEIC["Time"];
    CharacterVector pnEvt = PnEIC["Event"];
    NumericVector pnVals  = PnEIC["PnEIC"];
    int nPn = pnTime.size();

    // PnEIC_mat[jj, kk] = PnEIC for TargetEvent[jj] at TargetTime[kk]
    arma::mat PnEIC_mat(J, K, arma::fill::zeros);
    for (int jj = 0; jj < J; jj++) {
        std::string jname = as<std::string>(TargetEvent[jj]);
        for (int kk = 0; kk < K; kk++) {
            double tau = TargetTime(kk);
            bool found = false;
            for (int r = 0; r < nPn; r++) {
                if (std::abs(pnTime[r] - tau) < 1e-10 &&
                    std::string(pnEvt[r]) == jname) {
                    PnEIC_mat(jj, kk) = pnVals[r];
                    found = true;
                    break;
                }
            }
            if (!found) stop("PnEIC not found for Event=" + jname +
                             ", Time=" + std::to_string(tau));
        }
    }

    // ------------------------------------------------------------------
    // W[t,i] = GStar[i] * NuisanceWeight[t,i]
    // GStar is passed as N x 1 matrix; flatten to rowvec
    // ------------------------------------------------------------------
    arma::mat W = NuisanceWeight;
    arma::rowvec gs = arma::vectorise(GStar).t();  // 1 x N
    W.each_row() %= gs;

    // ------------------------------------------------------------------
    // Precompute F_j matrices for each target event j, and psum_j[t]
    // ------------------------------------------------------------------

    // Map TargetEvent names -> hazard indices
    std::vector<int> j_haz_idx(J, -1);
    for (int jj = 0; jj < J; jj++) {
        std::string jname = as<std::string>(TargetEvent[jj]);
        for (int l = 0; l < L; l++) {
            if (as<std::string>(hazNames[l]) == jname) { j_haz_idx[jj] = l; break; }
        }
        if (j_haz_idx[jj] < 0) stop("TargetEvent not found in Hazards");
    }

    // Fj_all[jj] = T x N cumulative incidence for event j
    std::vector<arma::mat> Fj_all(J);
    for (int jj = 0; jj < J; jj++) {
        arma::mat haz_j = as<arma::mat>(Hazards[j_haz_idx[jj]]);
        Fj_all[jj] = arma::cumsum(haz_j % TotalSurv, 0);
    }

    // psum[t, jj] = sum_{k: TargetTime[k] >= EvalTimes[t]} PnEIC[jj, k]
    // Computed for all t in 0..T-1
    arma::mat psum(T, J, arma::fill::zeros);
    for (int t = 0; t < T; t++) {
        for (int jj = 0; jj < J; jj++) {
            double s = 0.0;
            for (int kk = 0; kk < K; kk++) {
                if (EvalTimes(t) <= TargetTime(kk)) s += PnEIC_mat(jj, kk);
            }
            psum(t, jj) = s;
        }
    }

    // ------------------------------------------------------------------
    // Fjsum[t,i] = sum_{j,k: tau_k >= EvalTimes[t]} Fj(tau_k, i) * PnEIC[j,k]
    // Built by reverse scan: as t decreases, add Fj(tau_k, i)*P_{j,k} when
    // EvalTimes[t] first reaches tau_k from above.
    //
    // Equivalently: start running = 0; scan t from T-1 down to 0;
    // whenever EvalTimes[t] == tau_k, add Fj[t,i] * PnEIC[j,k] to running.
    // ------------------------------------------------------------------
    arma::mat Fjsum(T, N, arma::fill::zeros);
    {
        arma::rowvec running(N, arma::fill::zeros);
        for (int t = T - 1; t >= 0; t--) {
            for (int kk = 0; kk < K; kk++) {
                if (std::abs(EvalTimes(t) - TargetTime(kk)) < 1e-10) {
                    for (int jj = 0; jj < J; jj++) {
                        running += Fj_all[jj].row(t) * PnEIC_mat(jj, kk);
                    }
                    break;  // each EvalTimes[t] matches at most one TargetTime
                }
            }
            Fjsum.row(t) = running;
        }
    }

    // Fpsum[t,i] = sum_j Fj[t,i] * psum[t,j]
    arma::mat Fpsum(T, N, arma::fill::zeros);
    for (int jj = 0; jj < J; jj++) {
        for (int t = 0; t < T; t++) {
            Fpsum.row(t) += Fj_all[jj].row(t) * psum(t, jj);
        }
    }

    // B[t,i] = (Fjsum[t,i] - Fpsum[t,i]) / S[t,i]  — shared across all l
    arma::mat B = (Fjsum - Fpsum) / TotalSurv;

    // ------------------------------------------------------------------
    // Per-hazard update: update_l[t,i] = W[t,i] * (psum_l[t] - B[t,i])
    // ------------------------------------------------------------------
    List NewHazards(L);
    NewHazards.names() = hazNames;

    for (int l = 0; l < L; l++) {
        arma::mat haz_l = as<arma::mat>(Hazards[l]);
        std::string lname = as<std::string>(hazNames[l]);

        // psum_l_col[t] = psum[t, jj] if l is TargetEvent jj, else 0
        arma::vec psum_l(T, arma::fill::zeros);
        for (int jj = 0; jj < J; jj++) {
            if (as<std::string>(TargetEvent[jj]) == lname) {
                psum_l = psum.col(jj);
                break;
            }
        }

        // Broadcast psum_l across columns: update_l[t,i] = W[t,i] * (psum_l[t] - B[t,i])
        arma::mat psum_l_mat(T, N);
        psum_l_mat.each_col() = psum_l;

        arma::mat update_l = W % (psum_l_mat - B);
        NewHazards[l] = haz_l % arma::exp(update_l * StepSize);
    }

    return NewHazards;
}

// ============================================================================
// Legacy helpers — kept for reference/testing but no longer called by getIC
// ============================================================================

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
