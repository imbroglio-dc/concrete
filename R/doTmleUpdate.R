# Pseudocode
# Variables:
#     N = number of observations
#     T = number of time points
#     J = number of target Events
#     L = number of Events
#     K = number of target time points
# Data:
#     For a in A:
#         pi(A | W)     : N x 1 vector
#         Sc(t | A, W)  : T x N matrix
# updated S(t | A, W)   : T x N matrix
# updated h(t | A, W)   : L (T x N) matrices
# Do:
#     For a in A:
#         for l in L:
#             for j in J and tau in K:
#                 for i in N:
#                     calculate       h_jk      = (T x (J x K)) matrix          # clever covariate
#                     fetch           PnEIC     = (J x K) x 1 vector            # emprcl mean EIC
#                     calculate  <PnEIC, h_jk>  = T x 1 vector                  # update direction
#           store     recalculate  h(t | A, W)  = T x 1 vector                  # updated hazard
#
#           store     recalculate S(t | A, W)   = (T x N) matrices
#
#         for j in J , tau in K: recalculate EIC :
#
#     undo or commit updated h() and F() based on PnEIC

#' Perform TMLE Targeting Step
#'
#' Iteratively updates the initial hazard estimates along the least favorable
#' submodel to solve the efficient influence curve equation. This is the core
#' "targeting" step that distinguishes TMLE from plug-in estimators.
#'
#' @param Estimates list; initial estimates from getInitialEstimate(), containing
#'   Hazards, EvntFreeSurv, PropScore, NuisanceWeight for each regime
#' @param SummEIC data.table; summarized EIC with PnEIC and stopping criterion
#' @param Data data.table; observed data with attribute column names
#' @param TargetEvent numeric vector; event types to target
#' @param TargetTime numeric vector; time points to target
#' @param MaxUpdateIter integer; maximum number of update iterations
#' @param OneStepEps numeric in (0, 1]; initial step size for updates
#' @param NormPnEIC numeric; initial ||PnEIC|| norm
#' @param Verbose logical; whether to print convergence diagnostics
#'
#' @return Updated Estimates list with:
#'   \itemize{
#'     \item Updated Hazards, EvntFreeSurv, SummEIC, IC for each regime
#'     \item Attribute TmleConverged: list(converged = TRUE/FALSE, step = integer)
#'     \item Attribute NormPnEICs: numeric vector of ||PnEIC|| trajectory
#'   }
#'
#' @details
#' The update uses a multiplicative perturbation of the hazards:
#'   h_new(t) = h_old(t) * exp(eps * sum_jk H_jk(t) * PnEIC_jk / ||PnEIC||)
#'
#' where H_jk is the clever covariate for event j at time k.
#'
#' Convergence is declared when: |PnEIC_jk| <= seEIC_jk / (sqrt(n) * log(n)) for all j, k
#'
#' If an update increases ||PnEIC||, the step size is halved and the iteration
#' repeats. This ensures monotonic decrease of ||PnEIC|| (up to numerical precision).
#'
#' @section Numerical Stability:
#' Several safeguards prevent numerical issues:
#' \itemize{
#'   \item Survival probabilities bounded away from 0 (minimum 1e-12)
#'   \item NA hazard values replaced with 0 (with warning)
#'   \item Step size halved if update increases ||PnEIC||
#'   \item Maximum iteration limit prevents infinite loops
#' }
#'
#' @importFrom nleqslv nleqslv
#' @keywords internal
doTmleUpdate <- function(Estimates, SummEIC, Data, TargetEvent, TargetTime,
                         MaxUpdateIter, OneStepEps, NormPnEIC, Verbose) {
  Time <- Event <- `seEIC/(sqrt(n)log(n))` <- PnEIC <- NULL
  EvalTimes <- attr(Estimates, "Times")
  T.tilde <- Data[[attr(Data, "EventTime")]]
  Delta <- Data[[attr(Data, "EventType")]]
  WorkingEps <- OneStepEps
  NormPnEICs <- NormPnEIC

  ## one-step tmle loop starts here ----
  StepNum <- 1
  IterNum <- 1
  # if (!Verbose) {
  #     Progress <- utils::txtProgressBar(min = 0, max = 5, initial = 0, style = 3)
  #     Q <- 0
  # }
  while (StepNum <= MaxUpdateIter & IterNum <= MaxUpdateIter * 2) {
    IterNum <- IterNum + 1
    if (Verbose) {
      cat("Starting step", StepNum, "with update epsilon =", WorkingEps, "\n")
      # } else {
      #     # add progress toward PnEIC cutoff ?
      #     MaxUpdateIterQuantile <- floor(quantile(1:MaxUpdateIter, probs = .2*(1:5)))
      #     if (StepNum %in% MaxUpdateIterQuantile) {
      #         Q <- Q + 1
      #         utils::setTxtProgressBar(Progress, value = Q)
      #     }
    }

    ## Get updated hazards and EICs
    newEsts <- lapply(Estimates, function(est.a) {
      ## deprecated in favor of the c++ version of the updateHazards function
      ## checked on pbc: max element-wise diff between R and C++ Hazards was less than .Machine$double.eps
      # NewHazards <- updateHazard(
      #   GStar = attr(est.a[["PropScore"]], "g.star.obs"),
      #   Hazards = est.a[["Hazards"]],
      #   TotalSurv = est.a[["EvntFreeSurv"]],
      #   NuisanceWeight = est.a[["NuisanceWeight"]],
      #   EvalTimes = EvalTimes, T.tilde = T.tilde,
      #   Delta = Delta, PnEIC = est.a[["SummEIC"]],
      #   NormPnEIC = NormPnEIC, OneStepEps = WorkingEps,
      #   TargetEvent = TargetEvent, TargetTime = TargetTime
      # )
      NewHazards <- updateHazardsCpp(
        Hazards = est.a[["Hazards"]],
        TotalSurv = est.a[["EvntFreeSurv"]],
        GStar = as.matrix(attr(est.a[["PropScore"]], "g.star.obs")),
        NuisanceWeight = est.a[["NuisanceWeight"]],
        EvalTimes = EvalTimes,
        TargetEvent = TargetEvent,
        TargetTime = TargetTime,
        PnEIC = est.a[["SummEIC"]],
        OneStepEps = WorkingEps,
        NormPnEIC = NormPnEIC
      )
      NewHazards <- lapply(NewHazards, function(hazards) {
        if (anyNA(hazards)) {
          warning("NA/NaN values detected in updated hazard estimates; setting to 0. ",
                  "This may indicate numerical instability. Consider increasing MinNuisance.",
                  call. = FALSE)
          hazards[is.na(hazards) | is.nan(hazards)] <- 0
        }
        return(hazards)
      })

      NewSurv <- apply(Reduce(`+`, NewHazards), 2, function(haz) exp(-cumsum(haz)))
      NewSurv[NewSurv < 1e-12 | is.na(NewSurv) | is.nan(NewSurv)] <- 1e-12
      NewIC <- getIC(
        GStar = attr(est.a[["PropScore"]], "g.star.obs"),
        Hazards = NewHazards, TotalSurv = NewSurv,
        NuisanceWeight = est.a[["NuisanceWeight"]],
        TargetEvent = TargetEvent, TargetTime = TargetTime,
        T.tilde = T.tilde, Delta = Delta,
        EvalTimes = EvalTimes, GComp = FALSE
      )
      return(list(
        "Hazards" = NewHazards, "EvntFreeSurv" = NewSurv,
        "SummEIC" = summarizeIC(NewIC), "IC" = NewIC
      ))
    })

    ## Check for improvement
    NewSummEIC <- do.call(rbind, lapply(seq_along(newEsts), function(a) {
      cbind("Trt" = names(newEsts)[a], newEsts[[a]][["SummEIC"]])
    }))
    NewNormPnEIC <- getNormPnEIC(NewSummEIC[Time %in% TargetTime & Event %in% TargetEvent, PnEIC])

    # TMLE update may produce NA if survival -> 0 or numerical instability
    if (anyNA(NewNormPnEIC)) {
      stop("NormPnEIC contains NA values; TMLE update may be numerically unstable. ",
              "Consider increasing MinNuisance or checking data for extreme values.")
    }

    if (NormPnEIC < NewNormPnEIC) {
      if (Verbose) cat("Update increased ||PnEIC||, halving OneStepEps\n")
      WorkingEps <- WorkingEps / 2
      next
    }
    StepNum <- StepNum + 1

    for (a in seq_along(Estimates)) {
      Estimates[[a]][["Hazards"]] <- newEsts[[a]][["Hazards"]]
      Estimates[[a]][["EvntFreeSurv"]] <- newEsts[[a]][["EvntFreeSurv"]]
      Estimates[[a]][["SummEIC"]] <- newEsts[[a]][["SummEIC"]]
      Estimates[[a]][["IC"]] <- newEsts[[a]][["IC"]]
    }

    SummEIC <- NewSummEIC
    NormPnEIC <- NewNormPnEIC
    NormPnEICs <- c(NormPnEICs, NewNormPnEIC)
    OneStepStop <- NewSummEIC[, list(
      "check" = abs(PnEIC) <= `seEIC/(sqrt(n)log(n))`,
      "ratio" = abs(PnEIC) / `seEIC/(sqrt(n)log(n))`
    ),
    by = c("Trt", "Time", "Event")
    ]

    if (Verbose) printOneStepDiagnostics(OneStepStop, NormPnEIC)

    if (all(sapply(OneStepStop[["check"]], isTRUE))) {
      attr(Estimates, "TmleConverged") <- list("converged" = TRUE, "step" = StepNum)
      attr(Estimates, "NormPnEICs") <- NormPnEICs
      return(Estimates)
    }
  }
  warning(
    "TMLE has not converged by step ", MaxUpdateIter, " - Estimates may not have ",
    "the desired asymptotic properties"
  )
  attr(Estimates, "TmleConverged") <- list("converged" = FALSE, "step" = StepNum)
  attr(Estimates, "NormPnEICs") <- NormPnEICs
  return(Estimates)
}

#' Update Hazards Along Least Favorable Submodel (R Implementation)
#'
#' Computes the one-step TMLE update of cause-specific hazards. This is the R
#' reference implementation; the C++ version (updateHazardsCpp) is used in production.
#'
#' @param GStar numeric vector; intervention probabilities g*(A|W) for each subject
#' @param Hazards list; cause-specific hazard matrices (times x subjects)
#' @param TotalSurv matrix; overall survival S(t) (times x subjects)
#' @param NuisanceWeight matrix; inverse of g-related weights (times x subjects)
#' @param EvalTimes numeric vector; time points for hazard evaluation
#' @param T.tilde numeric vector; observed event/censoring times
#' @param Delta numeric vector; observed event types (0 = censored)
#' @param PnEIC data.table; empirical mean EIC with columns Time, Event, PnEIC
#' @param NormPnEIC numeric; current ||PnEIC|| norm for scaling updates
#' @param OneStepEps numeric; step size multiplier
#' @param TargetEvent numeric vector; event types being targeted
#' @param TargetTime numeric vector; time points being targeted
#'
#' @return list of updated hazard matrices with same structure as input Hazards
#'
#' @details
#' For each cause l, the hazard update is:
#'   h_l,new(t) = h_l,old(t) * exp(eps / ||PnEIC|| * sum_\{j,tau\} H_\{l,j,tau\}(t) * PnEIC_\{j,tau\})
#'
#' where H_\{l,j,tau\} is the clever covariate measuring the sensitivity of the
#' risk F_j(tau) to perturbations in hazard h_l.
#'
#' @keywords internal
updateHazard <- function(GStar, Hazards, TotalSurv, NuisanceWeight, EvalTimes, T.tilde,
                         Delta, PnEIC, NormPnEIC, OneStepEps, TargetEvent, TargetTime, 
                        Iterative = FALSE) {
  eps <- Time <- Event <- NULL
  GStar <- as.numeric(unlist(GStar))
  if (min(TotalSurv) == 0) {
    stop("Survival probability reached exactly 0, causing undefined clever covariate. ",
         "This typically occurs when:\n",
         "  1. TargetTime extends beyond most observed events\n",
         "  2. Hazard estimates are unstable (try different model specification)\n",
         "  3. Sample size is insufficient for the target time horizon\n",
         "Consider targeting an earlier time point or checking data quality.",
         call. = FALSE)
  }
  if (Iterative) {
    warning("Iterative TMLE not yet implemented. Performing one-step TMLE instead.")
  }

  NewHazards <- lapply(Hazards, function(haz.al) { # loop over L
    l <- attr(haz.al, "j")

    update.l <-
      Reduce("+", x = lapply(TargetEvent, function(j) {
        F.j.t <- apply(Hazards[[as.character(j)]] * TotalSurv, 2, cumsum)
        Reduce("+", x = lapply(TargetTime, function(tau) {
          ClevCov <- h.FS <- matrix(0, nrow = nrow(F.j.t), ncol = ncol(F.j.t))
          h.FS[EvalTimes <= tau, ] <-
            (matrix(F.j.t[EvalTimes == tau, ],
              ncol = ncol(F.j.t),
              nrow = nrow(F.j.t[EvalTimes <= tau, ]),
              byrow = TRUE
            ) -
              F.j.t[EvalTimes <= tau, ]) /
              TotalSurv[EvalTimes <= tau, ]

          ClevCov[EvalTimes <= tau, ] <-
            getCleverCovariate(
              GStar = GStar,
              NuisanceWeight = NuisanceWeight[EvalTimes <= tau, ],
              hFS = h.FS[EvalTimes <= tau, ],
              LeqJ = as.integer(l == j)
            )

          return(ClevCov * PnEIC[Time == tau & Event == j, PnEIC])
        }))
      }))
    newhaz.al <- haz.al * exp(update.l * OneStepEps / NormPnEIC)
    attr(newhaz.al, "j") <- l
    return(newhaz.al)
  })

  return(NewHazards)
}

# getFluctPnEIC <- function(GStar, Hazards, TotalSurv, NuisanceWeight, TargetEvent, TargetTime,
#                           T.tilde, Delta, EvalTimes, GComp, l = NULL, fluct.eps = NULL) {
#     IC <- F.j.tau <- eps <- NULL
#     Target <- expand.grid("Time" = TargetTime, "Event" = TargetEvent)
#     IC.a <- do.call(rbind, lapply(1:ncol(NuisanceWeight), function(i) {
#         Nuisance.i <- NuisanceWeight[, i]
#         Surv.i <- TotalSurv[, i]
#         Hazards.i <- lapply(Hazards, function(haz) haz[, i])
#         Risks.i <- lapply(Hazards.i, function(haz.i) cumsum(Surv.i * haz.i))
#
#         if (GStar[i] == 0) # 1(A != a*)
#             return(cbind("ID" = i, Target,  "IC" = 0,
#                          "F.j.tau" = apply(Target,  1, function(target) {
#                              tau <- target[["Time"]]
#                              j <- target[["Event"]]
#                              return(Risks.i[[as.character(j)]][EvalTimes == tau])
#                          })))
#
#         IC.jk <- t(apply(Target,  1, function(target) {
#             j <- target[["Event"]]
#             tau <- target[["Time"]]
#             t.tilde <- T.tilde[i]
#             TimeIndices.ik <- EvalTimes <= min(tau, t.tilde) ## 1(t \leq tau) * 1(t \leq t.tilde)
#             F.j.tau <- Risks.i[[as.character(j)]][EvalTimes == tau]
#
#             s.ik <- EvalTimes[TimeIndices.ik]
#             Nuisance.ik <- Nuisance.i[TimeIndices.ik]
#             Surv.ik <- Surv.i[TimeIndices.ik]
#             haz.J.ik <- lapply(Hazards.i, function(r) r[TimeIndices.ik])
#             F.j.t <- Risks.i[[as.character(j)]][TimeIndices.ik]
#
#             h.jk <- GStar[i] * Nuisance.ik * ((l == j) - (F.j.tau - F.j.t) / Surv.ik)
#             IC.ljk <- sum(h.jk * ((s.ik == t.tilde) * (Delta[i] == l) - eps * haz.J.ik[[as.character(l)]]))
#             return(c("IC" = IC.ljk, "F.j.tau" = F.j.tau))
#         }))
#         IC.jk <- cbind("ID" = i, Target, IC.jk)
#         return(IC.jk)
#     }))
#     ## the second EIC component ( ... + F_j(tau | a, L) - Psi )
#     IC.a <- as.data.table(IC.a)[, list(mean(IC))]
#     return(IC.a)
# }

#' Print TMLE Convergence Diagnostics
#'
#' Displays information about the current state of TMLE convergence, including
#' the targets with the highest |PnEIC| / stopping_criterion ratios.
#'
#' @param OneStepStop data.table with columns Trt, Time, Event, check (logical),
#'   and ratio (|PnEIC| / stopping_criterion)
#' @param NormPnEIC numeric; current Euclidean norm of PnEIC vector
#'
#' @return NULL (called for side effect of printing)
#'
#' @details
#' The ratio column indicates how far each target is from convergence:
#' - ratio < 1: Target has converged

#' - ratio = 1: Exactly at stopping criterion
#' - ratio > 1: Target has not yet converged; higher values indicate
#'   more targeting needed
#'
#' @keywords internal
printOneStepDiagnostics <- function(OneStepStop, NormPnEIC) {
  ratio <- NULL
  Worst <- OneStepStop[, !"check"][, ratio := round(ratio, 2)][order(ratio, decreasing = TRUE)]
  print(Worst[1:min(nrow(Worst), 3), ])
  cat("Norm PnEIC = ", NormPnEIC, "\n", sep = "")
}
