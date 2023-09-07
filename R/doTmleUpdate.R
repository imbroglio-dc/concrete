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

#' Title
#'
#' @param Estimates list
#' @param SummEIC data.table
#' @param Data data.table
#' @param TargetEvent numeric vector
#' @param TargetTime numeric vector
#' @param MaxUpdateIter numeric
#' @param OneStepEps numeric
#' @param NormPnEIC numeric
#' @param Verbose boolean
#'
#' @importFrom nleqslv nleqslv

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
            NewHazards <- updateHazard(GStar = attr(est.a[["PropScore"]], "g.star.obs"),
                                       Hazards = est.a[["Hazards"]],
                                       TotalSurv = est.a[["EvntFreeSurv"]],
                                       NuisanceWeight = est.a[["NuisanceWeight"]],
                                       EvalTimes = EvalTimes, T.tilde = T.tilde,
                                       Delta = Delta, PnEIC = est.a[["SummEIC"]],
                                       NormPnEIC = NormPnEIC, OneStepEps = WorkingEps,
                                       TargetEvent = TargetEvent, TargetTime = TargetTime)
            NewHazards <- lapply(NewHazards, function(hazards) {
                if (anyNA(hazards))
                    hazards[is.na(hazards) | is.nan(hazards)] <-  0
                return(hazards)
            })
            NewSurv <- apply(Reduce(`+`, NewHazards), 2, function(haz) exp(-cumsum(haz)))
            NewSurv[NewSurv < 1e-12 | is.na(NewSurv) | is.nan(NewSurv)] <- 1e-12
            NewIC <- getIC(GStar =  attr(est.a[["PropScore"]], "g.star.obs"),
                           Hazards = NewHazards, TotalSurv = NewSurv,
                           NuisanceWeight = est.a[["NuisanceWeight"]],
                           TargetEvent = TargetEvent, TargetTime = TargetTime,
                           T.tilde = T.tilde, Delta = Delta,
                           EvalTimes = EvalTimes, GComp = FALSE)
            return(list("Hazards" = NewHazards, "EvntFreeSurv" = NewSurv, 
                        "SummEIC" = summarizeIC(NewIC), "IC" = NewIC))
        })
        
        ## Check for improvement
        NewSummEIC <- do.call(rbind, lapply(seq_along(newEsts), function(a) {
            cbind("Trt" = names(newEsts)[a], newEsts[[a]][["SummEIC"]])}))
        NewNormPnEIC <- getNormPnEIC(NewSummEIC[Time %in% TargetTime & Event %in% TargetEvent, PnEIC])
        
        # TMLE update breaking because Survival -> 0?
        if(anyNA(NewNormPnEIC)) browser()
        
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
        OneStepStop <- NewSummEIC[, list("check" = abs(PnEIC) <= `seEIC/(sqrt(n)log(n))`,
                                         "ratio" = abs(PnEIC) / `seEIC/(sqrt(n)log(n))`),
                                  by = c("Trt", "Time", "Event")]
        
        if (Verbose) printOneStepDiagnostics(OneStepStop, NormPnEIC)
        
        if (all(sapply(OneStepStop[["check"]], isTRUE))) {
            attr(Estimates, "TmleConverged") <- list("converged" = TRUE, "step" = StepNum)
            attr(Estimates, "NormPnEICs") <- NormPnEICs
            return(Estimates)
        }
    }
    warning("TMLE has not converged by step ", MaxUpdateIter, " - Estimates may not have ",
            "the desired asymptotic properties")
    attr(Estimates, "TmleConverged") <- list("converged" = FALSE, "step" = StepNum)
    attr(Estimates, "NormPnEICs") <- NormPnEICs
    return(Estimates)
}

updateHazard <- function(GStar, Hazards, TotalSurv, NuisanceWeight, EvalTimes, T.tilde,
                         Delta, PnEIC, NormPnEIC, OneStepEps, TargetEvent, TargetTime) {
    eps <- Time <- Event <- NULL
    Iterative <- FALSE
    GStar <- as.numeric(unlist(GStar))
    if (min(TotalSurv) == 0)
        stop("People's survival probabilty -> 0 makes the clever covariate explode.")
    if (Iterative) {
        warning("Iterative TMLE not yet implemented. Performing one-step TMLE instead.")
    } 
    # tmp <- microbenchmark(
    #     "cpp" = {
    #         NewHaz <- updateHazardsCpp(J = TargetEvent,
    #                                    TargetTimes = TargetTime, 
    #                                    L = as.numeric(names(Hazards)), 
    #                                    Hazards = array(unlist(Hazards), dim = c(dim(Hazards[[1]]), length(Hazards))), 
    #                                    TotalSurv = TotalSurv, 
    #                                    EvalTimes = EvalTimes, 
    #                                    GStar = GStar, 
    #                                    NuisanceWeight = NuisanceWeight, 
    #                                    PnEIC = unlist(PnEIC[Event > 0, ][order(Event, Time), PnEIC]), 
    #                                    StepSize = OneStepEps / NormPnEIC)
    #         NewHaz <- lapply(seq(dim(NewHaz)[3]), function(l) NewHaz[, , l])}, 
    #     "apply" = {
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
                                byrow = TRUE) - 
                             F.j.t[EvalTimes <= tau, ]) / 
                        TotalSurv[EvalTimes <= tau, ]
                    
                    ClevCov[EvalTimes <= tau, ] <- 
                        getCleverCovariate(GStar = GStar, 
                                           NuisanceWeight = NuisanceWeight[EvalTimes <= tau, ], 
                                           hFS = h.FS[EvalTimes <= tau, ], 
                                           LeqJ = as.integer(l == j))
                    
                    return(ClevCov * PnEIC[Time == tau & Event == j, PnEIC])
                }))
            }))
        newhaz.al <- haz.al * exp(update.l * OneStepEps / NormPnEIC)
        attr(newhaz.al, "j") <- l
        return(newhaz.al)
    })
    # })
    
    #     eps.l <- nleqslv(0.01, function(eps) getFluctPnEIC(GStar = GStar, Hazards = Hazards,
    #                                                        TotalSurv = TotalSurv,
    #                                                        NuisanceWeight = NuisanceWeight,
    #                                                        TargetEvent = TargetEvent,
    #                                                        TargetTime = TargetTime, T.tilde = T.tilde,
    #                                                        Delta = Delta, EvalTimes = EvalTimes,
    #                                                        GComp = FALSE, l = attr(haz.al, "j"),
    #                                                        fluct.eps = eps))$x
    #     lapply(TargetEvent, function(j) {
    #         `F.j.t` <- apply(Hazards[[as.character(j)]] * TotalSurv, 2, cumsum)
    #         eps.lj <- nleqslv(0.01, function(eps) getFluctPnEIC(GStar = GStar, Hazards = Hazards,
    #                                                             TotalSurv = TotalSurv,
    #                                                             NuisanceWeight = NuisanceWeight,
    #                                                             TargetEvent = j,
    #                                                             TargetTime = TargetTime, T.tilde = T.tilde,
    #                                                             Delta = Delta, EvalTimes = EvalTimes,
    #                                                             GComp = FALSE, l = attr(haz.al, "j"),
    #                                                             fluct.eps = eps))$x
    #         update.j <- sapply(TargetTime, function(tau) {
    #             eps.ljk <- nleqslv(0.01, function(eps) getFluctPnEIC(GStar = GStar, Hazards = Hazards,
    #                                                                  TotalSurv = TotalSurv,
    #                                                                  NuisanceWeight = NuisanceWeight,
    #                                                                  TargetEvent = j,
    #                                                                  TargetTime = tau, T.tilde = T.tilde,
    #                                                                  Delta = Delta, EvalTimes = EvalTimes,
    #                                                                  GComp = FALSE, l = attr(haz.al, "j"),
    #                                                                  fluct.eps = eps))$x
    #             `F.j.tau` <- `F.j.t`[EvalTimes == tau, ]
    #             h.q <- (l == j) - t(apply(`F.j.t`, 1, function(r) `F.j.tau` - r)) / TotalSurv
    #             h.q[EvalTimes > tau, ] <- 0
    #             update <- h.q * NuisanceWeight %*% diag(GStar) * eps.ljk
    #             update.l <<- update.l + update
    #             return(NULL)
    #         })
    #         return(NULL)
    #     })
    # }
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

printOneStepDiagnostics <- function(OneStepStop, NormPnEIC) {
    ratio <- NULL
    Worst <- OneStepStop[, !"check"][, ratio := round(ratio, 2)][order(ratio, decreasing = TRUE)]
    print(Worst[1:min(nrow(Worst), 3), ])
    cat("Norm PnEIC = ", NormPnEIC, "\n", sep = "")
}