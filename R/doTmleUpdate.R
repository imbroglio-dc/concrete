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
    EvalTimes <- attr(Estimates, "times")
    T.tilde <- Data[[attr(Data, "EventTime")]]
    Delta <- Data[[attr(Data, "EventType")]]
    WorkingEps <- OneStepEps

    ## one-step tmle loop starts here ----
    StepNum <- 1
    IterNum <- 1
    while (StepNum <= MaxUpdateIter & IterNum <= MaxUpdateIter * 2) {
        IterNum <- IterNum + 1
        if (Verbose)
            cat("starting step", StepNum, "with update epsilon =", WorkingEps, "\n")

        ## Get updated hazards and EICs
        newEsts <- lapply(Estimates, function(est.a) {
            NewHazards <- updateHazard(GStar = attr(est.a[["PropScore"]], "g.star.intervention"),
                                       Hazards = est.a[["Hazards"]],
                                       TotalSurv = est.a[["EvntFreeSurv"]],
                                       NuisanceWeight = est.a[["NuisanceWeight"]],
                                       EvalTimes = EvalTimes, T.tilde = T.tilde,
                                       Delta = Delta, PnEIC = est.a[["SummEIC"]],
                                       NormPnEIC = NormPnEIC, OneStepEps = WorkingEps,
                                       TargetEvent = TargetEvent, TargetTime = TargetTime)
            NewSurv <- apply(Reduce(`+`, NewHazards), 2, function(haz) exp(-cumsum(haz)))
            NewIC <- summarizeIC(
                getIC(GStar =  attr(est.a[["PropScore"]], "g.star.obs"),
                      Hazards = NewHazards, TotalSurv = NewSurv,
                      NuisanceWeight = est.a[["NuisanceWeight"]],
                      TargetEvent = TargetEvent, TargetTime = TargetTime,
                      T.tilde = T.tilde, Delta = Delta,
                      EvalTimes = EvalTimes, GComp = FALSE)
            )
            return(list("Hazards" = NewHazards, "EvntFreeSurv" = NewSurv, "SummEIC" = NewIC))
        })

        ## Check for improvement
        NewSummEIC <- do.call(rbind, lapply(seq_along(newEsts), function(a) {
            cbind("Trt" = names(newEsts)[a], newEsts[[a]][["SummEIC"]])}))
        NewNormPnEIC <- getNormPnEIC(NewSummEIC[Time %in% TargetTime & Event %in% TargetEvent, PnEIC])

        if (NormPnEIC < NewNormPnEIC) {
            if (Verbose) {
                cat("Update increased ||PnEIC||, halving the OneStepEps\n")
            }
            WorkingEps <- WorkingEps / 2
            next
        }
        StepNum <- StepNum + 1

        for (a in seq_along(Estimates)) {
            Estimates[[a]][["Hazards"]] <- newEsts[[a]][["Hazards"]]
            Estimates[[a]][["EvntFreeSurv"]] <- newEsts[[a]][["EvntFreeSurv"]]
            Estimates[[a]][["SummEIC"]] <- newEsts[[a]][["SummEIC"]]
        }

        SummEIC <- NewSummEIC
        NormPnEIC <- NewNormPnEIC
        OneStepStop <- NewSummEIC[, list("check" = abs(PnEIC) <= `seEIC/(sqrt(n)log(n))`,
                                         "ratio" = abs(PnEIC) / `seEIC/(sqrt(n)log(n))`),
                                  by = c("Trt", "Time", "Event")]
        
        if (Verbose)  printOneStepDiagnostics(OneStepStop)
        
        if (all(sapply(OneStepStop[["check"]], isTRUE))) {
            attr(Estimates, "TmleConverged") <- list("converged" = TRUE, "step" = StepNum)
            return(Estimates)
        }
    }
    warning("TMLE has not converged by step ", MaxUpdateIter, " - Estimates may not have ",
            "the desired asymptotic properties")
    attr(Estimates, "TmleConverged") <- list("converged" = FALSE, "step" = StepNum)
    return(Estimates)
}

updateHazard <- function(GStar, Hazards, TotalSurv, NuisanceWeight, EvalTimes, T.tilde,
                         Delta, PnEIC, NormPnEIC, OneStepEps, TargetEvent, TargetTime) {
    eps <- Time <- Event <- NULL
    Iterative <- FALSE
    if (min(TotalSurv) == 0)
        stop("max(TargetTime) when people's survival probabilty -> 0. ",
             "This makes the clever covariate explode.")
    lapply(Hazards, function(haz.al) { # loop over L
        l <- attr(haz.al, "j")
        update.l <- matrix(0, nrow = nrow(haz.al), ncol = ncol(haz.al))
        lapply(TargetEvent, function(j) {
            `F.j.t` <- apply(Hazards[[as.character(j)]] * TotalSurv, 2, cumsum)
            update.j <- sapply(TargetTime, function(tau) {
                `F.j.tau` <- `F.j.t`[EvalTimes == tau, ]
                h.q <- (l == j) - t(apply(`F.j.t`, 1, function(r) `F.j.tau` - r)) / TotalSurv
                h.q[EvalTimes > tau, ] <- 0
                update <- h.q * NuisanceWeight %*% diag(GStar) * PnEIC[Time == tau & Event == j, PnEIC]
                update.l <<- update.l + update
                return(NULL)
            })
            return(NULL)
        })
        newhaz.al <- haz.al * exp(update.l * OneStepEps / NormPnEIC)

        if (Iterative) {
            eps.l <- nleqslv(0.01, function(eps) getFluctPnEIC(GStar = GStar, Hazards = Hazards,
                                                               TotalSurv = TotalSurv,
                                                               NuisanceWeight = NuisanceWeight,
                                                               TargetEvent = TargetEvent,
                                                               TargetTime = TargetTime, T.tilde = T.tilde,
                                                               Delta = Delta, EvalTimes = EvalTimes,
                                                               GComp = FALSE, l = attr(haz.al, "j"),
                                                               fluct.eps = eps))$x
            lapply(TargetEvent, function(j) {
                `F.j.t` <- apply(Hazards[[as.character(j)]] * TotalSurv, 2, cumsum)
                eps.lj <- nleqslv(0.01, function(eps) getFluctPnEIC(GStar = GStar, Hazards = Hazards,
                                                                    TotalSurv = TotalSurv,
                                                                    NuisanceWeight = NuisanceWeight,
                                                                    TargetEvent = j,
                                                                    TargetTime = TargetTime, T.tilde = T.tilde,
                                                                    Delta = Delta, EvalTimes = EvalTimes,
                                                                    GComp = FALSE, l = attr(haz.al, "j"),
                                                                    fluct.eps = eps))$x
                update.j <- sapply(TargetTime, function(tau) {
                    eps.ljk <- nleqslv(0.01, function(eps) getFluctPnEIC(GStar = GStar, Hazards = Hazards,
                                                                         TotalSurv = TotalSurv,
                                                                         NuisanceWeight = NuisanceWeight,
                                                                         TargetEvent = j,
                                                                         TargetTime = tau, T.tilde = T.tilde,
                                                                         Delta = Delta, EvalTimes = EvalTimes,
                                                                         GComp = FALSE, l = attr(haz.al, "j"),
                                                                         fluct.eps = eps))$x
                    `F.j.tau` <- `F.j.t`[EvalTimes == tau, ]
                    h.q <- (l == j) - t(apply(`F.j.t`, 1, function(r) `F.j.tau` - r)) / TotalSurv
                    h.q[EvalTimes > tau, ] <- 0
                    update <- h.q * NuisanceWeight %*% diag(GStar) * eps.ljk
                    update.l <<- update.l + update
                    return(NULL)
                })
                return(NULL)
            })
        }

        attr(newhaz.al, "j") <- l
        return(newhaz.al)

        # newhaz.al <- sapply(1:ncol(NuisanceWeight), function(i) {# loop over individuals
        #     if (GStar[i] == 0) {
        #         return(rep_len(0, nrow(NuisanceWeight)))
        #     }
        #     Nuisance.i <- NuisanceWeight[, i]
        #     Surv.i <- TotalSurv[, i]
        #     UpdateDir <- Reduce(`+`, lapply(TargetEvent, function(j) {
        #         haz <- Hazards[[as.character(j)]]
        #         `F.j.t/S.t` <- cumsum(Surv.i * haz[, i]) / Surv.i
        #         `diffF/S` <- outer(`F.j.t/S.t`[EvalTimes %in% TargetTime] , `F.j.t/S.t`, "-")
        #         clev.covs <- GStar[i] * apply(`diffF/S`, 1, function(r) ((l == attr(haz, "j")) - r) * Nuisance.i)
        #         update.dir <- rowSums(sapply(seq_along(TargetTime), function(k) {
        #             tau <- TargetTime[k]
        #             (EvalTimes <= tau) * clev.covs[, k] * PnEIC[Time == tau & Event == attr(haz, "j"), PnEIC]}))
        #         return(update.dir)}))
        #     return(haz.al[, i] * exp(UpdateDir * OneStepEps / NormPnEIC))
        # })
    })
}

getFluctPnEIC <- function(GStar, Hazards, TotalSurv, NuisanceWeight, TargetEvent, TargetTime,
                          T.tilde, Delta, EvalTimes, GComp, l = NULL, fluct.eps = NULL) {
    IC <- F.j.tau <- eps <- NULL
    Target <- expand.grid("Time" = TargetTime, "Event" = TargetEvent)
    IC.a <- do.call(rbind, lapply(1:ncol(NuisanceWeight), function(i) {
        Nuisance.i <- NuisanceWeight[, i]
        Surv.i <- TotalSurv[, i]
        Hazards.i <- lapply(Hazards, function(haz) haz[, i])
        Risks.i <- lapply(Hazards.i, function(haz.i) cumsum(Surv.i * haz.i))

        if (GStar[i] == 0) # 1(A != a*)
            return(cbind("ID" = i, Target,  "IC" = 0,
                         "F.j.tau" = apply(Target,  1, function(target) {
                             tau <- target[["Time"]]
                             j <- target[["Event"]]
                             return(Risks.i[[as.character(j)]][EvalTimes == tau])
                         })))

        IC.jk <- t(apply(Target,  1, function(target) {
            j <- target[["Event"]]
            tau <- target[["Time"]]
            t.tilde <- T.tilde[i]
            TimeIndices.ik <- EvalTimes <= min(tau, t.tilde) ## 1(t \leq tau) * 1(t \leq t.tilde)
            F.j.tau <- Risks.i[[as.character(j)]][EvalTimes == tau]

            s.ik <- EvalTimes[TimeIndices.ik]
            Nuisance.ik <- Nuisance.i[TimeIndices.ik]
            Surv.ik <- Surv.i[TimeIndices.ik]
            haz.J.ik <- lapply(Hazards.i, function(r) r[TimeIndices.ik])
            F.j.t <- Risks.i[[as.character(j)]][TimeIndices.ik]

            h.jk <- GStar[i] * Nuisance.ik * ((l == j) - (F.j.tau - F.j.t) / Surv.ik)
            IC.ljk <- sum(h.jk * ((s.ik == t.tilde) * (Delta[i] == l) - eps * haz.J.ik[[as.character(l)]]))
            return(c("IC" = IC.ljk, "F.j.tau" = F.j.tau))
        }))
        IC.jk <- cbind("ID" = i, Target, IC.jk)
        return(IC.jk)
    }))
    ## the second EIC component ( ... + F_j(tau | a, L) - Psi )
    IC.a <- as.data.table(IC.a)[, list(mean(IC))]
    return(IC.a)
}

printOneStepDiagnostics <- function(OneStepStop) {
    ratio <- NULL
    Verb <- OneStepStop[, !"check"][, ratio := round(ratio, 2)][order(ratio, decreasing = TRUE)]
    print(Verb[1:min(nrow(Verb), 3), ])
}