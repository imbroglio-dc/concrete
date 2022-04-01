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
#' @param Censored boolean
#' @param TargetEvents numeric vector
#' @param TargetTimes numeric vector
#' @param Events numeric vector
#' @param NumUpdateSteps numeric
#' @param OneStepEps numeric
#' @param NormPnEIC numeric
#' @param Verbose boolean
#'
#'

doTmleUpdate <- function(Estimates, SummEIC, Data, Censored, TargetEvents, TargetTimes, Events,
                         NumUpdateSteps, OneStepEps, NormPnEIC, Verbose) {
    Time <- Event <- `seEIC/(root(n)log(n))` <- NULL
    Targets <- expand.grid("Time" = TargetTimes, "Event" = TargetEvents)
    EvalTimes <- attr(Estimates, "times")
    T.tilde <- Data[["Time"]]
    Delta <- Data[["Event"]]

    ## one-step tmle loop starts here ----
    for (step in 1:NumUpdateSteps) {
        if (Verbose)
            cat("starting step", step, "with update epsilon =", OneStepEps, "\n")

        ## Get updated hazards and EICs
        newEsts <- lapply(Estimates, function(est.a) {
            NewHazards <- updateHazard(GStar = attr(est.a[["PropScore"]], "g.star.intervention"),
                                       Hazards = est.a[["Hazards"]],
                                       TotalSurv = est.a[["EvntFreeSurv"]],
                                       NuisanceWeight = est.a[["NuisanceWeight"]],
                                       Events = Events, EvalTimes = EvalTimes, T.tilde = T.tilde,
                                       Delta = Delta, PnEIC = est.a[["SummEIC"]],
                                       NormPnEIC = NormPnEIC, OneStepEps = OneStepEps,
                                       TargetEvents = TargetEvents, TargetTimes = TargetTimes)
            NewSurv <- apply(do.call(`+`, NewHazards), 2, function(haz) exp(-cumsum(haz)))
            NewIC <- summarizeIC(getIC(GStar =  attr(est.a[["PropScore"]], "g.star.obs"),
                                       Hazards = NewHazards, TotalSurv = NewSurv,
                                       NuisanceWeight = est.a[["NuisanceWeight"]], Targets = Targets,
                                       Events = Events, T.tilde = T.tilde, Delta = Delta,
                                       EvalTimes = EvalTimes, GComp = FALSE))
            return(list("Hazards" = NewHazards, "EvntFreeSurv" = NewSurv, "SummEIC" = NewIC))
        })

        ## Check for improvement
        NewSummEIC <- do.call(rbind, lapply(seq_along(newEsts), function(a) {
            cbind("Trt" = names(newEsts)[a], newEsts[[a]][["SummEIC"]])}))
        NewNormPnEIC <- getNormPnEIC(NewSummEIC[Time %in% TargetTimes & Event %in% TargetEvents,
                                                PnEIC])
        if (NormPnEIC < NewNormPnEIC) {
            print("Update increased ||PnEIC||, halving the OneStepEps")
            OneStepEps <- 0.5 * OneStepEps
            next
        }

        for (a in seq_along(Estimates)) {
            Estimates[[a]][["Hazards"]] <- newEsts[[a]][["Hazards"]]
            Estimates[[a]][["EvntFreeSurv"]] <- newEsts[[a]][["EvntFreeSurv"]]
            Estimates[[a]][["SummEIC"]] <- newEsts[[a]][["SummEIC"]]
        }

        SummEIC <- NewSummEIC
        NormPnEIC <- NewNormPnEIC
        OneStepStop <- NewSummEIC[, list("check" = abs(PnEIC) <= `seEIC/(root(n)log(n))`,
                                         "ratio" = abs(PnEIC) / `seEIC/(root(n)log(n))`),
                                  by = c("Trt", "Time", "Event")]

        if (Verbose) print(OneStepStop[["ratio"]])
        if (all(sapply(OneStepStop[["check"]], isTRUE))) {
            return(Estimates)
        }
    }
    warning("TMLE has not converged by step", NumUpdateSteps, " - Results may not be reliable")
    return(Estimates)
}

updateHazard <- function(GStar, Hazards, TotalSurv, NuisanceWeight, Events, EvalTimes, T.tilde,
                         Delta, PnEIC, NormPnEIC, OneStepEps, TargetEvents, TargetTimes) {
    Time <- Event <- NULL
    if (min(TotalSurv) == 0)
        stop("max(TargetTimes) goes past the point where some people's survival probabilty = 0.",
             " This makes the clever covariate explode")
    lapply(Hazards, function(haz.al) { # loop over L
        l <- attr(haz.al, "j")
        update.l <- matrix(0, nrow = nrow(haz.al), ncol = ncol(haz.al))
        lapply(TargetEvents, function(j) {
            `F.j.t/S.t` <- apply(Hazards[[as.character(j)]] / TotalSurv, 2, cumsum)
            update.j <- sapply(TargetTimes, function(tau) {
                `F.j.tau/S.t` <- `F.j.t/S.t`[EvalTimes == tau, ]
                h.q <- (l == j) - t(apply(`F.j.t/S.t`, 1, function(r) `F.j.tau/S.t` - r))
                h.q[EvalTimes > tau, ] <- 0
                update <- h.q * NuisanceWeight %*% diag(GStar) * PnEIC[Time == tau & Event == j, PnEIC]
                update.l <<- update.l + update
                return(NULL)
            })
            return(NULL)
        })
        newhaz.al <- haz.al * exp(update.l * OneStepEps / NormPnEIC)

        # newhaz.al <- sapply(1:ncol(NuisanceWeight), function(i) {# loop over individuals
        #     if (GStar[i] == 0) {
        #         return(rep_len(0, nrow(NuisanceWeight)))
        #     }
        #     Nuisance.i <- NuisanceWeight[, i]
        #     Surv.i <- TotalSurv[, i]
        #     UpdateDir <- do.call(`+`, lapply(TargetEvents, function(j) {
        #         haz <- Hazards[[as.character(j)]]
        #         `F.j.t/S.t` <- cumsum(Surv.i * haz[, i] / Surv.i)
        #         `diffF/S` <- outer(`F.j.t/S.t`[EvalTimes %in% TargetTimes] , `F.j.t/S.t`, "-")
        #         clev.covs <- GStar[i] * apply(`diffF/S`, 1, function(r) ((l == attr(haz, "j")) - r) * Nuisance.i)
        #         update.dir <- rowSums(sapply(seq_along(TargetTimes), function(k) {
        #             tau <- TargetTimes[k]
        #             (EvalTimes <= tau) * clev.covs[, k] * PnEIC[Time == tau & Event == attr(haz, "j"), PnEIC]}))
        #         return(update.dir)}))
        #     return(haz.al[, i] * exp(UpdateDir * OneStepEps / NormPnEIC))
        # })
        attr(newhaz.al, "j") <- l
        return(newhaz.al)
    })
}


