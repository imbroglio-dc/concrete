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

doTmleUpdate <- function(Estimates, SummEIC, Data, Censored, TargetEvents, TargetTimes, Events,
                         NumUpdateSteps, OneStepEps, NormPnEIC, Verbose) {
    Targets <- expand.grid("Time" = TargetTimes, "Event" = TargetEvents)
    EvalTimes <- attr(Estimates, "times")
    T.tilde <- Data[["Time"]]
    Delta <- Data[["Event"]]


## one-step tmle loop starts here ----
    for (step in 1:NumUpdateSteps) {
         if (Verbose)
            cat("starting step", step, "with update epsilon =", OneStepEps, "\n")

        ## Get updated hazards and EICs
        for (a in 1:length(Estimates)) {
            NuisanceWeight <- Estimates[[a]][["NuisanceWeight"]]
            GStar <- attr(Estimates[[a]][["PropScore"]], "g.star.intervention")
            Hazards <- Estimates[[a]][["OldHazards"]] <- Estimates[[a]][["Hazards"]]
            TotalSurv <- Estimates[[a]][["OldSurv"]] <- Estimates[[a]][["EvntFreeSurv"]]
            PnEIC <- Estimates[[a]][["SummEIC"]]

            NewHazards <- updateHazards(GStar, Hazards, TotalSurv, NuisanceWeight, Targets,
                                        Events, EvalTimes, T.tilde, Delta, PnEIC, OneStepEps)
            NewSurv <- apply(do.call(`+`, NewHazards, 2, function(haz) exp(-cumsum(haz))))

            GStar.obs <- attr(Estimates[[a]][["PropScore"]], "g.star.obs")
            IC.a <- getICs(GStar.obs, NewHazards, NewSurv, NuisanceWeight, Targets,
                           Events, T.tilde, Delta, EvalTimes, GComp = FALSE)

            Estimates[[a]][["Hazards"]] <- NewHazards
            Estimates[[a]][["EvntFreeSurv"]] <- NewSurv
            Estimates[[a]][["SummEIC"]] <- summarizeIC(IC.a)
        }

        ## Check for improvement
        OnestepStop <- SummEIC[, list("check" = abs(PnEIC) <= `seEIC/(root(n)log(n))`,
                                      "ratio" = abs(PnEIC) / `seEIC/(root(n)log(n))`),
                               by = c("Trt", "Time", "Event")]
        if (Verbose)
            print(OnestepStop[["ratio"]])
        if (all(sapply(OnestepStop[["check"]], isTRUE))) {

            return()
        }

        ## recalculate eic -----
        getEIC(EstEICs, Data, RegsOfInterest, Censored, TargetEvents,
               TargetTimes, Events, MinNuisance, GComp)
    }
}

updateHazards <- function(GStar, Hazards, TotalSurv, NuisanceWeight, Targets,
                          Events, EvalTimes, T.tilde, Delta, PnEIC, OneStepEps) {
    lapply(Hazards, function(haz.al) { # loop over L
        l <- attr(haz.al, "j")
        newhaz.al <- sapply(1:ncol(NuisanceWeight), function(i) {# loop over individuals
            Nuisance.i <- NuisanceWeight[, i]
            Surv.i <- TotalSurv[, i]
            Hazards.i <- lapply(Hazards, function(haz) haz[, i])
            Risks.i <- lapply(Hazards.i, function(haz.i) cumsum(Surv.i * haz.i))

            update.l <- rowSums(apply(Targets, 1, function(target) {
                j <- target[["Event"]]
                tau <- target[["Time"]]
                Nuisance.ik <- Nuisance.i * as.numeric(EvalTimes <= tau) ## 1(t \leq tau) * 1(t \leq t.tilde)
                F.j.tau <- Risks.i[[as.character(j)]][EvalTimes == tau]
                F.j.t <- Risks.i[[as.character(j)]]

                h.jk <- GStar[i] * Nuisance.ik * ((l == j) - (F.j.tau - F.j.t) / Surv.i)
                return(PnEIC[Time == tau & Event == j, PnEIC] * h.jk)
            }))
            return(haz.al[, i] * exp(update.l * OneStepEps / NormPnEIC))
        })
    })
}


