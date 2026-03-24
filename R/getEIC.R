#' Compute Efficient Influence Curves for TMLE
#'
#' Calculates the efficient influence curve (EIC) for each intervention regime,
#' target event, and target time combination. The EIC is used for both variance
#' estimation and the TMLE targeting step.
#'
#' @param Estimates list of estimation results for each intervention regime,
#'   containing Hazards, EvntFreeSurv, PropScore, and NuisanceWeight
#' @param Data data.table of observed data with EventTime and EventType attributes
#' @param Regime list of intervention regime specifications (currently unused but
#'   passed for consistency)
#' @param TargetEvent numeric vector of event types to target
#' @param TargetTime numeric vector of times at which to evaluate estimands
#' @param MinNuisance numeric; minimum value for propensity score truncation
#'   (passed for documentation; actual truncation happens earlier)
#' @param GComp logical; if TRUE, also compute g-computation plug-in estimates
#'
#' @return Modified Estimates list with added components for each regime:
#'   - SummEIC: data.table summarizing IC (PnEIC, seEIC) by Time and Event
#'   - IC: data.table of individual-level IC values
#'   - GCompEst: (if GComp=TRUE) g-computation estimates
#'
#' @details
#' The efficient influence curve for the cumulative incidence function under
#' intervention is derived from the nonparametric efficiency bound. It involves
#' the clever covariate (h) which depends on the cause-specific hazards and
#' overall survival function.
#'
#' @useDynLib concrete
#' @importFrom Rcpp evalCpp
#' @export
#' @importFrom stats var
#' @keywords internal
getEIC <- function(Estimates, Data, Regime, TargetEvent, TargetTime, MinNuisance, GComp = FALSE) {
    EvalTimes <- attr(Estimates, "Times")
    # Censored <- 0 %in% Data[[attr(Data, "EventType")]]
    T.tilde <- Data[[attr(Data, "EventTime")]]
    Delta <- Data[[attr(Data, "EventType")]]

    # Validate inputs
    if (length(TargetEvent) == 0) {
        stop("TargetEvent is empty. Specify at least one event type to target.",
             call. = FALSE)
    }
    if (length(TargetTime) == 0) {
        stop("TargetTime is empty. Specify at least one time point to target.",
             call. = FALSE)
    }

    for (a in seq_along(Estimates)) {
        NuisanceWeight <- Estimates[[a]][["NuisanceWeight"]]
        GStar <- attr(Estimates[[a]][["PropScore"]], "g.star.obs")
        Hazards <- Estimates[[a]][["Hazards"]]
        TotalSurv <- Estimates[[a]][["EvntFreeSurv"]]

        IC.a <- getIC(GStar = GStar, Hazards = Hazards, TotalSurv = TotalSurv,
                      NuisanceWeight = NuisanceWeight, TargetEvent = TargetEvent,
                      TargetTime = TargetTime, T.tilde = T.tilde,
                      Delta = Delta, EvalTimes = EvalTimes, GComp = GComp)

        if (GComp)
            Estimates[[a]][["GCompEst"]] <- getGComp(EvalTimes, Hazards, TotalSurv, TargetTime)

        Estimates[[a]][["SummEIC"]] <- summarizeIC(IC.a)
        Estimates[[a]][["IC"]] <- IC.a
    }
    return(Estimates)
}

#' Compute Individual-Level Efficient Influence Curves
#'
#' Computes the efficient influence curve (EIC) for each subject, target event,
#' and target time. This is the core computation for TMLE variance estimation.
#'
#' @param GStar numeric vector; intervention probabilities g*(A|W) for each subject
#' @param Hazards named list of hazard matrices (times x subjects) for each event type
#' @param TotalSurv matrix (times x subjects); overall survival probability S(t)
#' @param NuisanceWeight matrix (times x subjects); weights combining g(A|W) and
#'   censoring survival
#' @param TargetEvent numeric vector; event types to compute IC for
#' @param TargetTime numeric vector; time points to compute IC at
#' @param T.tilde numeric vector; observed event/censoring times per subject
#' @param Delta numeric vector; observed event types per subject (0=censored)
#' @param EvalTimes numeric vector; all evaluation time points
#' @param GComp logical; whether g-computation is being performed (for consistency)
#'
#' @return data.table with columns ID, Time, Event, IC containing the efficient
#'   influence curve value for each subject-time-event combination
#'
#' @details
#' The function uses matrix operations with strategic caching to avoid redundant
#' subsetting operations. For each target time tau, the following are pre-computed:
#' - Logical index for times <= tau
#' - Subset of cumulative incidence function F(t)
#' - Subset of survival function S(t)
#' - Subset of nuisance weights
#'
#' @keywords internal
getIC <- function(GStar, Hazards, TotalSurv, NuisanceWeight, TargetEvent, TargetTime,
                  T.tilde, Delta, EvalTimes, GComp) {
    GStar <- as.numeric(unlist(GStar))

    # Delegate all heavy computation to C++ (computeIC).
    # Math: for each (j, tau), IC[i] = point_mass[i] - integral[i] + F_j(tau,i) - mean(F_j(tau))
    # where the loop over L hazard causes is collapsed using h_total = sum_l h_l.
    # See rcpp_getCleverCovariate.cpp::computeIC for details.
    res <- computeIC(
        Hazards        = Hazards,
        TotalSurv      = TotalSurv,
        NuisanceWeight = NuisanceWeight,
        GStar          = GStar,
        T_tilde        = T.tilde,
        Delta          = as.integer(Delta),
        EvalTimes      = EvalTimes,
        TargetEvent    = as.integer(TargetEvent),
        TargetTime     = TargetTime
    )

    # res$IC is N x (J*K); res$Time and res$Event label each column
    N <- nrow(res$IC)
    IC.a <- rbindlist(lapply(seq_len(ncol(res$IC)), function(k) {
        ic_vals <- res$IC[, k]
        if (anyNA(ic_vals))
            stop("IC overflow: either increase MinNuisance or specify a target estimand",
                 " (Target Event, Target Time, & Intervention) with more support in the data.")
        data.table(ID = seq_len(N), Time = res$Time[k], Event = res$Event[k], IC = ic_vals)
    }))
    return(IC.a)
}

getGComp <- function(EvalTimes, Hazards, TotalSurv, TargetTime) {
    F.j.tau <- Event <- Time <- NULL
    Risks <- do.call(rbind, lapply(Hazards, function(haz.j) {
        Risk.a <- sapply(1:ncol(haz.j), function(i) {
            cumsum(TotalSurv[, i] * haz.j[, i])
        })
        Risk.a <- cbind("Event" = as.numeric(attr(haz.j, "j")),
                        "Time" = EvalTimes[EvalTimes %in% TargetTime],
                        "F.j.tau" = rowMeans(subset(Risk.a, EvalTimes %in% TargetTime)))
        return(Risk.a)
    }))
    Risks <- as.data.table(Risks)
    Risks <- rbind(Risks, Risks[, list("Event" = -1, "F.j.tau" = 1 - sum(F.j.tau)), by = "Time"])
    return(Risks[, list("Event" = Event, "Time" = Time, "Risk" = F.j.tau)])
}

summarizeIC <- function(IC.a) {
    IC <- NULL
    IC.a <- rbind(IC.a,
                  IC.a[, list("Event" = -1, "IC" = -sum(IC)), by = c("ID", "Time")])
    summIC <- IC.a[, list("PnEIC" = mean(IC), 
                          "seEIC" = sqrt(mean(IC^2)),
                          "seEIC/(sqrt(n)log(n))" = 
                            sqrt(mean(IC^2)) / (sqrt(.N) * log(.N))),
                   by = c("Time", "Event")]
    return(summIC)
}
