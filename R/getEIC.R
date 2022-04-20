#' get EICs
#'
#' @param Estimates list
#' @param Data data.table
#' @param Regime list
#' @param Censored boolean
#' @param TargetEvent numeric vector
#' @param TargetTime numeric vector
#' @param Events numeric vector
#' @param MinNuisance numeric
#' @param GComp boolean
#'
#'

getEIC <- function(Estimates, Data, Regime, Censored, TargetEvent, TargetTime,
                   Events, MinNuisance, GComp = FALSE) {
    EvalTimes <- attr(Estimates, "times")
    T.tilde <- Data[[attr(Data, "EventTime")]]
    Delta <- Data[[attr(Data, "EventType")]]

    for (a in seq_along(Estimates)) {
        NuisanceWeight <- Estimates[[a]][["NuisanceWeight"]]
        GStar <- attr(Estimates[[a]][["PropScore"]], "g.star.obs")
        Hazards <- Estimates[[a]][["Hazards"]]
        TotalSurv <- Estimates[[a]][["EvntFreeSurv"]]

        IC.a <- getIC(GStar = GStar, Hazards = Hazards, TotalSurv = TotalSurv,
                      NuisanceWeight = NuisanceWeight, TargetEvent = TargetEvent,
                      TargetTime = TargetTime, Events = Events, T.tilde = T.tilde,
                      Delta = Delta, EvalTimes = EvalTimes, GComp = GComp)

        if (GComp)
            Estimates[[a]][["GCompEst"]] <- getGComp(EvalTimes, Hazards, TotalSurv, TargetTime)

        Estimates[[a]][["SummEIC"]] <- summarizeIC(IC.a)
    }
    return(Estimates)
}

getIC <- function(GStar, Hazards, TotalSurv, NuisanceWeight, TargetEvent, TargetTime,
                  Events, T.tilde, Delta, EvalTimes, GComp) {
    Target <- expand.grid("Time" = TargetTime, "Event" = TargetEvent)
    IC <- F.j.tau <- NULL
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

            IC.jk <- sum(sapply(Events, function(l) {
                h.jk <- GStar[i] * Nuisance.ik * ((l == j) - (F.j.tau - F.j.t) / Surv.ik)
                IC.ljk <- sum(h.jk * ((s.ik == t.tilde) * (Delta[i] == l) - haz.J.ik[[as.character(l)]]))
                ## the second EIC component ( ... + F_j(tau | a, L) - Psi ) is done outside
                return(IC.ljk)
            }))
            return(c("IC" = IC.jk, "F.j.tau" = F.j.tau))
        }))
        IC.jk <- cbind("ID" = i, Target,  IC.jk)
        return(IC.jk)
    }))
    ## the second EIC component ( ... + F_j(tau | a, L) - Psi )
    IC.a <- as.data.table(IC.a)
    IC.a[, IC := IC + F.j.tau - mean(F.j.tau), by = c("Time", "Event")]
    return(IC.a[,.SD, .SDcols = !"F.j.tau"])
}

getGComp <- function(EvalTimes, Hazards, TotalSurv, TargetTime) {
    F.j.tau <- NULL
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
    return(as.data.table(Risks))
}

summarizeIC <- function(IC.a) {
    IC <- NULL
    IC.a <- rbind(IC.a,
                  IC.a[, list("Event" = -1, "IC" = -sum(IC)), by = c("ID", "Time")])
    return(IC.a[, list("PnEIC" = mean(IC), "seEIC" = sqrt(stats::var(IC)),
                       "seEIC/(sqrt(n)log(n))" = sqrt(stats::var(IC)/.N)/log(.N)),
                by = c("Time", "Event")])
}
