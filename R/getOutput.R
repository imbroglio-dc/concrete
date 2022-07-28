#' getOutput
#'
#' @param Estimate : list concrete Estimate object
#' @param Estimand : character or list of function(s)
#' @param TargetTime : numeric
#' @param TargetEvent : numeric
#' @param GComp : boolean
#'
#' @return tbd
#' @export getOutput
#'
#' @examples
#' # getOutput returns risk difference, relative risk, and treatment-specific risks
#' # concrete.out <- getOutput(Estimate = concrete.est, Estimand = c("rd", "rr", "risk"), 
#' #                           TargetTime = target.time, TargetEvent = target.event, GComp = TRUE)
#' # concrete.out$RD
#' # concrete.out$RR
#' # concrete.out$Risk

getOutput <- function(Estimate, Estimand = c("RD", "RR", "Risk"), TargetTime, TargetEvent, GComp) {
    if (!all(sapply(Estimand, function(e) any(is.function(e), grepl("(rd)|(rr)|(risk)", tolower(e)))))) {
        stop("Estimand must be in c('RD', 'RR', 'Risk'), or be a list of user-specified function(s) of",
             "`Estimate`, `Estimand`, `TargetEvent`, `TargetTime`, and `GComp`.")
    }
    output <- list()
    if (any(sapply(Estimand, is.function))) {
        output <- lapply(Estimand[which(sapply(Estimand, is.function))], function(estimand) {
            do.call(estimand, list("Estimate" = Estimate, "TargetTime" = TargetTime,
                                   "TargetEvent" = TargetEvent, "GComp" = GComp))
        })
        names(output) <- names(Estimand)[which(sapply(Estimand, is.function))]
    }
    if (any(tolower(Estimand) == "rd")) {
        output[["RD"]] <- getRD(Estimate = Estimate, TargetTime = TargetTime, TargetEvent = TargetEvent, GComp = GComp)
    }

    if (any(tolower(Estimand) == "rr")) {
        output[["RR"]] <- getRR(Estimate = Estimate, TargetTime = TargetTime, TargetEvent = TargetEvent, GComp = GComp)
    }

    if (any(tolower(Estimand) == "risk")) {
        output[["Risk"]] <- getRisk(Estimate = Estimate, TargetTime = TargetTime, TargetEvent = TargetEvent, GComp = GComp)
    }

    return(output)
}

getRD <- function(Estimate, TargetTime, TargetEvent, GComp) {
    Estimator <- Event <- Time <- seEIC <- Risk.x <- Risk.y <- se.x <- se.y <- NULL
    risk <- getRisk(Estimate = Estimate, TargetTime = TargetTime, TargetEvent = TargetEvent, GComp = GComp)
    rd.out <- as.data.table(merge(risk[[1]], risk[[2]], by = c("Estimator", "Event", "Time")))[order(Time)]
    rd.out <- rd.out[, list("Estimator" = Estimator, "Event" = Event, "Time" = Time, "RD" = Risk.x - Risk.y,
                            se = sqrt(se.x^2 + se.y^2))]
    attr(rd.out, "regimes") <- paste0(names(Estimate), collapse = " - ")
    return(rd.out)
}

getRR <- function(Estimate, TargetTime, TargetEvent, GComp) {
    Estimator <- Event <- Time <- seEIC <- Risk.x <- Risk.y <- se.x <- se.y <- NULL
    risk <- getRisk(Estimate = Estimate, TargetTime = TargetTime, TargetEvent = TargetEvent, GComp = GComp)
    rr.out <- as.data.table(merge(risk[[1]], risk[[2]], by = c("Estimator", "Event", "Time")))[order(Time)]
    rr.out <- rr.out[, list("Estimator" = Estimator, "Event" = Event, "Time" = Time, "RR" = Risk.x / Risk.y,
                            se = sqrt((se.x / Risk.y)^2 + (se.y * Risk.x / Risk.y^2)^2))]
    attr(rr.out, "regimes") <- paste0(names(Estimate), collapse = " / ")
    return(rr.out)
}

getRisk <- function(Estimate, TargetTime, TargetEvent, GComp) {
    Event <- Time <- seEIC <- NULL
    risk <- lapply(Estimate, function(est.a) {
        risk.a <- do.call(rbind, lapply(as.character(TargetEvent), function(j) {
            risks <- apply(est.a[["Hazards"]][[j]] * est.a[["EvntFreeSurv"]], 2, cumsum)
            Psi <- cbind("Time" = attr(Estimate, "times")[which(attr(Estimate, "times") %in% TargetTime)],
                         "Risk" = rowMeans(subset(risks, subset = attr(Estimate, "times") %in% TargetTime)))
            se <- subset(est.a$SummEIC[Event == j, ], select = c("Time", "seEIC"))
            Psi <- merge(Psi, se[, list("Time" = Time, "se" = seEIC / sqrt(ncol(risks)))], by = "Time")
            return(cbind("Estimator" = "tmle", "Event" = j, Psi))
        }))
        if (GComp)
            risk.a <- rbind(risk.a,
                            cbind("Estimator" = "gcomp",
                                  est.a$GCompEst[Event %in% TargetEvent, ],
                                  "se" = NA))
        return(risk.a)
    })
    return(risk)
}
