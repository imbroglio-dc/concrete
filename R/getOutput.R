#' getOutput
#'
#' @param Estimate list : a "ConcreteEst" object
#' @param Estimand character (default: c("RD", "RR", "Risk")) or a function or a list of functions 
# #' @param plot boolean (default = TRUE) : whether or not to plot the 
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

getOutput <- function(Estimate, Estimand = c("RD", "RR", "Risk")) {
    if (!all(sapply(Estimand, function(e) any(is.function(e), grepl("(rd)|(rr)|(risk)", tolower(e)))))) {
        stop("Estimand must be in c('RD', 'RR', 'Risk'), or be a list of user-specified function(s) of",
             "`Estimate`, `Estimand`, `TargetEvent`, `TargetTime`, and `GComp`.")
    }
    if (any(sapply(Estimand, function(e) is.function(e)))) {
        cat("User-specified Estimand functions are not checked or tested. ", 
            "Please verify the correctness of your custom Estimand function(s).")
    }
    TargetTime <- attr(Estimate, "TargetTime")
    TargetEvent <- attr(Estimate, "TargetEvent")
    GComp <- attr(Estimate, "GComp")
    
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
    
    if (any("risk" %in% tolower(Estimand))) {
        risks <- getRisk(Estimate = Estimate, TargetTime = TargetTime, TargetEvent = TargetEvent, GComp = GComp)
        risks <- do.call(rbind, lapply(seq_along(risks), function(a) {
            risk.a <- risks[[a]]
            a.name <- names(risks)[a]
            return(cbind("Intervention" = a.name, risk.a))
        }))
        risks <- structure(risks,
                           Estimand = "Absolute Risks", 
                           class = union("ConcreteOut", class(risks)))
        output[["Risk"]] <- risks
    }
    return(output)
}

getRD <- function(Estimate, TargetTime, TargetEvent, GComp) {
    Estimator <- Event <- Time <- seEIC <- Risk.x <- Risk.y <- se.x <- se.y <- NULL
    risk <- getRisk(Estimate = Estimate, TargetTime = TargetTime, TargetEvent = TargetEvent, GComp = GComp)
    rd.out <- as.data.table(merge(risk[[1]], risk[[2]], by = c("Estimator", "Event", "Time")))[order(Time)]
    rd.out <- rd.out[, list("Estimator" = Estimator, "Event" = Event, "Time" = Time, "RD" = Risk.x - Risk.y,
                            se = sqrt(se.x^2 + se.y^2))]
    rd.out <- structure(rd.out, 
                        Interventions = paste0(names(Estimate), collapse = " - "), 
                        Estimand = "Risk Difference", 
                        class = union("ConcreteOut", class(rd.out)))
    return(rd.out[order(Estimator, decreasing = TRUE)])
}

getRR <- function(Estimate, TargetTime, TargetEvent, GComp) {
    Estimator <- Event <- Time <- seEIC <- Risk.x <- Risk.y <- se.x <- se.y <- NULL
    risk <- getRisk(Estimate = Estimate, TargetTime = TargetTime, TargetEvent = TargetEvent, GComp = GComp)
    rr.out <- as.data.table(merge(risk[[1]], risk[[2]], by = c("Estimator", "Event", "Time")))[order(Time)]
    rr.out <- rr.out[, list("Estimator" = Estimator, "Event" = Event, "Time" = Time, "RR" = Risk.x / Risk.y,
                            se = sqrt((se.x / Risk.y)^2 + (se.y * Risk.x / Risk.y^2)^2))]
    rr.out <- structure(rr.out, 
                        Interventions = paste0(names(Estimate), collapse = " / "), 
                        Estimand = "Relative Risk", 
                        class = union("ConcreteOut", class(rr.out)))
    return(rr.out[order(Estimator, decreasing = TRUE)])
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
