#' getOutput
#'
#' @param Estimate list : a "ConcreteEst" object
#' @param Estimand character: a subset of c("RD", "RR", "Risk") specifying the target estimand
#' @param A1 default = 1, the ConcreteEst list element that is the "Treated" arm for risk comparisons
#' @param A0 default = 2, the ConcreteEst list element that is the "Control" arm for risk comparisons
#'
#' @return data.table of point estimates and standard deviations
#' @export getOutput
#'
#' @examples
#' # getOutput returns risk difference, relative risk, and treatment-specific risks
#' # concrete.out <- getOutput(Estimate = concrete.est, Estimand = c("rd", "rr", "risk"),
#' #                           TargetTime = target.time, TargetEvent = target.event, GComp = TRUE)
#' # concrete.out$RD
#' # concrete.out$RR
#' # concrete.out$Risk
#' 
#' @importFrom MASS mvrnorm
#' @importFrom stats cor cov

getOutput <- function(Estimate, Estimand = c("RD", "RR", "Risk"), A1 = 1, A0 = 2) {
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
    
    Output <- data.table()
    Risks <- getRisk(Estimate = Estimate, TargetTime = TargetTime, TargetEvent = TargetEvent, GComp = GComp)
    if (any(sapply(Estimand, is.function))) {
        Output <- lapply(Estimand[which(sapply(Estimand, is.function))], function(estimand) {
            do.call(estimand, list("Estimate" = Estimate, "TargetTime" = TargetTime,
                                   "TargetEvent" = TargetEvent, "GComp" = GComp))
        })
        names(Output) <- names(Estimand)[which(sapply(Estimand, is.function))]
    }
    
    if (any(grepl("risk", tolower(Estimand))))
        Output <- rbind(Output, Risks)
    
    if (any(grepl("rd", tolower(Estimand)))) 
        Output <- rbind(Output, 
                        getRD(Risks = Risks, A1 = names(Estimate)[A1], A0 = names(Estimate)[A0], 
                              TargetTime = TargetTime, TargetEvent = TargetEvent, GComp = GComp))
    
    if (any(grepl("rr", tolower(Estimand)))) 
        Output <- rbind(Output, 
                        getRR(Risks = Risks, A1 = names(Estimate)[A1], A0 = names(Estimate)[A0], 
                              TargetTime = TargetTime, TargetEvent = TargetEvent, GComp = GComp))
    
    return(Output)
}

getRisk <- function(Estimate, TargetTime, TargetEvent, GComp) {
    Event <- Time <- seEIC <- Estimand <- NULL
    risk <- lapply(seq_along(Estimate), function(a) {
        est.a <- Estimate[[a]]
        risk.a <- do.call(rbind, lapply(as.character(TargetEvent), function(j) {
            risks <- apply(est.a[["Hazards"]][[j]] * est.a[["EvntFreeSurv"]], 2, cumsum)
            Psi <- cbind("Time" = attr(Estimate, "Times")[which(attr(Estimate, "Times") %in% TargetTime)],
                         "Risk" = rowMeans(subset(risks, subset = attr(Estimate, "Times") %in% TargetTime)))
            se <- subset(est.a$SummEIC[Event == j, ], select = c("Time", "seEIC"))
            Psi <- merge(Psi, se[, list("Time" = Time, "se" = seEIC / sqrt(ncol(risks)))], by = "Time")
            return(cbind("Estimator" = "tmle", "Event" = as.numeric(j), Psi))
        }))
        if (GComp)
            risk.a <- rbind(risk.a,
                            cbind("Estimator" = "gcomp",
                                  est.a$GCompEst[Event %in% TargetEvent, ],
                                  "se" = NA))
        risk.a <- cbind(Intervention = names(Estimate)[a], 
                        risk.a)
        return(risk.a)
    })
    Risk <- setDT(do.call(rbind, risk))
    Risk[, Estimand := "Abs. Risk"]
    setnames(Risk, "Risk", "Pt Est")
    setcolorder(Risk, c("Intervention", "Estimand", "Estimator", "Event", "Time", "Pt Est", "se"))
    return(Risk)
}

getRD <- function(Risks, A1, A0, TargetTime, TargetEvent, GComp) {
    `Pt Est` <- se <- Intervention <- Estimand <- NULL
    RD <- Risks[, list("Pt Est" = `Pt Est`[Intervention == A1] - `Pt Est`[Intervention == A0], 
                       se = sqrt(se[Intervention == A1]^2 + se[Intervention == A0]^2)), 
                by = c("Estimator", "Event", "Time")]
    RD[, Intervention := paste0("[", A1, "] - [", A0, "]")]
    RD[, Estimand := "Risk Diff"]
    setcolorder(RD, c("Intervention", "Estimand", "Estimator", "Event", "Time", "Pt Est", "se"))
    return(RD)
}

getRR <- function(Risks, A1, A0, TargetTime, TargetEvent, GComp) {
    `Pt Est` <- se <- Intervention <- Estimand <- NULL
    RR <- Risks[, list("Pt Est" = `Pt Est`[Intervention == A1] / `Pt Est`[Intervention == A0], 
                       se = sqrt((se[Intervention == A1] / `Pt Est`[Intervention == A0])^2 + 
                                     se[Intervention == A0]^2 * 
                                     (`Pt Est`[Intervention == A1] / `Pt Est`[Intervention == A0]^2)^2)), 
                by = c("Estimator", "Event", "Time")]
    RR[, Intervention := paste0("[", A1, "] / [", A0, "]")]
    RR[, Estimand := "Rel Risk"]
    setcolorder(RR, c("Intervention", "Estimand", "Estimator", "Event", "Time", "Pt Est", "se"))
    return(RR)
}

getSimultaneous <- function(Estimate, A1 = 1, A0 = 2, Parameter = c("Risk", "RR", "RD"), 
                            SignifLevel = 0.05) {
    Intervention <- IC <- Risk <- NULL
    A1 <- names(Estimate)[A1]
    A0 <- names(Estimate)[A0]
    
    ICs <- do.call(rbind, 
                   lapply(seq_along(Estimate), function(a) {
                       est.a <- Estimate[[a]]
                       IC.a <- getIC(GStar =  attr(est.a[["PropScore"]], "g.star.obs"),
                                     Hazards = est.a[["Hazards"]], 
                                     TotalSurv = est.a[["EvntFreeSurv"]],
                                     NuisanceWeight = est.a[["NuisanceWeight"]],
                                     TargetEvent = attr(Estimate, "TargetEvent"), 
                                     TargetTime = attr(Estimate, "TargetTime"),
                                     T.tilde = attr(Estimate, "T.tilde"), 
                                     Delta = attr(Estimate, "Delta"),
                                     EvalTimes = attr(Estimate, "Times"), 
                                     GComp = FALSE)
                       return(cbind("Intervention" = names(Estimate)[a], IC.a))
                   }))
    Risks <- do.call(rbind, 
                     getRisk(Estimate = Estimate, TargetTime = attr(Estimate, "TargetTime"), 
                             TargetEvent = attr(Estimate, "TargetEvent"), GComp = FALSE))
    ICs <- merge(ICs, subset(Risks, select = c("Intervention", "Time", "Event", "Risk")), 
                 by = c("Intervention", "Time", "Event"))
    
    ICs <- rbind(subset(ICs, select = c("ID", "Time", "Event", "Intervention", "IC")), 
                 ICs[, list("Intervention" = "RD",
                            "IC" = IC[Intervention == A1] - IC[Intervention == A0]), 
                     by = c("ID", "Time", "Event")], 
                 ICs[, list("Intervention" = "RR",
                            "IC" = IC[Intervention == A1] / Risk[Intervention == A0] - 
                                IC[Intervention == A0] * Risk[Intervention == A1] / Risk[Intervention == A0]^2), 
                     by = c("ID", "Time", "Event")])
    WideICs <- dcast(ICs, ID ~ Intervention + Time + Event, value.var = "IC")[, !c("ID")]
    
    IndicesRR <- which(grepl("RR", colnames(WideICs)))
    IndicesRD <- which(grepl("RD", colnames(WideICs)))
    IndicesRisk <- setdiff(1:ncol(WideICs), c(IndicesRR, IndicesRD))
    IndicesKeep <- c()
    if (any(grepl("risk", tolower(Parameter))))
        IndicesKeep <- c(IndicesKeep, IndicesRisk)
    if (any(grepl("rd", tolower(Parameter))))
        IndicesKeep <- c(IndicesKeep, IndicesRD)
    if (any(grepl("rr", tolower(Parameter))))
        IndicesKeep <- c(IndicesKeep, IndicesRR)
    
    WideICs <- subset(WideICs, select = IndicesKeep)
    CovEIC <- cov(WideICs)
    CovEIC[is.na(CovEIC)] <- 1e-9
    CorrEIC <- cor(WideICs)
    n <- length(attr(Estimate, "T.tilde"))
    
    q <- apply(abs(MASS::mvrnorm(n = 1e3, mu = rep(0, nrow(CorrEIC)), Sigma = CorrEIC)), 1, max)
    q <- as.numeric(stats::quantile(q, 1 - SignifLevel))
    se <- data.table(names = rownames(CovEIC), se = sqrt(diag(CovEIC)) / sqrt(n) * q)
    se[, c("Intervention", "Time", "Event") := tstrsplit(names, "_")]
    se <- subset(se, select = c("Intervention", "Time", "Event", "se"))
    return(se)
}
