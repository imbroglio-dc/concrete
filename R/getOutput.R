#' getOutput
#'
#' @param Estimate list : a "ConcreteEst" object
#' @param Estimand character: a subset of c("RD", "RR", "Risk") specifying the target estimand
#' @param Intervention numeric (default = seq_along(ConcreteEst)): the ConcreteEst list element 
#'          corresponding to the target intervention. For comparison estimands such as RD and RR, 
#'          Intervention should be a numeric vector with length 2, the first term designating 
#'          "treatment" ConcreteEst list element and the second designating the "control".
#' @param GComp logical: to return g-formula point estimates based on initial nuisance parameter estimation
#' @param Simultaneous logical: to compute simultaneous confidence intervals
#' @param Signif numeric (default = 0.05): alpha for 2-tailed hypothesis testing
#'
#' @return data.table of point estimates and standard deviations
#' @export getOutput
#'
#' @examples
#' # getOutput returns risk difference, relative risk, and treatment-specific risks
#' # concrete.out <- getOutput(Estimate = concrete.est, Estimand = "RR", Intervention = 1:2)
#' # concrete.out$RD
#' # concrete.out$RR
#' # concrete.out$Risk
#' 
#' @importFrom MASS mvrnorm
#' @importFrom stats cor cov
#' 

getOutput <- function(Estimate, Estimand = c("RD", "RR", "Risk"), Intervention = seq_along(Estimate), 
                      GComp = NULL, Simultaneous = TRUE, Signif = 0.05) {
    `CI Low` <- `CI Hi` <- `Pt Est` <- `se` <- NULL
    if (!inherits(Estimate, "ConcreteEst")) 
        stop("Estimate must be a 'ConcreteEst' class object")
    if (!all(sapply(Estimand, function(e) ifelse(is.function(e), TRUE, 
                                                 grepl("(rd)|(rr)|(risk)|(risks)", tolower(e)))))) {
        stop("Estimand must be in c('RD', 'RR', 'Risk')")
    } else {
        RiskWanted <- any(sapply(Estimand, function(e) ifelse(is.function(e), FALSE, grepl("risk", tolower(e)))))
        RRWanted <- any(sapply(Estimand, function(e) ifelse(is.function(e), FALSE, "rr" == tolower(e))))
        RDWanted <- any(sapply(Estimand, function(e) ifelse(is.function(e), FALSE, "rd" == tolower(e))))
    }
    if (any(sapply(Estimand, function(e) is.function(e)))) {
        cat("User-specified Estimand functions are not checked or tested. ",
            "Please verify the correctness of your custom Estimand function(s).\n", sep = "")
    }
    if (!is.numeric(Intervention) & !is.integer(Intervention))
        stop("Intervention must be a numeric vector, specifying which 'Estimate' list elements", 
             "are of interest.")
    if (RRWanted | RDWanted) {
        if (length(Intervention) < 2)
            stop("Risk Ratios and Risk Differences can only be computed if Interventions is a vector of ", 
                 "length 2, i.e. c('treated' Estimate list index, 'control' Estimate list index). Either", 
                 "specify at least two indices in Intervention or remove 'RR' and 'RD' from Estimand.")
        if (length(Intervention) > 2)
            cat("Risk ratios and risk differences will be computed using only the first two ",
                "elements of the provided 'Intervention' argument.\n", sep = "")
    }
    if (!is.logical(Simultaneous))
        stop("Simultaneous must be a logical, TRUE or FALSE")
    if (!is.numeric(Signif) | Signif >= 1 | Signif <= 0)
        stop("Signif sets the alpha (significance level) for hypothesis testing and confidence intervals", 
             ", so should be small number that must be greater than 0")
    
    TargetTime <- attr(Estimate, "TargetTime")
    TargetEvent <- attr(Estimate, "TargetEvent")
    if (is.null(GComp))
        GComp <- attr(Estimate, "GComp")
    if (!is.logical(GComp))
        stop("GComp must be a logical, TRUE or FALSE")
    
    Output <- data.table()
    Risks <- getRisk(Estimate = Estimate, TargetTime = TargetTime, TargetEvent = TargetEvent, GComp = GComp)
    if (any(sapply(Estimand, is.function))) {
        Output <- lapply(Estimand[which(sapply(Estimand, is.function))], function(estimand) {
            do.call(estimand, list("Estimate" = Estimate, "TargetTime" = TargetTime,
                                   "TargetEvent" = TargetEvent, "GComp" = GComp))
        })
        names(Output) <- names(Estimand)[which(sapply(Estimand, is.function))]
    }
    
    if (RiskWanted) 
        Output <- rbind(Output, Risks)
    
    if (RDWanted) 
        Output <- rbind(Output, 
                        getRD(Risks = Risks, A1 = names(Estimate)[Intervention[1]], 
                              A0 = names(Estimate)[Intervention[2]], TargetTime = TargetTime, 
                              TargetEvent = TargetEvent, GComp = GComp))
    
    if (RRWanted) 
        Output <- rbind(Output, 
                        getRR(Risks = Risks, A1 = names(Estimate)[Intervention[1]], 
                              A0 = names(Estimate)[Intervention[2]], TargetTime = TargetTime, 
                              TargetEvent = TargetEvent, GComp = GComp))
    
    Output[, `CI Low` := `Pt Est` - qnorm(1 - Signif/2)*se]
    Output[, `CI Hi` := `Pt Est` + qnorm(1 - Signif/2)*se]
    
    if (Simultaneous)
        Output <- getSimultaneous(Estimate = Estimate, Risks = Output, RiskWanted = RiskWanted, 
                                  RDWanted = RDWanted, RRWanted = RRWanted, 
                                  Intervention = Intervention, Signif = Signif)
    class(Output) <- union("ConcreteOut", class(Output))
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

getSimultaneous <- function(Estimate, Risks, RiskWanted, RDWanted, RRWanted, Intervention, Signif) {
    Estimator <- `Pt Est` <- Event <- Time <- SimQ <- IC <- Risk <- NULL
    
    if (RRWanted | RDWanted) {
        A1 <- names(Estimate)[Intervention[1]]
        A0 <- names(Estimate)[Intervention[2]]
    }
    
    ICs <- do.call(rbind, 
                   lapply(Intervention, function(a) {
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
    
    ICs <- merge(ICs, getRisk(Estimate = Estimate, TargetTime = attr(Estimate, "TargetTime"), 
                              TargetEvent = attr(Estimate, "TargetEvent"), GComp = FALSE), 
                 by = c("Intervention", "Time", "Event"))
    if (RDWanted)
        RDICs <- ICs[, list("Intervention" = "Risk Diff",
                            "IC" = IC[Intervention == A1] - IC[Intervention == A0]), 
                     by = c("ID", "Time", "Event")]
    if (RRWanted)
        RRICs <- ICs[, list("Intervention" = "Rel Risk",
                            "IC" = IC[Intervention == A1] / `Pt Est`[Intervention == A0] - 
                                IC[Intervention == A0] * `Pt Est`[Intervention == A1] / 
                                `Pt Est`[Intervention == A0]^2), 
                     by = c("ID", "Time", "Event")]
    if (RiskWanted) {
        ICs <- subset(ICs, select = c("ID", "Time", "Event", "Intervention", "IC"))
    } else {
        ICs <- data.table()
    }
    if (exists("RDICs")) ICs <- rbind(ICs, RDICs)
    if (exists("RRICs")) ICs <- rbind(ICs, RRICs)
    
    ICs <- dcast(ICs, ID ~ Intervention + Time + Event, value.var = "IC")[, !c("ID")]
    CovEIC <- cov(ICs)
    CovEIC[is.na(CovEIC)] <- 1e-9
    CorrEIC <- cor(ICs)
    n <- length(attr(Estimate, "T.tilde"))
    
    q <- apply(abs(MASS::mvrnorm(n = 1e3, mu = rep(0, nrow(CorrEIC)), Sigma = CorrEIC)), 1, max)
    q <- as.numeric(stats::quantile(q, 1 - Signif))
    se <- data.table(names = rownames(CovEIC), SimQ = q)
    se[, c("Intervention", "Time", "Event") := tstrsplit(names, "_")]
    se[, Estimator := "tmle"][, Event := as.numeric(Event)][, Time := as.numeric(Time)]
    se[Intervention == "Rel Risk", Intervention := paste0("[", A1, "] / [", A0, "]")]
    se[Intervention == "Risk Diff", Intervention := paste0("[", A1, "] - [", A0, "]")]
    simCI <- merge(Risks, se, c("Intervention", "Estimator", "Event", "Time"), all.x = TRUE)
    simCI[, "SimCI Low" := `Pt Est` - se*SimQ][, "SimCI Hi" := `Pt Est` + se*SimQ]
    simCI <- subset(simCI, select = c("Intervention", "Estimand", "Estimator", "Event", "Time", 
                                      "Pt Est", "se","CI Low", "CI Hi", "SimCI Low", "SimCI Hi"))
    return(simCI)
}


#' @describeIn doConcrete plot.ConcreteOut plot method for "ConcreteOut" class
#' @param x a ConcreteOut object
#' @param Estimand character: "rr" to plot Relative Risks, "rd" to plot Risk Differences, and "risk" to plot absolute risks
#' @param NullLine logical: to add a red line at 1 for RR plots and at 0 for RD plots
#' @param GComp logical: to plot the g-comp point estimates
#' @param ask logical: to prompt for user input before each plot
#' @param ... additional arguments to be passed into plot methods
#' @exportS3Method plot ConcreteOut
plot.ConcreteOut <- function(x, Estimand = c("rr", "rd", "risk"), NullLine = TRUE, GComp = FALSE, 
                             ask = TRUE, ...) {
    Event <- Time <- `Pt Est` <- Estimator <- se <- NULL
    if(!requireNamespace("ggplot2", quietly = TRUE))
        stop("Plotting requires the 'ggplot2' package")
    if (ask) {
        oask <- devAskNewPage(TRUE)
        on.exit(devAskNewPage(oask))
    }
    fig <- list()
    if (any(grepl("risk", tolower(Estimand)))) {
        message("risk plots not yet implemented")
    }
    if (any(tolower(Estimand) == "rr")) {
        dev.hold()
        z <- x[Estimand == "Rel Risk", ][, Event := paste0("Event ", Event)]
        if (GComp) {
            fig[["rr"]] <- ggplot2::ggplot(z, ggplot2::aes(x = Time, y = `Pt Est`, shape = Estimator, 
                                                           ymin = `Pt Est` - 1.96*se, ymax = `Pt Est` + 1.96*se)) + 
                ggplot2::geom_point() + ggplot2::geom_errorbar() + ggplot2::scale_shape_manual(values = c(4, 16))
        } else {
            fig[["rr"]] <- ggplot2::ggplot(z[Estimator == "tmle", ], ggplot2::aes(x = Time, y = `Pt Est`, 
                                                                                  ymin = `Pt Est` - 1.96*se, 
                                                                                  ymax = `Pt Est` + 1.96*se)) + 
                ggplot2::geom_point() + ggplot2::geom_errorbar()
        }
        fig[["rr"]] <- fig[["rr"]] + 
            ggplot2::facet_wrap(~Event, scales = "free", ncol = 1) + 
            ggplot2::theme_minimal() +
            ggplot2::labs(y = "Relative Risk", x = "Time", 
                          title = "Relative Risk Point Estimates with 95% confidence intervals") 
        if (NullLine)
            fig[["rr"]] <- fig[["rr"]] + ggplot2::geom_hline(ggplot2::aes(yintercept = 1), colour = "red", alpha = 0.4)
        plot(fig[["rr"]])
        dev.flush()
    }
    if (any(tolower(Estimand) == "rd")) {
        dev.hold()
        z <- x[Estimand == "Risk Diff", ][, Event := paste0("Event ", Event)]
        if (GComp) {
            fig[["rd"]] <- ggplot2::ggplot(z, ggplot2::aes(x = Time, y = `Pt Est`, shape = Estimator, 
                                                           ymin = `Pt Est` - 1.96*se, ymax = `Pt Est` + 1.96*se)) + 
                ggplot2::geom_point() + ggplot2::geom_errorbar() + ggplot2::scale_shape_manual(values = c(4, 16))
        } else {
            fig[["rd"]] <- ggplot2::ggplot(z[Estimator == "tmle", ], ggplot2::aes(x = Time, y = `Pt Est`, 
                                                                                  ymin = `Pt Est` - 1.96*se, 
                                                                                  ymax = `Pt Est` + 1.96*se)) + 
                ggplot2::geom_point() + ggplot2::geom_errorbar()
        }
        fig[["rd"]] <- fig[["rd"]] + 
            ggplot2::facet_wrap(~Event, scales = "free", ncol = 1) + 
            ggplot2::theme_minimal() +
            ggplot2::labs(y = "Risk Difference", x = "Time", 
                          title = "Risk Difference Point Estimates with 95% Confidence Intervals") 
        if (NullLine)
            fig[["rd"]] <- fig[["rd"]] + ggplot2::geom_hline(ggplot2::aes(yintercept = 0), colour = "red", alpha = 0.4)
        plot(fig[["rd"]])
        dev.flush()
    }
    invisible(fig)
}
