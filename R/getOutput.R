#' getOutput
#'
#' @param ConcreteEst "ConcreteEst" object
#' @param Estimand character: "RR" for Relative Risks, "RD" for Risk Differences, and "Risk" for absolute risks
#' @param Intervention numeric (default = seq_along(ConcreteEst)): the ConcreteEst list element 
#'          corresponding to the target intervention. For comparison estimands such as RD and RR, 
#'          Intervention should be a numeric vector with length 2, the first term designating 
#'          "treatment" ConcreteEst list element and the second designating the "control".
#' @param GComp logical: return g-formula point estimates based on initial nuisance parameter estimation
#' @param Simultaneous logical: return simultaneous confidence intervals
#' @param Signif numeric (default = 0.05): alpha for 2-tailed hypothesis testing
#'
#' @return data.table of point estimates and standard deviations
#' @export getOutput
#'
#' @examples
#' library(data.table)
#' library(concrete)
#' 
#' data <- as.data.table(survival::pbc)
#' data <- data[1:200, .SD, .SDcols = c("id", "time", "status", "trt", "age", "sex")]
#' data[, trt := sample(0:1, nrow(data), TRUE)]
#' 
#' # formatArguments() returns correctly formatted arguments for doConcrete()
#' concrete.args <- formatArguments(DataTable = data,
#'                                  EventTime = "time",
#'                                  EventType = "status",
#'                                  Treatment = "trt",
#'                                  ID = "id",
#'                                  TargetTime = 2500,
#'                                  TargetEvent = c(1, 2),
#'                                  Intervention = makeITT(),
#'                                  CVArg = list(V = 2))
#'                                  
#' # doConcrete() returns tmle (and g-formula plug-in) estimates of targeted risks
#' \donttest{
#' concrete.est <- doConcrete(concrete.args)
#' 
#' # getOutput returns risk difference, relative risk, and treatment-specific risks
#' #  GComp=TRUE returns g-formula plug-in estimates
#' #  Simultaneous=TRUE computes simultaneous CI for all output TMLE estimates 
#' concrete.out <- getOutput(concrete.est, Estimand = c("RR", "RD", "Risk"),
#'                           GComp = TRUE, Simultaneous = TRUE)
#' print(concrete.out)
#' plot(concrete.out, ask = FALSE)
#' }
#' @importFrom MASS mvrnorm
#' @importFrom stats cor cov

getOutput <- function(ConcreteEst, Estimand = c("Risk"), Intervention = seq_along(ConcreteEst), 
                      GComp = NULL, Simultaneous = TRUE, Signif = 0.05) {
  `CI Low` <- `CI Hi` <- `Pt Est` <- `se` <- NULL
  if (!inherits(ConcreteEst, "ConcreteEst")) 
    stop("ConcreteEst must be a 'ConcreteEst' class object")
  TargetTime <- attr(ConcreteEst, "TargetTime")
  TargetEvent <- attr(ConcreteEst, "TargetEvent")
  if (is.null(GComp))
    GComp <- attr(ConcreteEst, "GComp")
  if (!is.logical(GComp))
    stop("GComp must be a logical, TRUE or FALSE")
  
  EstimandType <- sapply(Estimand, function(e) {
    if (is.function(e)) return("Function")
    if (grepl("rd", tolower(e))) return("RD")
    if (grepl("rr", tolower(e))) return("RR")
    if (grepl("risk", tolower(e))) return("Risk")
    else stop("Estimand inputs must be in c('RD', 'RR', 'Risk')")
  })
  if (any(EstimandType == "Function")) {
    warning("User-specified Estimand functions are not checked or tested. ",
            "Please verify the correctness of your custom Estimand function(s).\n", sep = "")
  }
  
  if (!is.numeric(Intervention) & !is.integer(Intervention))
    stop("Intervention must be a numeric vector, specifying which 'ConcreteEst' list elements", 
         "are of interest.")
  if (any(grepl("(RR)|(RD)", EstimandType))) {
    if (length(Intervention) < 2)
      stop("Risk Ratios and Risk Differences can only be computed if Interventions is a vector of ", 
           "length 2, i.e. c('treated' ConcreteEst list index, 'control' ConcreteEst list index). Either", 
           "specify at least two indices in Intervention or remove 'RR' and 'RD' from Estimand.")
    if (length(Intervention) > 2)
      message("Risk ratios and risk differences will be computed using only the first two ",
          "elements of the provided 'Intervention' argument.\n", sep = "")
  }
  if (!is.logical(Simultaneous))
    stop("Simultaneous must be a logical, TRUE or FALSE")
  if (!is.numeric(Signif) | Signif >= 1 | Signif <= 0)
    stop("Signif sets the alpha (significance level) for hypothesis testing and confidence intervals", 
         ", so should be small number that must be greater than 0")
  
  Output <- data.table()
  Risks <- getRisk(ConcreteEst = ConcreteEst, TargetTime = TargetTime, TargetEvent = TargetEvent, GComp = GComp)
  if (any(sapply(Estimand, is.function))) {
    Output <- lapply(Estimand[which(sapply(Estimand, is.function))], function(estimand) {
      do.call(estimand, list("ConcreteEst" = ConcreteEst, "TargetTime" = TargetTime,
                             "TargetEvent" = TargetEvent, "GComp" = GComp))
    })
    names(Output) <- names(Estimand)[which(sapply(Estimand, is.function))]
  }
  
  if (any(grepl("Risk", EstimandType))) 
    Output <- rbind(Output, Risks)
  
  if (any(grepl("RD", EstimandType))) 
    Output <- rbind(Output, 
                    getRD(Risks = Risks, A1 = names(ConcreteEst)[Intervention[1]], 
                          A0 = names(ConcreteEst)[Intervention[2]], TargetTime = TargetTime, 
                          TargetEvent = TargetEvent, GComp = GComp))
  
  if (any(grepl("RR", EstimandType))) 
    Output <- rbind(Output, 
                    getRR(Risks = Risks, A1 = names(ConcreteEst)[Intervention[1]], 
                          A0 = names(ConcreteEst)[Intervention[2]], TargetTime = TargetTime, 
                          TargetEvent = TargetEvent, GComp = GComp))
  
  Output[, `CI Low` := `Pt Est` - qnorm(1 - Signif/2)*se]
  Output[, `CI Hi` := `Pt Est` + qnorm(1 - Signif/2)*se]
  
  if (Simultaneous)
    Output <- getSimultaneous(ConcreteEst = ConcreteEst, Output = Output, EstimandType = EstimandType, 
                              Intervention = Intervention, Signif = Signif)
  
  setorderv(Output, cols = c("Time", "Event", "Estimand", "Intervention", "Estimator"), 
            order = c(1, 1, 1, 1, -1))
  setcolorder(Output, c("Time", "Event", "Estimand", "Intervention", "Estimator"))
  
  attr(Output,"Estimand") <- EstimandType
  attr(Output,"GComp") <- GComp
  attr(Output,"Simultaneous") <- Simultaneous
  attr(Output,"Signif") <- Signif 
  class(Output) <- union("ConcreteOut", class(Output))
  return(Output)
}

getRisk <- function(ConcreteEst, TargetTime, TargetEvent, GComp) {
  Event <- Time <- seEIC <- Estimand <- NULL
  risk <- lapply(seq_along(ConcreteEst), function(a) {
    est.a <- ConcreteEst[[a]]
    risk.a <- do.call(rbind, lapply(as.character(TargetEvent), function(j) {
      risks <- apply(est.a[["Hazards"]][[j]] * est.a[["EvntFreeSurv"]], 2, cumsum)
      Psi <- cbind("Time" = attr(ConcreteEst, "Times")[which(attr(ConcreteEst, "Times") %in% TargetTime)],
                   "Risk" = rowMeans(subset(risks, subset = attr(ConcreteEst, "Times") %in% TargetTime)))
      se <- subset(est.a$SummEIC[Event == j, ], select = c("Time", "seEIC"))
      Psi <- merge(Psi, se[, list("Time" = Time, "se" = seEIC / sqrt(ncol(risks)))], by = "Time")
      return(cbind("Estimator" = "tmle", "Event" = as.numeric(j), Psi))
    }))
    if (GComp)
      risk.a <- rbind(risk.a,
                      cbind("Estimator" = "gcomp",
                            est.a$GCompEst[Event %in% TargetEvent, ],
                            "se" = NA))
    risk.a <- cbind(Intervention = names(ConcreteEst)[a], 
                    risk.a)
    return(risk.a)
  })
  Risk <- setDT(do.call(rbind, risk))
  Risk[, Estimand := "Abs Risk"]
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

getSimultaneous <- function(ConcreteEst, Output, EstimandType, Intervention, Signif) {
  Estimator <- `Pt Est` <- Event <- Time <- SimQ <- IC <- Risk <- NULL
  
  if (any(grepl("(RR)|(RD)", EstimandType))) {
    A1 <- names(ConcreteEst)[Intervention[1]]
    A0 <- names(ConcreteEst)[Intervention[2]]
  }
  
  ICs <- do.call(rbind, lapply(Intervention, function(a) {
    cbind("Intervention" = names(ConcreteEst)[a], ConcreteEst[[a]][["IC"]])
  }))
  
  ICs <- merge(ICs, getRisk(ConcreteEst = ConcreteEst, TargetTime = attr(ConcreteEst, "TargetTime"), 
                            TargetEvent = attr(ConcreteEst, "TargetEvent"), GComp = FALSE), 
               by = c("Intervention", "Time", "Event"))
  if (any(grepl("RD", EstimandType)))
    RDICs <- ICs[, list("Intervention" = "Risk Diff",
                        "IC" = IC[Intervention == A1] - IC[Intervention == A0]), 
                 by = c("ID", "Time", "Event")]
  if (any(grepl("RR", EstimandType)))
    RRICs <- ICs[, list("Intervention" = "Rel Risk",
                        "IC" = IC[Intervention == A1] / `Pt Est`[Intervention == A0] - 
                          IC[Intervention == A0] * `Pt Est`[Intervention == A1] / 
                          `Pt Est`[Intervention == A0]^2), 
                 by = c("ID", "Time", "Event")]
  if (any(grepl("Risk", EstimandType))) {
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
  n <- length(attr(ConcreteEst, "T.tilde"))
  
  q <- apply(abs(MASS::mvrnorm(n = 1e3, mu = rep(0, nrow(CorrEIC)), Sigma = CorrEIC)), 1, max)
  q <- as.numeric(stats::quantile(q, 1 - Signif))
  se <- data.table(names = rownames(CovEIC), SimQ = q)
  se[, c("Intervention", "Time", "Event") := tstrsplit(names, "_")]
  se[, Estimator := "tmle"][, Event := as.numeric(Event)][, Time := as.numeric(Time)]
  se[Intervention == "Rel Risk", Intervention := paste0("[", A1, "] / [", A0, "]")]
  se[Intervention == "Risk Diff", Intervention := paste0("[", A1, "] - [", A0, "]")]
  simCI <- merge(Output, se, c("Intervention", "Estimator", "Event", "Time"), all.x = TRUE)
  simCI[, "SimCI Low" := `Pt Est` - se*SimQ][, "SimCI Hi" := `Pt Est` + se*SimQ]
  simCI <- subset(simCI, select = c("Intervention", "Estimand", "Estimator", "Event", "Time", 
                                    "Pt Est", "se","CI Low", "CI Hi", "SimCI Low", "SimCI Hi"))
  return(simCI)
}

#' @describeIn doConcrete print.ConcreteOut print method for "ConcreteOut" class
#' @param x a ConcreteOut object
#' @param ... additional arguments to be passed into print methods
#' @exportS3Method print ConcreteOut
print.ConcreteOut <- function(x, ...) {
  num.cols <- c("Pt Est", "se", "CI Low", "CI Hi", "SimCI Low", "SimCI Hi")
  dt <- x[, (num.cols) := lapply(.SD, function(y) signif(y, 2)), .SDcols = num.cols]
  NextMethod(generic = "print", object = dt)
}


#' @describeIn getOutput plot.ConcreteOut plot method for "ConcreteOut" class
#' @param x a ConcreteOut object
#' @param NullLine logical: to plot a red line at y=1 for RR plots and at y=0 for RD plots
#' @param ask logical: to prompt for user input before each plot
#' @param ... additional arguments to be passed into plot methods
#' @exportS3Method plot ConcreteOut
plot.ConcreteOut <- function(x, NullLine = TRUE, ask = TRUE, ...) {
  `CI Low` <- `CI Hi` <- `SimCI Low` <- `SimCI Hi` <- Intervention <- 
    Event <- Time <- `Pt Est` <- Estimator <- se <- NULL
  Signif <- attr(x, "Signif")
  Estimand <- attr(x,"Estimand")
  GComp <- attr(x,"GComp")
  Simultaneous <- attr(x,"Simultaneous")
  
  if(!requireNamespace("ggplot2", quietly = TRUE))
    stop("Plotting requires the 'ggplot2' package")
  if (ask) {
    oask <- devAskNewPage(TRUE)
    on.exit(devAskNewPage(oask))
  }
  fig <- list()
  
  if (any(grepl("RR", Estimand))) {
    dev.hold()
    z <- x[Estimand == "Rel Risk", ][, Event := paste0("Event ", Event)]
    if (GComp) {
      fig[["rr"]] <- ggplot2::ggplot(z, ggplot2::aes(x = Time, y = `Pt Est`, shape = Estimator, 
                                                     ymin = `CI Low`, ymax = `CI Hi`)) + 
        ggplot2::geom_point() + ggplot2::geom_errorbar() + ggplot2::scale_shape_manual(values = c(4, 16))
    } else {
      fig[["rr"]] <- ggplot2::ggplot(z[Estimator == "tmle", ], 
                                     ggplot2::aes(x = Time, y = `Pt Est`, ymin = `CI Low`, ymax = `CI Hi`)) + 
        ggplot2::geom_point() + ggplot2::geom_errorbar()
    }
    if (Simultaneous & length(unique(z[["Time"]])) > 1) {
      fig[["rr"]] <- fig[["rr"]] + 
        ggplot2::geom_ribbon(data = z[Estimator == "tmle", ], 
                             ggplot2::aes(ymin = `SimCI Low`, ymax = `SimCI Hi`), alpha = 0.06)
    }
    fig[["rr"]] <- fig[["rr"]] + 
      ggplot2::facet_wrap(~Event, scales = "free", nrow = 1) + 
      ggplot2::theme_minimal() +
      ggplot2::labs(y = "Relative Risk", x = "Time", 
                    title = paste0("Relative Risk Point Estimates with ", 
                                   round(100 * (1 - Signif), digits = 0), "% confidence intervals"))
    if (NullLine)
      fig[["rr"]] <- fig[["rr"]] + ggplot2::geom_hline(ggplot2::aes(yintercept = 1), colour = "red", alpha = 0.4)
    plot(fig[["rr"]])
    dev.flush()
  }
  if (any(grepl("RD", Estimand))) {
    dev.hold()
    z <- x[Estimand == "Risk Diff", ][, Event := paste0("Event ", Event)]
    if (GComp) {
      fig[["rd"]] <- ggplot2::ggplot(z, ggplot2::aes(x = Time, y = `Pt Est`, shape = Estimator, 
                                                     ymin = `CI Low`, ymax = `CI Hi`)) + 
        ggplot2::geom_point() + ggplot2::geom_errorbar() + ggplot2::scale_shape_manual(values = c(4, 16))
    } else {
      fig[["rd"]] <- ggplot2::ggplot(z[Estimator == "tmle", ], ggplot2::aes(x = Time, y = `Pt Est`, 
                                                                            ymin = `CI Low`, ymax = `CI Hi`)) + 
        ggplot2::geom_point() + ggplot2::geom_errorbar()
    }
    if (Simultaneous & length(unique(z[["Time"]])) > 1) {
      fig[["rd"]] <- fig[["rd"]] + 
        ggplot2::geom_ribbon(data = z[Estimator == "tmle", ], 
                             ggplot2::aes(ymin = `SimCI Low`, ymax = `SimCI Hi`), alpha = 0.06)
    }
    fig[["rd"]] <- fig[["rd"]] + 
      ggplot2::facet_wrap(~Event, scales = "free", nrow = 1) + 
      ggplot2::theme_minimal() +
      ggplot2::labs(y = "Risk Difference", x = "Time", 
                    title = paste0("Risk Difference Point Estimates with ", 
                                   round(100 * (1 - Signif), digits = 0), "% confidence intervals")) 
    if (NullLine)
      fig[["rd"]] <- fig[["rd"]] + ggplot2::geom_hline(ggplot2::aes(yintercept = 0), colour = "red", alpha = 0.4)
    plot(fig[["rd"]])
    dev.flush()
  }
  if (any(grepl("Risk", Estimand))) {
    dev.hold()
    z <- x[Estimand == "Abs Risk", ][, Event := paste0("Event ", Event)]
    if (GComp) {
      fig[["risk"]] <- ggplot2::ggplot(z, ggplot2::aes(x = Time, y = `Pt Est`, shape = Estimator, 
                                                       colour = Intervention, ymin = `CI Low`, ymax = `CI Hi`)) + 
        ggplot2::geom_point() + ggplot2::geom_errorbar(alpha = 0.5) + 
        ggplot2::scale_shape_manual(values = c(4, 16))
    } else {
      fig[["risk"]] <- ggplot2::ggplot(z[Estimator == "tmle", ], 
                                       ggplot2::aes(x = Time, y = `Pt Est`, colour = Intervention, 
                                                    ymin = `CI Low`, ymax = `CI Hi`)) + 
        ggplot2::geom_point() + ggplot2::geom_errorbar(alpha = 0.5)
    }
    if (Simultaneous & length(unique(z[["Time"]])) > 1) {
      fig[["risk"]] <- fig[["risk"]] + 
        ggplot2::geom_ribbon(ggplot2::aes(ymin = `SimCI Low`, ymax = `SimCI Hi`, fill = Intervention), 
                             data = z[Estimator == "tmle", ], alpha = 0.06)
    }
    fig[["risk"]] <- fig[["risk"]] + 
      ggplot2::facet_wrap(~Event, scales = "free", nrow = 1) + 
      ggplot2::theme_minimal() +
      ggplot2::labs(y = "Absolute Risk", x = "Time", 
                    title = paste0("Absolute Risk Point Estimates with ", 
                                   round(100 * (1 - Signif), digits = 0), "% confidence intervals")) 
    plot(fig[["risk"]])
    dev.flush()
  }
  if (any(grepl("Function", Estimand))) {
    warning("Plotting for non-default estimands is not currently available.")
  }
  invisible(fig)
}
