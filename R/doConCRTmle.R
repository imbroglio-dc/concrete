#' doConcrete
#'
#' @param ConcreteArgs : output from concrete::formatArguments
#'
# #' @param Data : data.table (N x ?)
# #' @param CovDataTable : data.table (N x ?)
# #' @param LongTime : numeric vector (?? x 1)
# #' @param ID : vector (N x 1)
# #' @param Events : numeric
# #' @param Censored : boolean
# #' @param TargetTime : numeric vector (length = K)
# #' @param TargetEvent : numeric vector \\subset EventType (length = J)
# #' @param Regime : list
# #' @param CVFolds : list
# #' @param Model : list of functions (length = L)
# #' @param PropScoreBackend : character
# #' @param HazEstBackend : character
# #' @param MaxUpdateIter : numeric
# #' @param OneStepEps : numeric
# #' @param MinNuisance : numeric
# #' @param Verbose : boolean
# #' @param GComp : boolean
#'
#' @import data.table
#'
#' @return tbd
#'
#' @export doConcrete
#'
#' @examples
#' library(data.table)
#' library(survival)
#' library(concrete)
#' data <- as.data.table(survival::pbc)
#' data[, trt := sample(0:1, nrow(data), TRUE)]
#' cols <- c("id", "time", "status", "trt",
#'           "age", "albumin", "sex", "bili")
#' data <- data[, .SD, .SDcols = cols]
#' 
#' intervention <- concrete:::ITT
#' target.time <- 2500
#' target.event <- 1:2
#' model <- list("trt" = c("SL.glm", "SL.glmnet"),
#'               "0" = list(Surv(time, status == 0) ~ .),
#'               "1" = list(Surv(time, status == 1) ~ .),
#'               "2" = list(Surv(time, status == 2) ~ .))
#' 
#' # formatArguments() returns correctly formatted arguments for doConcrete()
#' concrete.args <- formatArguments(DataTable = data,
#'                                  EventTime = "time",
#'                                  EventType = "status",
#'                                  Treatment = "trt",
#'                                  ID = "id",
#'                                  Intervention = intervention,
#'                                  TargetTime = target.time,
#'                                  TargetEvent = target.event,
#'                                  Model = model, Verbose = FALSE)
#' 
#' # doConcrete() returns tmle (and g-comp plug-in) estimates of targeted risks
#' concrete.est <- doConcrete(concrete.args)
#' 
#' # getOutput returns risk difference, relative risk, and treatment-specific risks
#' concrete.out <- getOutput(Estimate = concrete.est, Estimand = c("rd", "rr", "risk"), 
#'                           TargetTime = target.time, TargetEvent = target.event, GComp = TRUE)
#' concrete.out$RD
#' concrete.out$RR
#' concrete.out$Risk

doConcrete <- function(ConcreteArgs) {
  return(do.call(doConCRTmle, ConcreteArgs))
}

doConCRTmle <- function(Data, CovDataTable, LongTime, ID, Events, Censored,
                        TargetTime, TargetEvent, Regime,
                        CVFolds, Model, PropScoreBackend, HazEstBackend,
                        MaxUpdateIter, OneStepEps, MinNuisance,
                        Verbose, GComp)
{
  Time <- Event <- PnEIC <- `seEIC/(sqrt(n)log(n))` <- NULL # for data.table compatibility w/ global var binding check

  # initial estimation ------------------------------------------------------------------------
  Estimates <- getInitialEstimate(Data = Data, CovDataTable = CovDataTable,
                                  Model = Model, CVFolds = CVFolds, MinNuisance = MinNuisance,
                                  TargetEvent = TargetEvent, TargetTime = TargetTime, Regime = Regime,
                                  PropScoreBackend = PropScoreBackend, HazEstBackend = HazEstBackend,
                                  Censored = Censored)

  # get initial EIC (possibly with GComp plug-in estimate) ---------------------------------------------
  Estimates <- getEIC(Estimates = Estimates, Data = Data, Regime = Regime, Censored = Censored,
                      TargetEvent = TargetEvent, TargetTime = TargetTime, Events = Events,
                      MinNuisance = MinNuisance, GComp = GComp)


  # Update step -------------------------------------------------------------------------------
  ## Check if EIC is solved sufficienty and return outputs ----
  ## check PnEIC <= seEIC / (sqrt(n) log(n))
  SummEIC <- do.call(rbind, lapply(seq_along(Estimates), function(a) {
    cbind("Trt" = names(Estimates)[a], Estimates[[a]][["SummEIC"]])}))
  NormPnEIC <- getNormPnEIC(SummEIC[Time %in% TargetTime & Event %in% TargetEvent, PnEIC])
  OneStepStop <- SummEIC[, list("check" = abs(PnEIC) <= `seEIC/(sqrt(n)log(n))`,
                                "ratio" = abs(PnEIC) / `seEIC/(sqrt(n)log(n))`),
                         by = c("Trt", "Time", "Event")]
  if (Verbose)
    print(OneStepStop[["ratio"]])

  ## one-step tmle loop (one-step) ----
  if (!all(sapply(OneStepStop[["check"]], isTRUE))) {
    Estimates <- doTmleUpdate(Estimates = Estimates, SummEIC = SummEIC, Data = Data,
                              Censored = Censored, TargetEvent = TargetEvent,
                              TargetTime = TargetTime, Events = Events,
                              MaxUpdateIter = MaxUpdateIter, OneStepEps = OneStepEps,
                              NormPnEIC = NormPnEIC, Verbose = Verbose)
  }

  # format output --------------------------------------------------------------------------------------

  # g-comp (sl estimate)
  # unadjusted cox model
  # tmle & ic

  return(Estimates)
}

getNormPnEIC <- function(PnEIC, Sigma = NULL) {
  WeightedPnEIC <- PnEIC
  if (!is.null(Sigma)) {
    SigmaInv <- try(solve(Sigma))
    if (any(class(SigmaInv) == "try-error")) {
      SigmaInv <- solve(Sigma + diag(x = 1e-6, nrow = nrow(Sigma)))
      warning("regularization of Sigma needed for inversion")
    }
    WeightedPnEIC <- PnEIC %*% SigmaInv
  }
  return(sqrt(sum(unlist(PnEIC) * unlist(WeightedPnEIC))))
}
