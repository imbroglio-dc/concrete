
#' doConcrete
#'
#' @param ConcreteArgs : output from concrete::formatArguments
#'
# #' @param Data : data.table (N x ?)
# #' @param CovDataTable : data.table (N x ?)
# #' @param ID : vector (N x 1)
# #' @param LongTime : numeric vector (?? x 1)
# #' @param TargetTime : numeric vector (length = K)
# #' @param TargetEvent : numeric vector \\subset EventType (length = J)
# #' @param Model : list of functions (length = L)
# #' @param CVFolds : list
# #' @param MaxUpdateIter : numeric
# #' @param OneStepEps : numeric
# #' @param MinNuisance : numeric
# #' @param Verbose : boolean
# #' @param PropScoreBackend : character
# #' @param GComp : boolean
# #' @param Events : numeric
# #' @param Censored : boolean
# #' @param Regime : list
#'
#' @import data.table
#'
#' @return tbd
#' @export doConcrete
#'
#' @examples
#' "tbd"

doConcrete <- function(ConcreteArgs) {
  return(do.call(doConCRTmle, ConcreteArgs))
}

doConCRTmle <- function(Data, CovDataTable,
                        LongTime, ID, Events, Censored,
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

  # get initial EIC (possibly with GComp Estimate) ---------------------------------------------
  Estimates <- getEIC(Estimates = Estimates, Data = Data, Regime = Regime, Censored = Censored,
                      TargetEvent = TargetEvent, TargetTime = TargetTime, Events = Events,
                      MinNuisance = MinNuisance, GComp = GComp)

  SummEIC <- do.call(rbind, lapply(seq_along(Estimates), function(a) {
    cbind("Trt" = names(Estimates)[a], Estimates[[a]][["SummEIC"]])}))
  NormPnEIC <- getNormPnEIC(SummEIC[Time %in% TargetTime & Event %in% TargetEvent, PnEIC])

  ## initial estimator (g-computation) --------------------------------------------------------
  if (GComp)
    GCompEst <- do.call(rbind, lapply(seq_along(Estimates), function(a) {
      cbind("Trt" = names(Estimates)[a], Estimates[[a]][["GCompEst"]])}))

  # Update step -------------------------------------------------------------------------------
  ## Check if EIC is solved sufficienty and return outputs ----
  ## check PnEIC <= seEIC / (sqrt(n) log(n))
  OneStepStop <- SummEIC[, list("check" = abs(PnEIC) <= `seEIC/(sqrt(n)log(n))`,
                                "ratio" = abs(PnEIC) / `seEIC/(sqrt(n)log(n))`),
                         by = c("Trt", "Time", "Event")]
  if (Verbose)
    print(OneStepStop[["ratio"]])
  if (!all(sapply(OneStepStop[["check"]], isTRUE))) {
    ## one-step tmle loop (one-step) ----
    Estimates <- doTmleUpdate(Estimates = Estimates, SummEIC = SummEIC, Data = Data,
                              Censored = Censored, TargetEvent = TargetEvent,
                              TargetTime = TargetTime, Events = Events,
                              MaxUpdateIter = MaxUpdateIter, OneStepEps = OneStepEps,
                              NormPnEIC = NormPnEIC, Verbose = Verbose)
  }

  # format output --------------------------------------------------------------------------------------

  SummEIC <- do.call(rbind, lapply(seq_along(Estimates), function(a) {
    cbind("Trt" = names(Estimates)[a], Estimates[[a]][["SummEIC"]])}))

  NormPnEIC <- getNormPnEIC(SummEIC[Time %in% TargetTime & Event %in% TargetEvent, PnEIC])

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
