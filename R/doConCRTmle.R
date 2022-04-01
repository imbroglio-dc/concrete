
#' doConCRTmle
#'
#' @param EventTime : Numeric vector (N x 1)
#' @param EventType : Numeric (integer) vector (N x 1)
#' @param Treatment : Numeric vector (N x 1)
#' @param Intervention : list of function (length = A*)
#' @param CovDataTable : data.table (N x ?)
#' @param ID : vector (N x 1)
#' @param LongTime : numeric vector (?? x 1)
#' @param TargetTimes : numeric vector (length = K)
#' @param TargetEvents : numeric vector \\subset EventType (length = J)
#' @param Models : list of functions (length = L)
#' @param CVArgs : list
#' @param NumUpdateSteps : numeric
#' @param OneStepEps : numeric
#' @param MinNuisance : numeric
#' @param Verbose : boolean
#' @param PropScoreBackend : character
#' @param GComp : boolean
#'
#' @import data.table
#'
#' @return tbd
#' @export
#'
#' @examples
#' "tbd"

# To Do : return warning if targeting time before any observed events of type J


doConCRTmle <- function(EventTime, EventType, Treatment, Intervention, CovDataTable, LongTime = NULL,
                        ID = NULL, TargetTimes = sort(unique(EventTime)),
                        TargetEvents = NULL, Models, CVArgs = NULL, NumUpdateSteps = 25,
                        OneStepEps = 0.1, MinNuisance = 0.05, PropScoreBackend = "sl3",
                        Verbose = FALSE, GComp = FALSE)
{
  Time <- Event <- PnEIC <- `seEIC/(root(n)log(n))` <- NULL # for R CMD globar variable binding check
  Args <- formatArguments(EventTime, EventType, Treatment, Intervention, CovDataTable, ID, TargetTimes,
                          TargetEvents, Models, CVArgs, NumUpdateSteps, OneStepEps, MinNuisance,
                          PropScoreBackend, Verbose, GComp)
  Data <- Args[["Data"]]
  RegsOfInterest <- Args[["RegsOfInterest"]]
  Events <- Args[["Events"]]
  Censored <- Args[["Censored"]]

  # initial estimation ------------------------------------------------------------------------
  Estimates <- getInitialEstimate(Data, CovDataTable, Models, MinNuisance, TargetEvents,
                                   TargetTimes, RegsOfInterest, PropScoreBackend, Censored)

  # get initial EIC (possibly with GComp Estimate) ---------------------------------------------
  Estimates <- getEIC(Estimates, Data, RegsOfInterest, Censored, TargetEvents,
                      TargetTimes, Events, MinNuisance, GComp)

  SummEIC <- do.call(rbind, lapply(seq_along(Estimates), function(a) {
    cbind("Trt" = names(Estimates)[a], Estimates[[a]][["SummEIC"]])}))
  NormPnEIC <- getNormPnEIC(SummEIC[Time %in% TargetTimes & Event %in% TargetEvents, PnEIC])

  ## initial estimator (g-computation) --------------------------------------------------------
  if (GComp)
    GCompEst <- do.call(rbind, lapply(seq_along(Estimates), function(a) {
      cbind("Trt" = names(Estimates)[a], Estimates[[a]][["GCompEst"]])}))

  # Update step -------------------------------------------------------------------------------
  ## Check if EIC is solved sufficienty and return outputs ----
  ## check PnEIC <= seEIC / (sqrt(n) log(n))
  OneStepStop <- SummEIC[, list("check" = abs(PnEIC) <= `seEIC/(root(n)log(n))`,
                                "ratio" = abs(PnEIC) / `seEIC/(root(n)log(n))`),
                         by = c("Trt", "Time", "Event")]
  if (Verbose)
    print(OneStepStop[["ratio"]])
  if (!all(sapply(OneStepStop[["check"]], isTRUE))) {
    ## one-step tmle loop (one-step) ----
    Estimates <- doTmleUpdate(Estimates, SummEIC, Data, Censored, TargetEvents, TargetTimes, Events,
                              NumUpdateSteps, OneStepEps, NormPnEIC, Verbose)
  }

  # format output --------------------------------------------------------------------------------------

  SummEIC <- do.call(rbind, lapply(seq_along(Estimates), function(a) {
    cbind("Trt" = names(Estimates)[a], Estimates[[a]][["SummEIC"]])}))

  NormPnEIC <- getNormPnEIC(SummEIC[Time %in% TargetTimes & Event %in% TargetEvents, PnEIC])


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
