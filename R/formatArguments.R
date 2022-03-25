
#' Title
#' @param EventTime : Numeric vector (N x 1)
#' @param EventType : Numeric (integer) vector (N x 1)
#' @param Treatment : Numeric vector (N x 1)
#' @param Intervention : list of function (length = A*)
#' @param CovDataTable : data.table (N x ?)
#' @param ID : vector (N x 1)
# #' @param LongTime : numeric vector (?? x 1)
#' @param TargetTimes : numeric vector (length = K)
#' @param TargetEvents : numeric vector \\subset EventType (length = J)
#' @param Models : list of functions (length = L)
#' @param CVArgs : list
#' @param NumUpdateSteps : numeric
#' @param OneStepEps : numeric
#' @param MinNuisance : numeric
#' @param PropScoreBackend : character
#' @param Verbose : boolean
#' @param GComp : boolean
#'
#' @return tbd
#' @export
#'

formatArguments <- function(EventTime, EventType, Treatment, Intervention, CovDataTable, ID, TargetTimes,
                            TargetEvents, Models, CVArgs, NumUpdateSteps, OneStepEps, MinNuisance,
                            PropScoreBackend, Verbose, GComp)
{
    # For user input: data + formula ----
    # try use prodlim::EventHistory.frame
    # x=EventHistory.frame(Hist(time,Status,cens.code="censored")~age+sex+intervention(trt)+stage,data=pbc,specials="intervention")
    # names(x)

    # OBS: learners can either work with the original data or with the design matrix (dummies)

    Models <-

    ## observed data ----
    if (!data.table::is.data.table(CovDataTable)) {
        CovDataTable <- data.table::as.data.table(CovDataTable)
        warning("CovDataTable must be a data.table. We have attempted to convert ",
                "the CovDataTable argument into an object of data.table class.\n")
    }

    if (any(c("ID",  "Event", "Trt", "t") %in% colnames(CovDataTable)))
        stop("'ID', 'Event', 'Trt', and 't' are reserved column",
             "names. Rename covariate column names to avoid name collisions.\n")


    if (!all(sapply(list(EventTime, EventType, Treatment, TargetTimes, TargetEvents),
                    function(vec) is.numeric(vec) | is.null(vec))))
        stop("EventTime, EventType, and Treatment ",
             "arguments must be numeric vectors\n")

    if (is.null(ID))
        ID = seq_along(EventTime)

    Data <- try(
        data.table::data.table("ID" = ID,
                               "Time" = EventTime,
                               "Event" = EventType,
                               "Trt" = Treatment,
                               CovDataTable)
    )
    if ("try-error" %in% class(Data)) {
        warning("Failed to create data datatable. ",
                "Check data inputs; see function help page\n")
        return(Data)
    }

    ReservedColumns <- c('Time', 'Event', 'Trt')

    Events <- sort(unique(Data$Event))
    Censored <- 0 %in% Events
    if (length(Models) != length(Events) + 1)
        stop("Models must be provided for treatment and every observed event or censoring type")
    Events <- Events[Events > 0]

    ## target event(s) ----
    if (is.null(TargetEvents))
        TargetEvents <- Events
    else if (!all(TargetEvents %in% Events))
        stop("Some TargetEvents are not present in the data")

    ## target time(s) ----
    if (max(TargetTimes) >= max(Data["Event" != 0, "Time"])) {
        TargetTimes <- TargetTimes[TargetTimes < max(Data[["Time"]])]
        warning(paste0("No Observed events at max target time:",
                       " truncating target time(s) to be at or ",
                       "before the last observed event time"))
    }

    ## regimes of interest ----
    if (is.list(Intervention)) {
        RegsOfInterest <- lapply(Intervention, function(intervene) {
            if (is.function(intervene)) {
                Regime <- do.call(intervene, list(Treatment, CovDataTable))
                if (is.null(attr(Regime, "g.star"))) {
                    attr(Regime, "g.star") <- function(a) as.numeric(a == Regime)
                    warning("no g.star input, defaulting to the indicator that observed Treatment == desired RegName")
                }
                return(Regime)
            }
            else stop("Intervention must be a list of functions. See doConCRTmle documentation")
        })
    }

    # To do:
    # check if covariates too highly correlated with ID, Time, Event, or Trt
    # target event(s)
    # binary treatment
    # stopping criteria
    # ...


    return(list(Data = Data, Events = Events, RegsOfInterest = RegsOfInterest, Censored = Censored))
}


checkEventTimes <- function(x) {
    if (any(!is.vector(x), !is.numeric(x), !try(x > 0), is.infinite(x), is.list(x)))
        stop("EventTimes must be a numeric vector with positive, finite values")
}

checkEventTypes <- function(x) {
    if (any(!is.vector(x), !is.numeric(x), !try(x >= 0), is.list(x)))
        stop("EventTypes must be a numeric vector with non-negative values")
}

checkTreatment <- function(x) {
    if (any(!is.vector(x), !is.numeric(x), is.nan(x), is.infinite(x), is.list(x)))
        stop("Treatment must be a numeric vector with finite values")
}

checkCovDataTable <- function(x) {
    datatypes <- c("matrix", "data.frame", "data.table")
    if (!any(sapply(datatypes, function(class) inherits(x, class))))
        stop("CovDataTable must be a matrix, data.frame, or data.table")
    if (any(is.infinite(unlist(x))))
        warning("CovDataTable contains infinite values; regression models may break")
}

checkID <- function(x) {
    if (any(!is.vector(x), is.list(x), is.null(x), is.nan(x), is.infinite(x), is.na(x)))
        stop("ID must be a vector")
}

checkIntervention <- function(x) {
    if (!is.list(x) | is.null(names(x)))
        stop("Intervention must be a named list of intervention functions")
}

# RegOfInt <- setClass(Class = "RegOfInt", slots = c(intervention = "function", g.star = "function"))
# setMethod(f = "initialize", signature = "RegOfInt",
#           definition = function(.Object,
#                                 intervention = function(a.obs, L.obs) return(rep_len(1, length(a.obs))),
#                                 g.star = function(a, L) return(a == 1), ...) {
#               .Object <- callNextMethod(.Object, ...)
#               .Object@intervention <- intervention
#               .Object@g.star <- g.star
#               .Object
#           })

