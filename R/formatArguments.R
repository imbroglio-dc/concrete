
## check argument formats, types, and sizes ----
formatArguments <- function(EventTime, EventType, Treatment, CovDataTable, ID, TargetTimes,
                            TargetEvents, CVArgs, NumUpdateSteps, OneStepEps, MinNuisance,
                            PropScoreBackend, Verbose, GComp)
{
    # For user input: data + formula ----
    # try use prodlim::EventHistory.frame
    # x=EventHistory.frame(Hist(time,Status,cens.code="censored")~age+sex+intervention(trt)+stage,data=pbc,specials="intervention")
    # names(x)

    # OBS: learners can either work with the original data or with the design matrix (dummies)

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
