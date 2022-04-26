
#' Title
#' @param DataTable data.table (n x (d + (3:4)); data.table of the observed data. Must include:
#' \itemize{
#'   \item{"EventTime"}{: non-negative real numbers; the observed event or censoring time}
#'   \item{"EventType"}{: numeric; the observed event type (0 for censoring)}
#'   \item{"Treatment"}{: numeric; the observed treatment}
#' }
#' May include
#' \itemize{
#'   \item{"ID"}{: factor, character, or numberic; subject id - if missing, ID will be
#'   set to row number. For longitudinal data, ID must be provided}
#'   \item{"LongTime"}{: numeric; (Specifies monitoring times for longitudinal data structures}
#'   \item{"Treatment"}{: numeric}
#' }
#' @param DataStructure (not implemented yet) formula: e.g. Surv("EventTime", "EventType") ~ Intervention("Treatment") + ...
#' @param EventTime character: the column name of the observed time
#' @param EventType character: the column name of the observed event type
#' @param Treatment character: the column name of the observed treatment assignment
#' @param ID (optional) character: the column name of the observed subject id
#' @param LongTime (situational) character: the column name of the monitoring times for
#' longitudinal data structures
#' @param Intervention list: a list of desired interventions on the treatment variable.
#' Each intervention must be a list containing two named functions: 'intervention' = function(treatment vector, covariate data)
#' and 'gstar' = function(treatment vector, covariate data)
#' @param TargetTime numeric vector (length = K)
#' @param TargetEvent numeric vector that is a subset of EventType (length = J)
#' @param Target (not yet implemented) data.table / data.frame (?? x 2); a table containing all combinations of target events (column 1) and target times (column 2).
#' @param CVArg list of arguments to be passed into origami::make_folds. the default = list(n = nrow(DataTable), fold_fun = folds_vfold, cluster_ids = NULL, strata_ids = NULL)
#' @param Model list of models (length = L + Censoring + Treatment)
#' @param PropScoreBackend (currently must be `sl3`) character
#' @param HazEstBackend (currently must be `coxph`) character
#' @param MaxUpdateIter numeric: the number of one-step update steps
#' @param OneStepEps numeric: the one-step tmle step size
#' @param MinNuisance numeric: the minimum value of the nuisance parameter denominator in the clever covariate
#' @param Verbose boolean
#' @param GComp boolean
#'
#' @return tbd
#'
#' @importFrom stats model.matrix
#' @import origami
#' @export formatArguments
#'
#'

formatArguments <- function(DataTable, DataStructure = NULL, EventTime, EventType, Treatment, ID = NULL, LongTime = NULL,
                            Intervention, TargetTime, TargetEvent, Target = NULL,
                            CVArg = list(n = nrow(DataTable), fold_fun = folds_vfold,
                                         cluster_ids = NULL, strata_ids = NULL),
                            Model, PropScoreBackend = "sl3", HazEstBackend = "coxph",
                            MaxUpdateIter = 100, OneStepEps = 0.1, MinNuisance = 0.05,
                            Verbose = TRUE, GComp = TRUE)
{
    ## Data Structure ----
    # incorporate prodlim::EventHistory.frame?
    data.tbl <- formatDataTable(x = DataTable, EventTime = EventTime, EventType = EventType,
                                Treatment = Treatment, ID = ID, LongTime = LongTime)
    event.time <- getEventTime(x = EventTime, DataTable = data.tbl)
    event.type <- getEventType(x = EventType, DataTable = data.tbl)
    treatment <- getTreatment(Treatment, DataTable = data.tbl)
    id <- getID(ID, DataTable = DataTable)
    long.time <- getLongTime(x = LongTime, DataTable = data.tbl)
    cov.data.table <- getCovDataTable(DataTable = data.tbl, EventTime = EventTime,
                                      EventType = EventType, Treatment = Treatment, ID = ID,
                                      LongTime = LongTime, Verbose = Verbose)
    censored <- 0 %in% event.type
    event.type.unique <- setdiff(sort(unique(event.type)), 0)

    ## Interventions & Targets ----
    regime <- getRegime(Intervention = Intervention, Treatment = treatment,
                        CovDataTable = cov.data.table)

    target.event <- getTargetEvent(TargetEvent = TargetEvent, UniqueEvents = event.type.unique)
    checkTargetTime(TargetTime = TargetTime, EventTime = event.time, TargetEvent = TargetEvent,
                    EventType = event.type)

    ## Estimation Paramters ----
    cv.folds <- getCVFolds(CVArg)
    checkModel(Model = Model, UniqueEvent = event.type.unique,
               Censored = censored, HazEstBackend = HazEstBackend)
    checkPropScoreBackend(PropScoreBackend)
    checkHazEstBackend(HazEstBackend)

    ## TMLE Update Parameters ----
    max.update.iter <- getMaxUpdateIter(MaxUpdateIter)
    checkOneStepEps(OneStepEps)

    # min.nuisance <- getMinNuisance()

    ## Misc. Parameters
    checkVerbose(Verbose)

    checkGComp(GComp)

    attr(data.tbl, "EventTime") <- EventTime
    attr(data.tbl, "EventType") <- EventType
    attr(data.tbl, "Treatment") <- Treatment
    attr(data.tbl, "LongTime") <- LongTime
    attr(data.tbl, "ID") <- ID

    return(list(Data = data.tbl,
                CovDataTable = cov.data.table,
                LongTime = long.time,
                ID = id,
                Events = event.type.unique,
                Censored = censored,
                TargetTime = TargetTime,
                TargetEvent = target.event,
                Regime = regime,
                CVFolds = cv.folds,
                Model = Model,
                PropScoreBackend = PropScoreBackend,
                HazEstBackend = HazEstBackend,
                MaxUpdateIter = max.update.iter,
                OneStepEps = OneStepEps,
                MinNuisance = MinNuisance,
                Verbose = Verbose,
                GComp = GComp))
}

formatDataTable <- function(x, EventTime, EventType, Treatment, ID, LongTime) {
    data_types <- c("matrix", "data.frame", "data.table")
    if (!any(sapply(data_types, function(class) inherits(x, class))))
        stop("CovDataTable must be a data.table or a matrix / data.frame coercible into data.table.")
    if (!inherits(x, "data.table"))
        x <- data.table::as.data.table(x)
    if (any(is.infinite(unlist(x)), anyNA(unlist(x))))
        stop("CovDataTable contains infinite or missing values; regression models may break")
    return(x)
}

getEventTime <- function(x, DataTable = NULL) {
    if (is.character(x)) {
        tmp <- try(DataTable[[x]])
        names(tmp) <- x
        x <- tmp
    }
    if (any(!is.vector(x), !is.numeric(x), try(x <= 0), inherits(try(x <= 0), "try-error"),
            is.infinite(x), is.list(x)))
        stop("EventTime must be a numeric vector with positive, finite values")
    return(x)
}

getEventType <- function(x, DataTable = NULL) {
    if (is.character(x)) {
        tmp <- try(DataTable[[x]])
        names(tmp) <- x
        x <- tmp
    }
    if (any(!is.vector(x), !is.numeric(x), try(x < 0), inherits(try(x < 0), "try-error"), is.list(x)))
        stop("EventType must be a numeric vector with non-negative values (0 indicating censoring)")
    return(x)
}

getTreatment <- function(x, DataTable = NULL) {
    if (is.character(x)) {
        tmp <- try(DataTable[[x]])
        names(tmp) <- x
        x <- tmp
    }
    if (any(!is.vector(x), !is.numeric(x), is.nan(x), is.infinite(x), is.list(x)))
        stop("Treatment must be a numeric vector with finite values")
    return(x)
}

getID <- function(x, DataTable = NULL) {
    if (is.null(x)) {
        x <- 1:nrow(DataTable)
        message("No ID column specified. DataTable row numbers will be used as subject IDs, ",
                "which will not be appropriate for longitudinal data structures.")
    }
    if (is.character(x)) {
        tmp <- try(DataTable[[x]])
        names(tmp) <- x
        x <- tmp
    }
    if (any(!is.vector(x), is.list(x), is.null(x), is.nan(x), is.na(x)))
        stop("ID column must not include missing values")
    return(x)
}

getLongTime <- function(x, DataTable = NULL) {
    warning("checks not yet implemented for longitudinal monitoring time, LongTime")
    return(NULL)
}

getCovDataTable <- function(DataTable, EventTime, EventType, Treatment, ID, LongTime, Verbose) {
    cov.names <- setdiff(colnames(DataTable), c(EventTime, EventType, Treatment, ID, LongTime))
    cov.dt <- DataTable[, .SD ,.SDcols = cov.names]
    non.num.ind <- sapply(cov.names, function(cov.name) {!inherits(cov.dt[[cov.name]], c("numeric", "integer"))})
    non.num.covs <- cov.names[non.num.ind]
    if (length(cov.names[!non.num.ind]) == 0) {
        cov.dt.mod.matrixed <- data.table::data.table("dummy" = rep_len(1, nrow(cov.dt)))
    } else {
        cov.dt.mod.matrixed <- cov.dt[, .SD, .SDcols = cov.names[!non.num.ind]]
    }
    cov.names.mod.matrixed <- vector("list", length = length(non.num.covs))
    for (cov.name in non.num.covs) {
        cov.mod.matrixed <- data.table::as.data.table(model.matrix(~-1 + ., subset(cov.dt, select = cov.name)))

        cov.dt.mod.matrixed <- cbind(cov.dt.mod.matrixed, cov.mod.matrixed)
        cov.names.mod.matrixed[[cov.name]] <- colnames(cov.mod.matrixed)
    }
    if (length(cov.names[!non.num.ind]) == 0) {
        cov.dt.mod.matrixed <- cov.dt.mod.matrixed[, .SD,
                                                   .SDcols = setdiff(colnames(cov.dt.mod.matrixed),
                                                                     "dummy")]
    }
    attr(cov.dt.mod.matrixed, "cov.names.mod.matrixed") <- cov.names.mod.matrixed
    if (Verbose)
        # try(superheat::superheat(cov(scale(model.matrix(~., cov.dt.mod.matrixed)))))
        return(cov.dt.mod.matrixed)
}

getRegime <- function(Intervention, Treatment, CovDataTable) {
    if (is.list(Intervention)) {
        Regime <- lapply(seq_along(Intervention), function(i) {
            regime <- Intervention[[i]]
            reg_name <- names(Intervention)[i]
            if (is.function(regime$intervention)) {
                Regime <- try(do.call(regime$intervention, list(Treatment, CovDataTable)))
                if (inherits(Regime, "try-error"))
                    stop("Intervention must be a list of regimes specificed as list(intervention",
                         " = f(A, L), g.star = g(A, L)), and intervention function f(A, L) must ",
                         "be a function of treatment and covariates that returns numeric desired",
                         " treatment assignments (a*) with the same dimensions as the observed ",
                         "treatment. Amend Intervention[[", reg_name, "]] and try again")
                if (any(length(unlist(Regime)) != length(unlist(Treatment)), !is.numeric(unlist(Regime))))
                    stop("The intervention functions f(A, L) must be functions of treatment and ",
                         "covariates that return numeric desired treatment assignments (a*) with ",
                         "the same dimensions as the observed treatment. Amend Intervention[[",
                         reg_name, "]] and try again")
                if (mean(unlist(Regime) == unlist(Treatment)) < 0.05)
                    warning("The intervention function f(A, L) specified in Intervention[[", reg_name,
                            "]] is likely to cause near-positivity issues, resulting in inflated ",
                            "variance and perhaps instability. Highly recommend reconsidering this",
                            "choice of intervention.")
                if (is.null(regime$g.star)) {
                    attr(Regime, "g.star") <- function(a) as.numeric(a == Regime)
                    message("No g.star function specified, defaulting to the indicator that observed",
                            "treatment equals the intervention output treatment level, 1(A = a*).")
                } else {
                    attr(Regime, "g.star") <- regime$g.star
                    testRegGStar <- try(do.call(regime$g.star, list(Treatment, CovDataTable)))
                    if (inherits(testRegGStar, "try-error"))
                        stop("Intervention must be a list of regimes specificed as list(intervention",
                             " = f(A, L), g.star = g(A, L)), and the g.star function f(A, L) must ",
                             "be a function of treatment and covariates that returns numeric ",
                             "probabilities bounded in [0, 1]. Amend Intervention[[", reg_name,
                             "]] and try again")
                }
                return(Regime)
            }
            else
                stop("Intervention must be a list of named regimes, each specificed as ",
                     "list(intervention = f(A, L), g.star = g(A, L)). The intervention ",
                     "function f(A, L) must be a function of treatment and covariates that ",
                     "returns numeric desired treatment assignments (a*) with the same dimensions",
                     " as the observed treatment. The g.star function f(A, L) must be a function ",
                     "of treatment and covariates that returns numeric probabilities bounded in",
                     "[0, 1]. Amend Intervention[[", reg_name, "]] and try again.")
        })
    } else if (all(try(as.character(Intervention)) %in% c("0", "1"))) {
        Regime <- list()
        if ("1" %in% try(as.character(Intervention)))
            Regime[["A=1"]] <- list("intervention" = function(a, L) {rep_len(1, length(a))},
                                    "g.star" = function(a, L) {as.numeric(a == 1)})
        if ("0" %in% try(as.character(Intervention)))
            Regime[["A=0"]] <- list("intervention" = function(a, L) {rep_len(0, length(a))},
                                    "g.star" = function(a, L) {as.numeric(a == 0)})
    } else
        stop("Intervention must be \"0\", \"1\", c(\"0\", \"1\"), or a list of named regimes, ",
             "each specificed as list(intervention = f(A, L), g.star = g(A, L)). The ",
             "intervention function f(A, L) must be a function of treatment and covariates that ",
             "returns numeric desired treatment assignments (a*) with the same dimensions",
             " as the observed treatment. The g.star function f(A, L) must be a function ",
             "of treatment and covariates that returns numeric probabilities bounded in",
             "[0, 1].")
    names(Regime) <- names(Intervention)
    return(Regime)
}

getTargetEvent <- function(TargetEvent, UniqueEvents) {
    if (is.null(TargetEvent))
        TargetEvent <- UniqueEvents
    if (any(!is.vector(TargetEvent), !is.numeric(TargetEvent), is.list(TargetEvent),
            length(setdiff(TargetEvent, UniqueEvents)) > 0))
        stop("TargetEvent must be a numeric vector that is a subset of observed event types,",
             " DataTable[[EventType]], not including 0 (i.e. censoring)")
    return(TargetEvent)
}

checkTargetTime <- function(TargetTime, EventTime, TargetEvent, EventType) {
    if (any(!is.vector(TargetTime), !is.numeric(TargetTime), is.list(TargetTime),
            try(TargetTime <= 0)))
        stop("TargetTime must be a positive numeric vector.")
    tm.evnt <- data.table::data.table("EventTime" = EventTime,
                                      "EventType" = EventType)
    max.time <- tm.evnt[EventType > 0, ][, max(EventTime)]
    if (max(TargetTime) > max.time)
        stop("TargetTime must not target times after which all individuals are censored, ",
             max.time)
}

getCVFolds <- function(CVArg) {
    ## cross validation setup ----
    # stratifying cv so that folds are balanced for treatment assignment & outcomes
    # theory? but regressions may fail in practice with rare events otherwise ### make efficient CV representation ----
    cv.folds <- try(do.call(origami::make_folds, CVArg))
    if (inherits(cv.folds, "try-error"))
        stop("CVArg must be a list of arguments to be passed into do.call(origami::make_folds, ...)")
    return(cv.folds)
}

checkModel <- function(Model, UniqueEvent, Censored, HazEstBackend) {
    ## check that every event has a model
    if (!all(is.list(Model), length(Model) == length(UniqueEvent) + 1 + Censored))
        stop("Model must be a named list, one for each event type observed in the dataset")
    ## check model specifications are not obviously broken
    if (length(setdiff(as.character(UniqueEvent), names(Model))) > 0)
        stop("Model must be a named list, one for each event type observed in the dataset. ",
             "Model names must be `Trt` for the intervention, `0` for censoring, or the ",
             "corresponding numeric value for observed event type(s)")

    ## tag every model with an event attribute
    ## tag censoring as event = 0
    ## tag treatment model as treatment model


    # Model <- lapply(seq_along(Model), function(j) {
    #     haz.model <- Model[[j]]
    #     attr(haz.model, "j") <- as.numeric(gsub(".*(\\d+).*", "\\1", names(Model[j])))
    #     return(haz.model)
    # })
    warning("model checks not yet complete")
}

checkPropScoreBackend <- function(PropScoreBackend) {
    PropScoreBackendOK <- try(length(setdiff(PropScoreBackend, c("sl3"))) == 0)
    if (any(inherits(PropScoreBackendOK, "try-error"), !PropScoreBackendOK)) {
        stop("PropScoreBackend must now be `sl3`. Other options can be implemented in the future.")
    }
}

checkHazEstBackend <- function(HazEstBackend) {
    HazEstBackendOK <- try(length(setdiff(HazEstBackend, c("coxph"))) == 0)
    if (any(inherits(HazEstBackendOK, "try-error"), !HazEstBackendOK)) {
        stop("HazEstBackend must now be `coxph`. Other options can be implemented in the future.")
    }
}

getMaxUpdateIter <- function(MaxUpdateIter) {
    MaxUpdateIterOK <- try(all(is.numeric(MaxUpdateIter), length(MaxUpdateIter) == 1,
                               MaxUpdateIter > 0, !is.infinite(MaxUpdateIter)))
    if (any(inherits(MaxUpdateIterOK, "try-error"), !MaxUpdateIterOK)) {
        stop("MaxUpdateIter must a positive, finite whole number")
    }
    return(ceiling(MaxUpdateIter))
}

checkOneStepEps <- function(OneStepEps) {
    OneStepEpsOK <- try(all(is.numeric(OneStepEps), length(OneStepEps) == 1,
                            OneStepEps > 0, OneStepEps <= 1))
    if (any(inherits(OneStepEpsOK, "try-error"), !OneStepEpsOK)) {
        stop("OneStepEps must a positive number between (0, 1]")
    }
}

checkVerbose <- function(Verbose) {
    VerboseOK <- try(all(is.logical(Verbose), length(Verbose) == 1))
    if (any(inherits(VerboseOK, "try-error"), !VerboseOK)) {
        stop("Verbose must either be TRUE or FALSE")
    }
}

checkGComp <- function(GComp) {
    GCompOK <- try(all(is.logical(GComp), length(GComp) == 1))
    if (any(inherits(GCompOK, "try-error"), !GCompOK)) {
        stop("GComp must either be TRUE or FALSE")
    }
}

IntentToTreat <- list("A == 1" = list("intervention" = function(a, L) {rep_len(1, length(a))},
                                      "g.star" = function(a, L) {as.numeric(a == 1)}),
                      "A == 0" = list("intervention" = function(a, L) {rep_len(0, length(a))},
                                      "g.star" = function(a, L) {as.numeric(a == 0)}))

