
#' Title
#' @param ConcreteArgs list (default: NULL, not yet ready) : Use to recheck amended output from previous formatArguments() 
#'                                            calls. A non-NULL input will cause all other arguments to be ignored.
#' @param DataTable data.table (n x (d + (3:5)); data.table of the observed data, with rows n = 
#' the number of observations and d = the number of baseline covariates. DataTable must include 
#' the following columns:
#' \itemize{
#'   \item{"EventTime"}{: non-negative real numbers; the observed event or censoring time}
#'   \item{"EventType"}{: numeric; the observed event type, with 0 indicating censoring. There is 
#'   no separate column for indicating censoring}
#'   \item{"Treatment"}{: numeric; the observed treatment}
#' }
#' May include
#' \itemize{
#'   \item{"ID"}{: factor, character, or numeric; unique subject id. If ID column is missing, row 
#'   numbers will be used as ID. For longitudinal data, ID must be provided}
#'   \item{"LongTime"}{: numeric; Specifies monitoring times for longitudinal data structures}
#'   \item{"Baseline Covariates"}{: factor, character, or numeric; }
#' }
#' @param DataStructure formula (not ready): e.g. Surv(time, type) ~ Intervention(trt) + ...
#' @param EventTime character: the column name of the observed event or censoring time
#' @param EventType character: the column name of the observed event type. (0 indicating censoring)
#' @param Treatment character: the column name of the observed treatment assignment
#' @param ID character (default: NULL): the column name of the observed subject id
#' @param LongTime character (not used): the column name of the monitoring times for
#'                                       longitudinal data structures
#' @param Intervention list: a list of desired interventions on the treatment variable.
#'                           Each intervention must be a list containing two named functions: 
#'                             'intervention' = function(treatment vector, covariate data) and 
#'                             'gstar' = function(treatment vector, covariate data)
#'                           concrete:::ITT can be used to specify an intent-to-treat analysis for a
#'                           binary intervention variable 
#' @param TargetTime numeric: vector of target times
#' @param TargetEvent numeric: vector of target events - some subset of unique EventTypes. 
#' @param Target (not yet implemented) data.table / data.frame (?? x 2); a table containing all 
#' combinations of target events (column 1) and target times (column 2).
#' @param CVArg list: arguments to be passed into do.call(origami::make_folds). The default is 
#'                    list(n = nrow(DataTable), fold_fun = folds_vfold, cluster_ids = NULL, strata_ids = NULL)
#' @param Model list (default: NULL): named list of models, one for each failure or censoring event
#'                                    and one for the 'Treatment' variable. If Model = NULL, then  
#'                                    a template will be generated for the user to amend. 
#' @param PropScoreBackend character (default: "Superlearner"): currently must be either "sl3" or "Superlearner"
#' @param HazEstBackend character (default: "coxph"): currently must be "coxph"
#' @param MaxUpdateIter numeric: the number of one-step update steps
#' @param OneStepEps numeric: the one-step tmle step size
#' @param MinNuisance numeric: the minimum value of the nuisance parameter denominator in the 
#' clever covariate
#' @param Verbose boolean
#' @param GComp boolean
#' @param ... ...
#'
#' @return a list of class "concrete.args"
#' \itemize{
#'   \item{"Data"}{: data.table containing EventTime, EventType, Treatment, and baseline covariates}
#'   \item{"CovDataTable"}{: data.table containing 1-hot encoded baseline covariates}
#'   \item{"LongTime"}{: numeric vector of longitudinal monitoring times (not yet used)}
#'   \item{"ID"}{: vector of unique subject ids}
#'   \item{"Events"}{: numeric vector encoding unique failure event types}
#'   \item{"Censored"}{: indicator of the presence of right-censoring}
#'   \item{"TargetTime"}{: numeric vector of target times to evaluate risk/survival}
#'   \item{"TargetEvent"}{: numeric vector of target events}
#'   \item{"Regime"}{: named list of interventions, comprised of two functions}
#'     \itemize{
#'       \item{"intervention"}{: function of Treatment and Covariates, outputting a vector of desired treatment assignments}
#'       \item{"g.star"}{: function of Treatment and Covariates, outputting a vector of desired treatment assignment probabilities}
#'     }
#'   \item{"CVFolds"}{: list of cross-validation fold assignments in the structure as output by origami::make_folds()}
#'   \item{"Model"}{: list of cross-validation fold assignments in the structure as output by origami::make_folds()}
#'   \item{"PropScoreBackend"}{: list of cross-validation fold assignments in the structure as output by origami::make_folds()}
#'   \item{"HazEstBackend"}{: list of cross-validation fold assignments in the structure as output by origami::make_folds()}
#'   \item{"MaxUpdateIter"}{: list of cross-validation fold assignments in the structure as output by origami::make_folds()}
#'   \item{"OneStepEps"}{: list of cross-validation fold assignments in the structure as output by origami::make_folds()}
#'   \item{"MinNuisance"}{: numeric cutoff}
#'   \item{"Verbose"}{: boolean}
#'   \item{"GComp"}{: boolean, to return g-computation formula plug-in estimates or not}
#' }
#' May include
#' \itemize{
#'   \item{"ID"}{: factor, character, or numeric; unique subject id. If ID column is missing, row 
#'   numbers will be used as ID. For longitudinal data, ID must be provided}
#'   \item{"LongTime"}{: numeric; Specifies monitoring times for longitudinal data structures}
#'   \item{"Baseline Covariates"}{: factor, character, or numeric; }
#' }
#'
#' @importFrom stats model.matrix as.formula
#' @importFrom utils tail
#' @import origami
#'
#' @examples 
#' library(data.table)
#' library(concrete)
#' 
#' data <- as.data.table(survival::pbc)
#' data[, trt := sample(0:1, nrow(data), TRUE)]
#' cols <- c("id", "time", "status", "trt",
#'           "age", "albumin", "sex", "bili")
#' data <- data[, .SD, .SDcols = cols]
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
#'                                  Model = model)
#' 
#' # if formatArguments(Model = NULL), a model template will be returned for the user to amend.
#' # examples of editing models for censoring and failure events
#' concrete.args[["Model"]][["0"]] <- list("model1" = Surv(time, status == 0) ~ trt:sex + age + bili)
#' concrete.args[["Model"]][["1"]] <- list(Surv(time, status == 1) ~ trt, 
#'                                        Surv(time, status == 1) ~ .)
#' 
#' # examples of editing models for binary treatment, using PropScoreBackend = "Superlearner"
#' concrete.args[["Model"]][["trt"]] <- c("SL.glm", "SL.glmnet", "SL.bayesglm")
#' 
#' # examples of editing models for binary treatment, using PropScoreBackend = "sl3"
#' library(sl3)
#' concrete.args[["Model"]][["trt"]] <- make_learner(Stack, Lrnr_hal9001$new(), 
#'                                                  Lrnr_glmnet$new(), Lrnr_glm$new())
#' 
#' 
#' @export formatArguments

formatArguments <- function(DataTable, DataStructure = NULL, EventTime, EventType, Treatment, ID = NULL, 
                            LongTime = NULL, Intervention, TargetTime, TargetEvent = NULL, Target = NULL,
                            CVArg = list(n = nrow(DataTable), fold_fun = folds_vfold,
                                         cluster_ids = NULL, strata_ids = NULL),
                            Model = NULL, PropScoreBackend = "SuperLearner", HazEstBackend = "coxph",
                            MaxUpdateIter = 100, OneStepEps = 0.1, MinNuisance = 0.05,
                            Verbose = TRUE, GComp = TRUE, ConcreteArgs = NULL, ...)
{
    ## Data Structure ----
    # incorporate prodlim::EventHistory.frame?
    if (!is.null(ConcreteArgs)) {
        if (!inherits(ConcreteArgs, "concrete.args"))
            stop("ConcreteArgs must be of class 'concrete.args', the output of ", 
                 "concrete::formatArguments()")
        DataTable <- ConcreteArgs[["Data"]]
        EventTime <- attr(ConcreteArgs[["Data"]], "EventTime")
        EventType <- attr(ConcreteArgs[["Data"]], "EventType")
        Treatment <- attr(ConcreteArgs[["Data"]], "Treatment")
        LongTime <- attr(ConcreteArgs[["Data"]], "LongTime")
        ID <- attr(ConcreteArgs[["Data"]], "ID")
        Intervention <- ConcreteArgs[["Regime"]]
        TargetTime <- ConcreteArgs[["TargetTime"]]
        CVArg <- attr(ConcreteArgs[["CVFolds"]], "CVArg")
        CVSeed <- attr(ConcreteArgs[["CVFolds"]], "CVSeed")
        PropScoreBackend <- ConcreteArgs[["PropScoreBackend"]]
        HazEstBackend <- ConcreteArgs[["HazEstBackend"]]
        Model <- ConcreteArgs[["Model"]]
        MaxUpdateIter <- ConcreteArgs[["MaxUpdateIter"]]
        OneStepEps <- ConcreteArgs[["OneStepEps"]]
        MinNuisanc <- ConcreteArgs[["MinNuisanc"]]
        Verbose <- ConcreteArgs[["Verbose"]]
        GComp <- ConcreteArgs[["GComp"]]
    }
    data.tbl <- formatDataTable(x = DataTable, EventTime = EventTime, EventType = EventType, ID = ID, 
                                Treatment = Treatment, LongTime = LongTime, Verbose = Verbose)
    
    event.time <- data.tbl[[EventTime]]
    event.type <- data.tbl[[EventType]]
    treatment <- data.tbl[[Treatment]]
    id <- data.tbl[[ID]]
    # long.time <- getLongTime(x = LongTime, DataTable = data.tbl)
    cov.data.table <- getCovDataTable(DataTable = data.tbl, EventTime = EventTime,
                                      EventType = EventType, Treatment = Treatment, ID = ID,
                                      LongTime = LongTime, Verbose = Verbose)
    
    censored <- 0 %in% event.type
    event.type.unique <- setdiff(sort(unique(event.type)), 0)
    
    ## Interventions & Targets ----
    regime <- getRegime(Intervention = Intervention, Treatment = treatment,
                        CovDataTable = cov.data.table)
    
    target.event <- getTargetEvent(TargetEvent = TargetEvent, UniqueEvents = event.type.unique)
    checkTargetTime(TargetTime = TargetTime, EventTime = event.time, TargetEvent = target.event,
                    EventType = event.type)
    
    ## Estimation Paramters ----
    if (!exists("CVSeed")) CVSeed <- sample(0:1e8, 1)
    cv.folds <- getCVFolds(CVArg = CVArg, CVSeed = CVSeed)
    checkPropScoreBackend(PropScoreBackend)
    checkHazEstBackend(HazEstBackend)
    model <- getModel(Model = Model, UniqueEvent = event.type.unique, Censored = censored, 
                      HazEstBackend = HazEstBackend, PropScoreBackend = PropScoreBackend,
                      EventTime = EventTime, EventType = EventType, Treatment = Treatment, 
                      CovDataTable = cov.data.table)
    
    ## TMLE Update Parameters ----
    max.update.iter <- getMaxUpdateIter(MaxUpdateIter)
    checkOneStepEps(OneStepEps)
    min.nuisance <- getMinNuisance(MinNuisance)
    
    ## Misc. Parameters
    checkVerbose(Verbose)
    checkGComp(GComp)
    
    ## return
    concrete.args <- list(Data = data.tbl,
                          CovDataTable = cov.data.table,
                          LongTime = long.time,
                          ID = id,
                          Events = event.type.unique,
                          Censored = censored,
                          TargetTime = TargetTime,
                          TargetEvent = target.event,
                          Regime = regime,
                          CVFolds = cv.folds,
                          Model = model,
                          PropScoreBackend = PropScoreBackend,
                          HazEstBackend = HazEstBackend,
                          MaxUpdateIter = max.update.iter,
                          OneStepEps = OneStepEps,
                          MinNuisance = min.nuisance,
                          Verbose = Verbose,
                          GComp = GComp)
    class(concrete.args) <- union(class(concrete.args), "concrete.args")
    return(concrete.args)
}

formatDataTable <- function(x, EventTime, EventType, Treatment, ID, LongTime, Verbose) {
    if (!inherits(x, "data.table"))
        x <- try(data.table::as.data.table(x))
    if (!inherits(x, "data.table"))
        stop("CovDataTable must be a data.table or coercible into a data.table.")
    if (any(is.infinite(unlist(x)), anyNA(unlist(x))))
        stop("CovDataTable contains infinite or missing values; regression models may break")
    
    checkEventTime(EventTime = EventTime, DataTable = x)
    checkEventType(EventType = EventType, DataTable = x)
    checkTreatment(Treatment = Treatment, DataTable = x)
    id <- getID(ID = ID, DataTable = x)
    ID <- id[["ID"]]
    x[[ID]] <- id[["id"]]
    # long.time <- getLongTime(LongTime = LongTime, DataTable = x)
    
    cov.data.table <- getCovDataTable(DataTable = x, EventTime = EventTime, EventType = EventType, 
                                      Treatment = Treatment, ID = ID, LongTime = LongTime, 
                                      Verbose = Verbose)
    
    x <- cbind(subset(x, select = c(EventTime, EventType, Treatment, ID, LongTime)), cov.data.table)
    
    attr(x, "cov.names") <- attr(cov.data.table, "cov.names")
    attr(x, "EventTime") <- EventTime
    attr(x, "EventType") <- EventType
    attr(x, "Treatment") <- Treatment
    attr(x, "LongTime") <- LongTime
    attr(x, "ID") <- ID
    return(x)
}

checkEventTime <- function(EventTime, DataTable = NULL) {
    if (is.character(EventTime)) {
        tmp <- try(DataTable[[EventTime]])
        if (inherits(tmp, "try-error"))
            stop("No column named '", EventTime, "' was found in the supplied data. Check the ", 
                 "`EventTime=` and the `DataTable=` argument inputs")
        attr(tmp, "var.name") <- EventTime
        EventTime <- tmp
        if (any(!is.numeric(EventTime), try(EventTime <= 0), is.infinite(EventTime), 
                inherits(try(EventTime <= 0), "try-error"), is.list(EventTime)))
            stop("EventTime must be a numeric vector with positive, finite values")
    } else 
        stop("`EventTime` must be the name of the column containing the observed event times.")
}

checkEventType <- function(EventType, DataTable = NULL) {
    if (is.character(EventType)) {
        tmp <- try(DataTable[[EventType]])
        if (inherits(tmp, "try-error"))
            stop("No column named '", EventType, "' was found in the supplied data. Check spelling ", 
                 "or input argument into DataTable")
        attr(tmp, "var.name") <- EventType
        EventType <- tmp
        if (any(!is.numeric(EventType), try(EventType < 0), 
                inherits(try(EventType < 0), "try-error"), is.list(EventType)))
            stop("EventType must be a numeric vector with non-negative values (0 indicating censoring)")
    } else 
        stop("`EventType` must be the name of the column containing the observed event types (", 
             "with 0 indicating censoring).")
}

checkTreatment <- function(Treatment, DataTable = NULL) {
    if (is.character(Treatment)) {
        tmp <- try(DataTable[[Treatment]])
        if (inherits(tmp, "try-error"))
            stop("No column named '", Treatment, "' was found in the supplied data. Check spelling ", 
                 "or input argument into DataTable")
        attr(tmp, "var.name") <- Treatment
        Treatment <- tmp
        if (any(!is.numeric(Treatment), is.nan(Treatment), is.infinite(Treatment), is.list(Treatment)))
            stop("Treatment must be a numeric vector with finite values")
    } else 
        stop("`Treatment` must be the name of the column containing the intervention variable.")
    return(Treatment)
}

getID <- function(ID, DataTable = NULL) {
    if (is.null(ID)) {
        ID <- "ID"
        id <- 1:nrow(DataTable)
        message("No ID column specified. DataTable row numbers will be used as subject IDs, ",
                "which will not be appropriate for longitudinal data structures.")
    }
    if (is.character(ID)) {
        id <- try(DataTable[[ID]])
        if (inherits(id, "try-error"))
            stop("No column named '", ID, "' was found in the supplied data. Check spelling ", 
                 "or input argument into DataTable")
        attr(id, "var.name") <- ID
    }
    if (any(is.list(id), is.null(id), is.nan(id), is.na(id)))
        stop("ID column must not include missing values")
    return(list(id = id, ID = ID))
}

getLongTime <- function(LongTime, DataTable = NULL) {
    if (!is.null(LongTime)) stop("Longitudinal data structures not yet supported.")
    return(NULL)
}

getCovDataTable <- function(DataTable, EventTime, EventType, Treatment, ID, LongTime, Verbose) {
    "(Intercept)" <- NULL
    cov.names <- setdiff(colnames(DataTable), c(EventTime, EventType, Treatment, ID, LongTime))
    cov.dt <- DataTable[, .SD ,.SDcols = cov.names]
    non.num.ind <- sapply(cov.names, function(cov.name) {!inherits(cov.dt[[cov.name]], c("numeric", "integer"))})
    cov.names.1hot <- data.table()
    
    if (length(cov.names[non.num.ind]) == 0) return(cov.dt)
    
    if (length(cov.names[!non.num.ind]) == 0) {
        cov.dt.1hot <- data.table()
        l <- 0
    } else {
        cov.dt.1hot <- cov.dt[, .SD, .SDcols = cov.names[!non.num.ind]]
        cov.names.1hot <- cov.names.1hot[, list(col.name = paste0("L", 1:ncol(cov.dt.1hot)),
                                                cov.name = colnames(cov.dt.1hot), 
                                                cov.val = colnames(cov.dt.1hot))]
        setnames(cov.dt.1hot, colnames(cov.dt.1hot), paste0("L", 1:ncol(cov.dt.1hot)))
        l <- ncol(cov.dt.1hot)
    }
    
    for (cov.name in cov.names[non.num.ind]) {
        cov.1hot <- as.data.table(model.matrix(~., subset(cov.dt, select = cov.name)))[, `(Intercept)` := NULL]
        cov.names.1hot <- rbind(cov.names.1hot,
                                data.table(col.name = paste0("L", l + 1:ncol(cov.1hot)),
                                           cov.name = cov.name,
                                           cov.val = colnames(cov.1hot)))
        setnames(cov.1hot, colnames(cov.1hot), paste0("L", l + 1:ncol(cov.1hot)))
        l <- l + ncol(cov.1hot)
        
        cov.dt.1hot <- cbind(cov.dt.1hot, cov.1hot)
    }
    
    attr(cov.dt.1hot, "cov.names") <- cov.names.1hot
    
    if (Verbose) {
        # try(superheat::superheat(cov(scale(model.matrix(~., cov.dt.1hot)))))
    }
    return(cov.dt.1hot)
}

getRegime <- function(Intervention, Treatment, CovDataTable) {
    if (is.list(Intervention)) {
        Regime <- lapply(seq_along(Intervention), function(i) {
            regime <- Intervention[[i]]
            if (is.null(names(Intervention))) {
                reg_name <- paste0("intervention", i)
            } else if (names(Intervention)[i] == "") {
                reg_name <- paste0("intervention", i)
            } else 
                reg_name <- names(Intervention)[i]
            
            # regime ----
            if (is.numeric(unlist(regime))) {
                Regime <- regime
            } else if (is.function(regime[["intervention"]])) {
                Regime <- try(do.call(regime[["intervention"]], list(Treatment, CovDataTable)))
                if (inherits(Regime, "try-error"))
                    stop("Intervention must be a list of regimes specificed as list(intervention",
                         " = f(A, L), g.star = g(A, L)), and intervention function f(A, L) must ",
                         "be a function of treatment and covariates that returns numeric desired",
                         " treatment assignments (a*) with the same dimensions as the observed ",
                         "treatment. Amend Intervention[[", reg_name, "]] and try again")
            }
            
            if (!all(length(unlist(Regime)) == length(unlist(Treatment)),
                     min(unlist(Regime)) >= min(unlist(Treatment)), 
                     max(unlist(Regime)) <= max(unlist(Treatment))))
                stop("If providing an intervention function, the function output must return a ",
                     "numeric vector of desired treatment assignments (a*) with the same ",
                     "dimensions as the observed treatment and with values within the observed ", 
                     "range. If providing intervention values, then the input must have the same ", 
                     "dimensions as the observed treatment withvalues within the observed range.",
                     "Amend Intervention[[", reg_name, "]] and try again")
            
            if (mean(unlist(Regime) == unlist(Treatment)) < 0.05)
                warning("The intervention function f(A, L) specified in Intervention[[", reg_name,
                        "]] is likely to cause near-positivity issues, resulting in inflated ",
                        "variance and perhaps instability. Recommend specifying an intervention",
                        "better supported in the observed data.")
            
            # g.star ----
            if (!is.null(attr(regime, "g.star"))) {
                attr(Regime, "g.star") <- attr(regime, "g.star")
            } else if (is.list(regime)) {
                attr(Regime, "g.star") <- regime[["g.star"]]
            } 
            if (is.null(attr(Regime, "g.star"))) {
                attr(Regime, "g.star") <- function(a) as.numeric(a == Regime)
                message("No g.star function specified, defaulting to the indicator that observed",
                        "treatment equals the desired treatment assignment, 1(A = a*).")
            } 
            
            GStarOK <- try(do.call(attr(Regime, "g.star"), list(Treatment, CovDataTable)))
            if (inherits(GStarOK, "try-error") | !is.numeric(GStarOK)) {
                stop("Intervention must be a list of regimes specificed as list(intervention",
                     " = f(A, L), g.star = g(A, L)), and the g.star function f(A, L) must ",
                     "be a function of treatment and covariates that returns numeric ",
                     "probabilities bounded in [0, 1]. Amend Intervention[[", reg_name,
                     "]] and try again")
            } else {
                if (min(unlist(GStarOK)) < 0 | max(unlist(GStarOK)) > 1)
                    stop("Intervention must be a list of regimes specificed as list(intervention",
                         " = f(A, L), g.star = g(A, L)), and the g.star function f(A, L) must ",
                         "be a function of treatment and covariates that returns numeric ",
                         "probabilities bounded in [0, 1]. Amend Intervention[[", reg_name,
                         "]] and try again")
            }
            return(Regime)
        }) 
    }
    else if (all(try(as.character(Intervention)) %in% c("0", "1"))) {
        Regime <- list()
        if ("1" %in% try(as.character(Intervention)))
            Regime[["Treated"]] <- list("intervention" = function(a, L) {rep_len(1, length(a))},
                                        "g.star" = function(a, L) {as.numeric(a == 1)})
        if ("0" %in% try(as.character(Intervention)))
            Regime[["Control"]] <- list("intervention" = function(a, L) {rep_len(0, length(a))},
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
    class(Regime) <- union(class(Regime), "concrete.regime")
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
    if (any(!is.vector(TargetTime), !is.numeric(TargetTime), is.list(TargetTime), try(TargetTime <= 0)))
        stop("TargetTime must be a positive numeric vector.")
    tm.evnt <- data.table::data.table("EventTime" = EventTime, "EventType" = EventType)
    max.time <- tm.evnt[EventType > 0, ][, max(EventTime)]
    min.time <- tm.evnt[EventType > 0, ][, list(EventTime = min(EventTime)), by = "EventType"]
    min.time.event <- min.time[["EventType"]]
    min.time <- min.time[["EventTime"]]
    
    if (max(TargetTime) > max.time)
        stop("TargetTime must not target times after which all individuals are censored, ", max.time)
    
    if (any(min(TargetTime) < min.time))
        warning("TargetTime is targeting times before any events of type(s): ", 
                paste(min.time.event[min(TargetTime) < min.time], collapse = ", "))
}

getCVFolds <- function(CVArg, CVSeed = sample(0:1e8, 1)) {
    ## cross validation setup ----
    # stratifying cv so that folds are balanced for treatment assignment & outcomes
    # theory? but regressions may fail in practice with rare events otherwise ### make efficient CV representation ----
    set.seed(CVSeed)
    cv.folds <- try(do.call(origami::make_folds, CVArg))
    if (inherits(cv.folds, "try-error"))
        stop("CVArg must be a list of arguments to be passed into do.call(origami::make_folds, ...)")
    attr(cv.folds, "CVArg") <- CVArg
    attr(cv.folds, "CVSeed") <- CVSeed
    return(cv.folds)
}

getModel <- function(Model, UniqueEvent, Censored, PropScoreBackend, HazEstBackend, 
                     EventTime, EventType, Treatment, CovDataTable) {
    if (is.null(Model)) {
        message("Model input missing. An example template will be returned but should be amended to",  
                " suit your application. See examples in the concrete::formatArguments() documentation.")
        return(getModelTemplate(Treatment = Treatment, UniqueEvent = UniqueEvent, Censored = Censored, 
                                EventTime = EventTime, EventType = EventType))
    }
    ## check that treatment and every event (including censoring) has a model
    if (!all(is.list(Model), length(Model) == length(UniqueEvent) + 1 + Censored))
        stop("Model must be a named list, one for each event type observed in the dataset, ", 
             "including censoring, and one for the treatment variable.")
    
    ## check that model specifications are named correctly
    if (!(Treatment %in% names(Model)))
        stop("A named list must be provided, containg model specifications for the treatment ", 
             "variable. This list must be named with the treatment variable name (i.e.", 
             "formatArguments(Treatment = ...). Run formatArguments(Models=NULL) to get ", 
             "an example of the required formatting.")
    if (Censored & !("0" %in% names(Model)))
        stop("Data includes an EventType = 0, indicating the presence of right-censoring, so a ", 
             "list named `0` containing model specifications for censoring must be provided. Run ", 
             "formatArguments(Models=NULL) to get an example of the required formatting.")
    if (!all(as.character(UniqueEvent) %in% names(Model)))
        stop("For every unique value of EventType, a list named with the corresponding numeric ", 
             "value must be provided, containing the model specifications for the time-to-event. ", 
             "Run formatArguments(..., Models=NULL) to get an example of the required formatting.")
    
    ## check trt model fits with backend
    if (PropScoreBackend == "sl3") {
        if (!inherits(Model[[Treatment]], "R6") | !inherits(Model[[Treatment]], "Lrnr_base"))
            stop("For PropScoreBackend = `sl3`, the model(s) for Treatment must be R6 objects ", 
                 "produced by sl3::make_learner() or related functions. See examples in the ", 
                 "formatArguments() documentation or the sl3 chapter of the tlverse handbook (", 
                 "https://tlverse.org/tlverse-handbook/sl3.html)")
    } else if (PropScoreBackend == "SuperLearner") {
        message("Superlearner model specifications are not checked here. The input must be a valid", 
                " argument into the `sl.lib = ` argument of Superlearner::Superlearner()")
    } else 
        stop("PropScoreBackend must be either `sl3` or `SuperLearner`.")
    
    for (i in 1:length(Model)) {
        if (names(Model)[i] == Treatment) {
            message("no model specification check for treatment.")
        } else if (grepl("\\d+", names(Model)[i])) {
            if (is.list(Model[[i]])) {
                if (is.null(names(Model[[i]]))) {
                    names(Model[[i]]) <- paste0("model", 1:length(Model[[i]]))
                } else if (any(names(Model[[i]]) == "")) {
                    j <- which(names(Model[[i]]) == "")
                    names(Model[[i]])[j] <- paste0("model", j)
                }
            } else 
                Model[[i]] <- list("model1" = Model[[i]])
            
            if (HazEstBackend == "coxph") {
                cox.left <- paste0("Surv(", EventTime, ", ", EventType, 
                                   " == ", names(Model)[i], ") ~ ")
                for (j in seq_along(Model[[i]])) {
                    form.j <- as.character(Model[[i]][j])
                    if (!grepl(cox.left, form.j)) {
                        cox.right.j <- regmatches(form.j, regexpr("^.*~", form.j), invert = TRUE)
                        cox.right.j <- tail(unlist(cox.right.j), 1)
                        
                        # rename covariates ----
                        cov.names <- attr(CovDataTable, "cov.names")
                        if (!is.null(cov.names)) {
                            for (covar in unique(cov.names[["cov.name"]])) {
                                new.col <- cov.names[cov.name == covar, ][["col.name"]]
                                if (length(new.col) > 1)
                                    new.col <- paste0(new.col, collapse = "+")
                                new.col <- paste0("(", new.col, ")")
                                old.col.ind <- regexpr(covar, cox.right.j)
                                regmatches(cox.right.j, old.col.ind, invert = FALSE) <- new.col
                            }
                        }
                    }
                    Model[[i]][[j]] <- as.formula(paste0(cox.left, cox.right.j))
                }
            } else 
                stop("Models must be named for the treatment variable, or the numeric value ", 
                     "representing the failure or censoring event type")
        }
    }
    
    ## tag every model/model list as an target event, censoring, or treatment
    ## check model specifications match data set. Translate model specifications with factor variable 
    ##  names into model specs using 1-hot encoded variable names.
    
    warning("model checks not yet complete")
    
    return(Model)
}

getModelTemplate <- function(Treatment, UniqueEvent, Censored, EventType, EventTime) {
    trt.model <- "SL.glmnet"
    events <- utils::tail(c(0, UniqueEvent), length(UniqueEvent) + Censored)
    event.model <- lapply(sort(events), function(j) {
        j.model <- list("model1" = as.formula(paste0("Surv(", EventTime, ", ", 
                                                     EventType, " == ", j, ") ~ .")))
        return(j.model)
    })
    
    model <- c(trt.model, event.model) ;
    names(model) <- c(Treatment, sort(events))
    return(model)
}

checkPropScoreBackend <- function(PropScoreBackend) {
    PropScoreBackendOK <- try(length(setdiff(PropScoreBackend, c("sl3", "SuperLearner"))) == 0)
    if (any(inherits(PropScoreBackendOK, "try-error"), !PropScoreBackendOK)) {
        stop("Currently PropScoreBackend can only be `sl3` or 'SuperLearner'.",
             " Other options may be implemented in the future.")
    }
}

checkHazEstBackend <- function(HazEstBackend) {
    HazEstBackendOK <- try(length(setdiff(HazEstBackend, c("coxph"))) == 0)
    if (any(inherits(HazEstBackendOK, "try-error"), !HazEstBackendOK)) {
        stop("Currently HazEstBackend can only be `coxph`.",
             " Other options may be implemented in the future.")
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
getMinNuisance <- function(MinNuisance) {
    MinNuisanceOK <- try(all(is.numeric(MinNuisance), length(MinNuisance) == 1,
                             MinNuisance > 0, MinNuisance <= 1))
    if (any(inherits(MinNuisanceOK, "try-error"), !MinNuisanceOK)) {
        stop("MinNuisance must a positive number between (0, 1]")
    }
    return(MinNuisance)
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

ITT <- list("A == 1" = list("intervention" = function(trt, covars) {rep_len(1, length(trt))},
                            "g.star" = function(trt, covars) {as.numeric(trt == 1)}),
            "A == 0" = list("intervention" = function(trt, covars) {rep_len(0, length(trt))},
                            "g.star" = function(trt, covars) {as.numeric(trt == 0)}))


