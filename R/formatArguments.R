
#' formatArguments
#' @description formatArguments() checks and reformats inputs into a form that can be interpreted by doConcrete(). 
#'              makeITT() returns an Intervention list for a single, binary, point-treatment variable
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
# #'   \item{"LongTime"}{: numeric; Specifies monitoring times for longitudinal data structures}
#'   \item{"Baseline Covariates"}{: factor, character, or numeric; }
#' }
# #' @param DataStructure formula (not ready): e.g. Surv(time, type) ~ Intervention(trt) + ...
#' @param EventTime character: the column name of the observed event or censoring time
#' @param EventType character: the column name of the observed event type. (0 indicating censoring)
#' @param Treatment character: the column name of the observed treatment assignment
#' @param ID character (default: NULL): the column name of the observed subject id
# #' @param LongTime character (not used): the column name of the monitoring times for
#'                                       longitudinal data structures
#' @param Intervention list: a list of desired interventions on the treatment variable.
#'                           Each intervention must be a list containing two named functions: 
#'                             'intervention' = function(treatment vector, covariate data) and 
#'                             'gstar' = function(treatment vector, covariate data)
#'                           concrete::makeITT() can be used to specify an intent-to-treat analysis for a
#'                           binary intervention variable 
#' @param TargetTime numeric: vector of target times. If NULL, the last observed non-censoring event
#'                            time will be targeted.
#' @param TargetEvent numeric: vector of target events - some subset of unique EventTypes. If NULL, 
#'                             all non-censoring observed event types will be targeted.
# #' @param Target (not yet implemented) data.table / data.frame (?? x 2); a table containing all 
# #' combinations of target events (column 1) and target times (column 2).
#' @param CVArg list: arguments to be passed into do.call(origami::make_folds). If NULL, the default is 
#'                    list(n = nrow(DataTable), fold_fun = folds_vfold, cluster_ids = NULL, strata_ids = NULL)
#' @param Model list (default: NULL): named list of models, one for each failure or censoring event
#'                                    and one for the 'Treatment' variable. If Model = NULL, then  
#'                                    a template will be generated for the user to amend. 
#' @param PropScoreBackend character (default: "Superlearner"): currently must be either "sl3" or "Superlearner"
#' @param HazEstBackend character (default: "coxph"): currently must be "coxph"
#' @param MaxUpdateIter numeric (default: 100): the number of one-step update steps
#' @param OneStepEps numeric (default: 1): the one-step tmle step size
#' @param MinNuisance numeric (default: 0.05): the minimum value of the nuisance parameter denominator in the 
#' clever covariate
#' @param Verbose boolean
#' @param GComp boolean (default: TRUE): return g-computation formula plug-in estimates 
#' @param ReturnModels boolean (default: TRUE): return fitted models from the initial estimation stage
#' @param RenameCovs boolean (default: TRUE): whether or not to rename covariates
#' @param ... ...
#' 
#' @return a list of class "ConcreteArgs"
#' \itemize{
#'   \item{Data}{: data.table containing EventTime, EventType, Treatment, and potentially ID and baseline covariates. Has the following attributes}
#'   \itemize{
#'     \item{EventTime}{: the column name of the observed event or censoring time}
#'     \item{EventType}{: the column name of the observed event type. (0 indicating censoring)}
#'     \item{Treatment}{: the column name of the observed treatment assignment}
#'     \item{ID}{: the column name of the observed subject id}
# #'     \item{LongTime}{: the column name of the observed event type. (0 indicating censoring)}
#'     \item{RenameCovs}{: boolean whether or not covariates are renamed}
#'   }
#'   \item{TargetTime}{: numeric vector of target times to evaluate risk/survival}
#'   \item{TargetEvent}{: numeric vector of target events}
#'   \item{Regime}{: named list of desired regimes, each tagged with a 'g.star' attribute function}
#'     \itemize{
#'       \item{Regime\[\[i\]\]}{: a vector of desired treatment assignments}
#'       \item{attr(Regime\[\[i\]\], "g.star")}{: function of Treatment and Covariates, outputting a vector of desired treatment assignment probabilities}
#'     }
#'   \item{CVFolds}{: list of cross-validation fold assignments in the structure as output by origami::make_folds()}
#'   \item{Model}{: named list of model specifications, one for each unique 'EventType' and one for the 'Treatment' variable.}
#'   \item{PropScoreBackend}{: either "sl3" or "Superlearner"}
#'   \item{HazEstBackend}{: "coxph"}
#'   \item{MaxUpdateIter}{: the number of one-step update steps}
#'   \item{OneStepEps}{: list of cross-validation fold assignments in the structure as output by origami::make_folds()}
#'   \item{MinNuisance}{: numeric lower bound for the propensity score denominator in the efficient influence function}
#'   \item{Verbose}{: boolean to print additional information}
#'   \item{GComp}{: boolean to return g-computation formula plug-in estimates}
#'   \item{ReturnModels}{: boolean to return fitted models from the initial estimation stage}
#' }
#'
#' @importFrom stats model.matrix as.formula
#' @importFrom utils tail
#' @importFrom survival Surv coxph
#' @import origami data.table
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
#' 
#' # makeITT() creates a list of functions to specify intent-to-treat 
#' #   regimes for a binary, single, point treatment variable
#' intervention <- makeITT()
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
#' # if formatArguments(Model = NULL), a template will be returned for the user to modify
#' # examples of modifying/adding models for censoring and failure events
#' concrete.args[["Model"]][["0"]] <-
#'     list(Surv(time, status == 0) ~ trt:sex + age + bili)
#' concrete.args[["Model"]][["1"]] <- 
#'     list("mod1" = Surv(time, status == 1) ~ trt,
#'          "mod2" = Surv(time, status == 1) ~ .)
#' 
#' # 
#' # examples of modifying "Superlearner" treatment models
#' concrete.args[["Model"]][["trt"]] <- c("SL.glm", "SL.glmnet", "SL.bayesglm")
#' 
#' # examples of modifying "sl3" treatment models
#' library(sl3)
#' concrete.args[["Model"]][["trt"]] <-
#'     make_learner(Stack, Lrnr_hal9001$new(), Lrnr_glmnet$new(), Lrnr_glm$new())
#' 
#' @export formatArguments
#' @export makeITT

formatArguments <- function(DataTable, 
                            # DataStructure = NULL, 
                            EventTime, 
                            EventType, 
                            Treatment, 
                            Intervention, 
                            TargetTime = NULL, 
                            TargetEvent = NULL, 
                            # Target = NULL,
                            CVArg = NULL, 
                            Model = NULL, 
                            PropScoreBackend = "SuperLearner", 
                            HazEstBackend = "coxph",
                            MaxUpdateIter = 100, 
                            OneStepEps = 0.1, 
                            MinNuisance = 0.05,
                            Verbose = TRUE, 
                            GComp = TRUE, 
                            ReturnModels = TRUE, 
                            ConcreteArgs = NULL, 
                            ID = NULL, 
                            # LongTime = NULL, 
                            RenameCovs = TRUE, 
                            ...)
{
    ## Data Structure - incorporate prodlim::EventHistory.frame?
    if (!is.null(ConcreteArgs) | isTRUE(try(inherits(DataTable, "ConcreteArgs"), silent = TRUE))) {
        if (isTRUE(try(inherits(DataTable, "ConcreteArgs"), silent = TRUE)))
            ConcreteArgs <- DataTable
        if (!inherits(ConcreteArgs, "ConcreteArgs"))
            stop("ConcreteArgs must be of class 'ConcreteArgs', the output of ", 
                 "concrete::formatArguments()")
        DataTable <- ConcreteArgs[["Data"]]
        EventTime <- attr(DataTable, "EventTime")
        EventType <- attr(DataTable, "EventType")
        Treatment <- attr(DataTable, "Treatment")
        LongTime <- attr(DataTable, "LongTime")
        ID <- attr(DataTable, "ID")
        RenameCovs <- attr(DataTable, "RenameCovs")
        Intervention <- ConcreteArgs[["Regime"]]
        TargetEvent <- ConcreteArgs[["TargetEvent"]]
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
        ReturnModels <- ConcreteArgs[["ReturnModels"]]
    }
    LongTime <- NULL ## replace for future longitudinal compatibility ----
    DataTable <- formatDataTable(DT = DataTable, 
                                 EventTime = EventTime, 
                                 EventType = EventType, 
                                 Treatment = Treatment, 
                                 ID = ID, 
                                 LongTime = LongTime, 
                                 Verbose = Verbose, 
                                 RenameCovs = RenameCovs)
    
    TimeVal <- DataTable[[EventTime]]
    TypeVal <- DataTable[[EventType]]
    TrtVal <- DataTable[[Treatment]]
    CovNames <- setdiff(colnames(DataTable), c(EventTime, EventType, Treatment, ID, LongTime))
    CovDT <- subset(DataTable, select = CovNames)
    attr(CovDT, "CovNames") <- attr(DataTable, "CovNames")
    attr(CovDT, "RenameCovs") <- attr(DataTable, "RenameCovs")
    Censored <- 0 %in% TypeVal
    UniqueEvents <- setdiff(sort(unique(TypeVal)), 0)
    
    ## Interventions & Targets ----
    Regime <- getRegime(Intervention = Intervention, TrtVal = TrtVal, CovDT = CovDT)
    
    TargetEvent <- getTargetEvent(TargetEvent = TargetEvent, UniqueEvents = UniqueEvents)
    TargetTime <- getTargetTime(TargetTime = TargetTime, 
                                TimeVal = TimeVal, 
                                TargetEvent = TargetEvent,
                                TypeVal = TypeVal)
    
    ## Estimation Paramters ----
    if (!exists("CVSeed")) CVSeed <- sample(0:1e8, 1)
    CVFolds <- getCVFolds(CVArg = CVArg, DataTable = DataTable, CVSeed = CVSeed)
    checkPropScoreBackend(PropScoreBackend)
    checkHazEstBackend(HazEstBackend)
    Model <- getModel(Model = Model, 
                      UniqueEvents = UniqueEvents, 
                      Censored = Censored, 
                      HazEstBackend = HazEstBackend, 
                      PropScoreBackend = PropScoreBackend,
                      EventTime = EventTime, 
                      EventType = EventType, 
                      Treatment = Treatment, 
                      CovDT = CovDT)
    
    ## TMLE Update Parameters ----
    MaxUpdateIter <- getMaxUpdateIter(MaxUpdateIter)
    checkOneStepEps(OneStepEps)
    MinNuisance <- getMinNuisance(MinNuisance)
    
    ## Misc. Parameters
    checkVerbose(Verbose)
    checkGComp(GComp)
    
    ## return
    ConcreteArgs <- list(Data = DataTable,
                         TargetTime = TargetTime,
                         TargetEvent = TargetEvent,
                         Regime = Regime,
                         CVFolds = CVFolds,
                         Model = Model,
                         PropScoreBackend = PropScoreBackend,
                         HazEstBackend = HazEstBackend,
                         MaxUpdateIter = MaxUpdateIter,
                         OneStepEps = OneStepEps,
                         MinNuisance = MinNuisance,
                         Verbose = Verbose,
                         GComp = GComp, 
                         ReturnModels = ReturnModels)
    class(ConcreteArgs) <- union(class(ConcreteArgs), "ConcreteArgs")
    return(ConcreteArgs)
}

formatDataTable <- function(DT, EventTime, EventType, Treatment, ID, LongTime, Verbose, RenameCovs) {
    if (!inherits(DT, "data.table"))
        DT <- try(data.table::as.data.table(DT))
    if (!inherits(DT, "data.table"))
        stop("CovDataTable must be a data.table or coercible into a data.table.")
    if (any(is.infinite(unlist(DT)), anyNA(unlist(DT))))
        stop("CovDataTable contains infinite or missing values; regression models may break")
    
    checkEventTime(EventTime = EventTime, DataTable = DT)
    checkEventType(EventType = EventType, DataTable = DT)
    checkTreatment(Treatment = Treatment, DataTable = DT)
    ID <- getID(ID = ID, DataTable = DT)
    DT[[ID[["IDName"]]]] <- ID[["IDVal"]]
    ID <- ID[["IDName"]]
    LongTime <- NULL # LongTime <- getLongTime(LongTime = LongTime, DataTable = DT)
    
    CovNames <- setdiff(colnames(DT), c(EventTime, EventType, Treatment, ID, LongTime))
    if (identical(paste0(gsub("L\\d+", "", CovNames), collapse = ""), character(0)) | !RenameCovs) { 
        if (is.null(attr(DT, "CovNames")))
            attr(DT, "CovNames") <- data.table(ColName = colnames(DT),
                                               CovName = colnames(DT),
                                               CovVal = rep_len(".", ncol(DT)))
    } else {
        CovDT <- getCovDataTable(DataTable = DT, 
                                 EventTime = EventTime, 
                                 EventType = EventType, 
                                 Treatment = Treatment, 
                                 ID = ID, 
                                 LongTime = LongTime, 
                                 Verbose = Verbose)
        DT <- cbind(subset(DT, select = c(EventTime, EventType, Treatment, ID, LongTime)), CovDT)
        attr(DT, "CovNames") <- attr(CovDT, "CovNames")
    }
    attr(DT, "EventTime") <- EventTime
    attr(DT, "EventType") <- EventType
    attr(DT, "Treatment") <- Treatment
    attr(DT, "LongTime") <- LongTime
    attr(DT, "ID") <- ID
    attr(DT, "RenameCovs") <- RenameCovs
    return(DT)
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
        IDVal <- 1:nrow(DataTable)
        message("No ID column specified. DataTable row numbers will be used as subject IDs, ",
                "which will not be appropriate for longitudinal data structures.")
    }
    if (is.character(ID)) {
        IDVal <- try(DataTable[[ID]])
        if (inherits(IDVal, "try-error"))
            stop("No column named '", ID, "' was found in the supplied data. Check spelling ", 
                 "or input argument into DataTable")
    }
    if (any(is.list(IDVal), is.null(IDVal), is.nan(IDVal), is.na(IDVal)))
        stop("ID column must not include missing values")
    return(list(IDVal = IDVal, IDName = ID))
}

getLongTime <- function(LongTime, DataTable = NULL) {
    if (!is.null(LongTime)) stop("Longitudinal data structures not yet supported.")
    return(NULL)
}

getCovDataTable <- function(DataTable, EventTime, EventType, Treatment, ID, LongTime, Verbose) {
    `(Intercept)` <- NULL
    CovNames <- setdiff(colnames(DataTable), c(EventTime, EventType, Treatment, ID, LongTime))
    CovDT <- DataTable[, .SD ,.SDcols = CovNames]
    NonNumInd <- sapply(CovNames, function(CovName) {!inherits(CovDT[[CovName]], c("numeric", "integer"))})
    CovNames1Hot <- data.table()
    
    if (length(CovNames[!NonNumInd]) == 0) {
        CovDT1Hot <- data.table()
        l <- 0
    } else {
        CovDT1Hot <- CovDT[, .SD, .SDcols = CovNames[!NonNumInd]]
        CovNames1Hot <- CovNames1Hot[, list(ColName = paste0("L", 1:ncol(CovDT1Hot)),
                                            CovName = colnames(CovDT1Hot), 
                                            CovVal = rep_len(".", ncol(CovDT1Hot)))]
        setnames(CovDT1Hot, colnames(CovDT1Hot), paste0("L", 1:ncol(CovDT1Hot)))
        l <- ncol(CovDT1Hot)
        
        if (length(CovNames[NonNumInd]) == 0) {
            attr(CovDT1Hot, "CovNames") <- CovNames1Hot
            return(CovDT1Hot)
        } else 
            message("Categorical covariates detected: DataTable will be 1-hot encoded. New columns can ", 
                    "be linked with original columns through attr(.[['Data']], 'CovNames')")
    }
    
    for (CovName in CovNames[NonNumInd]) {
        Cov1Hot <- as.data.table(model.matrix(~., subset(CovDT, select = CovName)))[, `(Intercept)` := NULL]
        setnames(Cov1Hot, colnames(Cov1Hot), sub(CovName, "", colnames(Cov1Hot)))
        CovNames1Hot <- rbind(CovNames1Hot,
                              data.table(ColName = paste0("L", l + 1:ncol(Cov1Hot)),
                                         CovName = CovName,
                                         CovVal = colnames(Cov1Hot)))
        setnames(Cov1Hot, colnames(Cov1Hot), paste0("L", l + 1:ncol(Cov1Hot)))
        l <- l + ncol(Cov1Hot)
        CovDT1Hot <- cbind(CovDT1Hot, Cov1Hot)
    }
    
    # if (Verbose) try(superheat::superheat(cov(scale(model.matrix(~., CovDT1Hot)))))
    
    attr(CovDT1Hot, "CovNames") <- CovNames1Hot
    return(CovDT1Hot)
}

getRegime <- function(Intervention, TrtVal, CovDT) {
    if (is.list(Intervention)) {
        Regimes <- lapply(seq_along(Intervention), function(i) {
            Regime <- Intervention[[i]]
            if (is.null(names(Intervention))) {
                RegName <- paste0("intervention", i)
            } else if (names(Intervention)[i] == "") {
                RegName <- paste0("intervention", i)
            } else 
                RegName <- names(Intervention)[i]
            
            # regime ----
            if (is.numeric(unlist(Regime))) {
                RegimeVal <- Regime
            } else if (is.function(Regime[["intervention"]])) {
                RegimeVal <- try(do.call(Regime[["intervention"]], list(TrtVal, CovDT)))
                if (inherits(RegimeVal, "try-error"))
                    stop("Intervention must be a list of regimes specificed as list(intervention",
                         " = f(A, L), g.star = g(A, L)), and intervention function f(A, L) must ",
                         "be a function of treatment and covariates that returns numeric desired",
                         " treatment assignments (a*) with the same dimensions as the observed ",
                         "treatment. Amend Intervention[[", RegName, "]] and try again")
            }
            
            if (!all(length(unlist(RegimeVal)) == length(unlist(TrtVal)),
                     min(unlist(RegimeVal)) >= min(unlist(TrtVal)), 
                     max(unlist(RegimeVal)) <= max(unlist(TrtVal))))
                stop("If providing an intervention function, the function output must return a ",
                     "numeric vector of desired treatment assignments (a*) with the same ",
                     "dimensions as the observed treatment and with values within the observed ", 
                     "range. If providing intervention values, then the input must have the same ", 
                     "dimensions as the observed treatment withvalues within the observed range.",
                     "Amend Intervention[[", RegName, "]] and try again")
            
            if (mean(unlist(RegimeVal) == unlist(TrtVal)) < 0.05)
                warning("The intervention function f(A, L) specified in Intervention[[", RegName,
                        "]] is likely to cause near-positivity issues, resulting in inflated ",
                        "variance and perhaps instability. Recommend specifying an intervention",
                        "better supported in the observed data.")
            
            # g.star ----
            if (!is.null(attr(Regime, "g.star"))) {
                attr(RegimeVal, "g.star") <- attr(Regime, "g.star")
            } else if (is.list(Regime)) {
                attr(RegimeVal, "g.star") <- Regime[["g.star"]]
            } 
            if (is.null(attr(RegimeVal, "g.star"))) {
                attr(RegimeVal, "g.star") <- function(a) as.numeric(a == RegimeVal)
                message("No g.star function specified, defaulting to the indicator that observed",
                        "treatment equals the desired treatment assignment, 1(A = a*).")
            } 
            
            GStarOK <- try(do.call(attr(RegimeVal, "g.star"), list(TrtVal, CovDT)))
            if (inherits(GStarOK, "try-error") | !is.numeric(GStarOK)) {
                stop("Intervention must be a list of regimes specificed as list(intervention",
                     " = f(A, L), g.star = g(A, L)), and the g.star function f(A, L) must ",
                     "be a function of treatment and covariates that returns numeric ",
                     "probabilities bounded in [0, 1]. Amend Intervention[[", RegName,
                     "]] and try again")
            } else {
                if (min(unlist(GStarOK)) < 0 | max(unlist(GStarOK)) > 1)
                    stop("Intervention must be a list of regimes specificed as list(intervention",
                         " = f(A, L), g.star = g(A, L)), and the g.star function f(A, L) must ",
                         "be a function of treatment and covariates that returns numeric ",
                         "probabilities bounded in [0, 1]. Amend Intervention[[", RegName,
                         "]] and try again")
            }
            return(RegimeVal)
        }) 
    }
    else if (all(try(as.character(Intervention)) %in% c("0", "1"))) {
        Regimes <- list()
        if ("1" %in% try(as.character(Intervention))) {
            Regimes[["Treated"]] <- rep_len(1, length(TrtVal))
            attr(Regimes[["Treated"]], "g.star") <- function(a, L) {as.numeric(a == 1)}
        }
        if ("0" %in% try(as.character(Intervention))) {
            Regimes[["Control"]] <- rep_len(0, length(TrtVal))
            attr(Regimes[["Control"]], "g.star") <- function(a, L) {as.numeric(a == 0)}
        }
    } else
        stop("Intervention must be \"0\", \"1\", c(\"0\", \"1\"), or a list of named regimes, ",
             "each specificed as list(intervention = f(A, L), g.star = g(A, L)). The ",
             "intervention function f(A, L) must be a function of treatment and covariates that ",
             "returns numeric desired treatment assignments (a*) with the same dimensions",
             " as the observed treatment. The g.star function f(A, L) must be a function ",
             "of treatment and covariates that returns numeric probabilities bounded in",
             "[0, 1].")
    names(Regimes) <- names(Intervention)
    return(Regimes)
}

getTargetEvent <- function(TargetEvent, UniqueEvents) {
    if (is.null(TargetEvent))
        message("No TargetEvent specified; targeting all observed event types except for censoring")
    TargetEvent <- UniqueEvents
    if (any(!is.vector(TargetEvent), !is.numeric(TargetEvent), is.list(TargetEvent),
            length(setdiff(TargetEvent, UniqueEvents)) > 0))
        stop("TargetEvent must be a numeric vector that is a subset of observed event types,",
             " DataTable[[EventType]], not including 0 (i.e. censoring)")
    return(TargetEvent)
}

getTargetTime <- function(TargetTime, TimeVal, TargetEvent, TypeVal) {
    
    Times <- data.table::data.table("TimeVal" = TimeVal, "TypeVal" = TypeVal)
    MaxTime <- Times[TypeVal > 0, ][, max(TimeVal)]
    MinTime <- Times[TypeVal > 0, ][, list(TimeVal = min(TimeVal)), by = "TypeVal"]
    MinTimeEvents <- MinTime[["TypeVal"]]
    MinTime <- MinTime[["TimeVal"]]
    
    if (!is.null(TargetTime)) {
        if (any(!is.vector(TargetTime), !is.numeric(TargetTime), is.list(TargetTime), try(TargetTime <= 0)))
            stop("TargetTime must be a positive numeric vector.")
        if (max(TargetTime) > MaxTime)
            stop("TargetTime must not target times after which all individuals are Censored, ", MaxTime)
        if (any(min(TargetTime) < MinTime))
            message("TargetTime includes a time at which some events have not yet occurred - ", 
                    paste0(paste0("Event=", MinTimeEvents, ": ", MinTime), collapse = ", "))
    } else{
        TargetTime <- MaxTime
        message("No TargetTime provided; targeting the last observed event time by default, which may ", 
                "result in estimates with high variance if most subjects have been censored by that time")
    }
    return(TargetTime)
}

getCVFolds <- function(CVArg, DataTable, CVSeed = sample(0:1e8, 1)) {
    if (is.null(CVArg))
        CVArg <- list(n = nrow(DataTable), fold_fun = folds_vfold, cluster_ids = NULL, strata_ids = NULL)
    ## cross validation setup ----
    # stratifying cv so that folds are balanced for treatment assignment & outcomes
    # theory? but regressions may fail in practice with rare events otherwise ### make efficient CV representation ----
    set.seed(CVSeed)
    CVFolds <- try(do.call(origami::make_folds, CVArg))
    if (inherits(CVFolds, "try-error"))
        stop("CVArg must be a list of arguments to be passed into do.call(origami::make_folds, ...)")
    attr(CVFolds, "CVArg") <- CVArg
    attr(CVFolds, "CVSeed") <- CVSeed
    return(CVFolds)
}

getModel <- function(Model, UniqueEvents, Censored, PropScoreBackend, HazEstBackend, 
                     EventTime, EventType, Treatment, CovDT) {
    CovName <- NULL
    CovNames <- attr(CovDT, "CovNames")
    RenameCovs <- attr(CovDT, "RenameCovs")
    
    if (is.null(Model)) {
        message("Model input missing. An example template will be returned but should be amended to",  
                " suit your application. See examples in the concrete::formatArguments() documentation.")
        return(getModelTemplate(Treatment = Treatment, UniqueEvents = UniqueEvents, Censored = Censored, 
                                EventTime = EventTime, EventType = EventType))
    }
    ## check that treatment and every event (including censoring) has a model
    if (!all(is.list(Model), length(Model) == length(UniqueEvents) + 1 + Censored))
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
    if (!all(as.character(UniqueEvents) %in% names(Model)))
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
    
    CovNamesChanged <- FALSE
    for (FitVar in names(Model)) {
        if (grepl("\\d+", FitVar)) {
            if (HazEstBackend == "coxph") {
                if (is.null(Model[[FitVar]])) {
                    Model[[FitVar]] <- list("model1" = ~ .)
                    message("No model was provided for event ", FitVar, ", so by default ", 
                            "a main terms cox model with treatment and all covariates will be used.")
                }
                if (is.list(Model[[FitVar]])) {
                    if (is.null(names(Model[[FitVar]]))) {
                        names(Model[[FitVar]]) <- paste0("model", seq_along(Model[[FitVar]]))
                    } else if (any(names(Model[[FitVar]]) == "")) {
                        j <- which(names(Model[[FitVar]]) == "")
                        names(Model[[FitVar]])[j] <- paste0("model", j)
                    }
                } else {
                    Model[[FitVar]] <- list("model1" = Model[[FitVar]])
                }
                
                CoxLeft <- paste0("Surv(", EventTime, ", ", EventType, " == ", FitVar, ") ~ ")
                CoxLeftRegex <- paste0("^Surv\\(\\s*", EventTime, "\\s*,\\s*", EventType, 
                                       "\\s*==\\s*", FitVar, "\\s*\\)\\s*~\\s*")
                for (j in seq_along(Model[[FitVar]])) {
                    Formula <- as.character(Model[[FitVar]][j])
                    if (!grepl(CoxLeftRegex, Formula)) {
                        message("The left hand side of the Cox model formula for Model[[\"", 
                                FitVar, "\"]][[", j, "]] has been corrected to ", CoxLeft)
                    }
                    
                    CoxRight <- sub("^.*~", "", Formula)  
                    
                    # rename covariates ----
                    
                    if (!isTRUE(attr(Model[[FitVar]][[j]], "NameChecked")) & RenameCovs) {
                        for (covar in unique(CovNames[["CovName"]])) {
                            NewCol <- CovNames[CovName == covar, ][["ColName"]]
                            if (length(NewCol) > 1)
                                NewCol <- paste0(NewCol, collapse = "+")
                            CoxRight <- gsub(covar, paste0("(", NewCol, ")"), CoxRight)
                            CovNamesChanged <- TRUE
                        }
                    }
                    Model[[FitVar]][[j]] <- as.formula(paste0(CoxLeft, CoxRight))
                    attr(Model[[FitVar]][[j]], "NameChecked") <- TRUE
                }
            } else 
                stop("Models must be named for the treatment variable, or a numeric value ", 
                     "representing the failure or censoring event type")
        }
        
        
    }
    if (CovNamesChanged) {
        message("Cox model specifications have been renamed where necessary to reflect", 
                " changed covariate names. Model specifications in .[['Model']] can be ", 
                "checked against the covariate names in attr(.[['Data']], 'CovNames')")
    }
    # warning("model checks not yet complete")
    return(Model)
}

getModelTemplate <- function(Treatment, UniqueEvents, Censored, EventType, EventTime) {
    TrtModel <- list(c("SL.xgboost", "SL.glmnet"))
    Events <- utils::tail(c(0, UniqueEvents), length(UniqueEvents) + Censored)
    EventModels <- lapply(sort(Events), function(j) {
        EventModel <- list("model1" = as.formula(paste0("Surv(", EventTime, ", ", 
                                                        EventType, " == ", j, ") ~ ", Treatment)),
                           "model2" = as.formula(paste0("Surv(", EventTime, ", ", 
                                                        EventType, " == ", j, ") ~ .")))
        attr(EventModel[["model1"]], "NameChecked") <- TRUE
        attr(EventModel[["model2"]], "NameChecked") <- TRUE
        return(EventModel)
    })
    
    Model <- c(TrtModel, EventModels) ;
    names(Model) <- c(Treatment, sort(Events))
    return(Model)
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
getMinNuisance <- function(MinNuisance = 0.05) {
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

checkReturnModels <- function(ReturnModels) {
    ReturnModelsOK <- try(all(is.logical(ReturnModels), length(ReturnModels) == 1))
    if (any(inherits(ReturnModelsOK, "try-error"), !ReturnModelsOK)) {
        stop("ReturnModels must either be TRUE or FALSE")
    }
}

#' @describeIn formatArguments 
makeITT <- function() {
    ITT <- list("A=1" = list("intervention" = function(ObservedTreatment, Covariates) {
        IntervenedAssignment <- rep_len(1, length(ObservedTreatment))
        return(IntervenedAssignment)
    },
    "g.star" = function(Treatment, Covariates) {
        IntervenedProbability <- as.numeric(Treatment == 1)
        return(IntervenedProbability)
    }),
    "A=0" = list("intervention" = function(ObservedTreatment, Covariates) {
        IntervenedAssignment <- rep_len(0, length(ObservedTreatment))
        return(IntervenedAssignment)
    },
    "g.star" = function(Treatment, Covariates) {
        IntervenedProbability <- as.numeric(Treatment == 0)
        return(IntervenedProbability)
    }))
    return(ITT)
}

