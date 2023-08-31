
#' formatArguments
#' @description formatArguments() checks and reformats inputs into a form that can be interpreted by doConcrete().
#'              makeITT() returns an Intervention list for a single, binary, point-treatment variable
#' @param ConcreteArgs list (default: NULL, not yet ready) : Use to recheck amended output from previous formatArguments()
#'                                            calls. A non-NULL input will cause all other arguments to be ignored.
#' @param DataTable data.table (n x (d + (3:5)); data.table of the observed data, with rows n =
#' the number of observations and d = the number of baseline covariates. DataTable must include
#' the following columns:
#' \itemize{
#'   \item{"EventTime"}{: numeric; real numbers > 0, the observed event or censoring time}
#'   \item{"EventType"}{: numeric; the observed event type, censoring events indicated by integers <= 0}
#'   \item{"Treatment"}{: numeric; the observed treatment value. Binary treatments must be coded as 0, 1}
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
#' @param MaxUpdateIter numeric (default: 500): the number of one-step update steps
#' @param OneStepEps numeric (default: 1): the one-step tmle step size
#' @param MinNuisance numeric (default: 5/log(n)/sqrt(n)): value between (0, 1) for truncating the g-related denominator of the clever covariate
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
#'   \item{MaxUpdateIter}{: the number of one-step update steps}
#'   \item{OneStepEps}{: list of cross-validation fold assignments in the structure as output by origami::make_folds()}
#'   \item{MinNuisance}{: numeric lower bound for the propensity score denominator in the efficient influence function}
#'   \item{Verbose}{: boolean to print additional information}
#'   \item{GComp}{: boolean to return g-computation formula plug-in estimates}
#'   \item{ReturnModels}{: boolean to return fitted models from the initial estimation stage}
#' }
#'
#' @importFrom stats model.matrix as.formula quantile
#' @importFrom utils tail head capture.output
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
#' # library(sl3)
#' # concrete.args[["Model"]][["trt"]] <-
#' #     make_learner(Stack, Lrnr_hal9001$new(), Lrnr_glmnet$new(), Lrnr_glm$new())
#'
#' @export formatArguments
#' @export makeITT

formatArguments <- function(DataTable,
                            # DataStructure = NULL,
                            EventTime,
                            EventType,
                            Treatment,
                            ID = NULL,
                            # LongTime = NULL,
                            # DataStructure = NULL,
                            TargetTime = NULL,
                            TargetEvent = NULL,
                            Intervention,
                            # Target = NULL,
                            CVArg = NULL,
                            Model = NULL,
                            MaxUpdateIter = 500,
                            OneStepEps = 0.1,
                            MinNuisance = 5/sqrt(nrow(DataTable))/log(nrow(DataTable)),
                            Verbose = TRUE,
                            GComp = TRUE,
                            ReturnModels = TRUE,
                            ConcreteArgs = NULL,
                            # LongTime = NULL,
                            RenameCovs = TRUE,
                            ...)
{
    ## Data Structure - incorporate prodlim::EventHistory.frame?
    if (!is.null(ConcreteArgs) | isTRUE(try(inherits(DataTable, "ConcreteArgs"), silent = TRUE))) {
        if (isTRUE(try(inherits(DataTable, "ConcreteArgs"), silent = TRUE)))
            ConcreteArgs <- DataTable
        if (!inherits(ConcreteArgs, "ConcreteArgs"))
            stop("ConcreteArgs must be of class 'ConcreteArgs', the output of formatArguments()")
    } else {
        ConcreteArgs <- makeConcreteArgs(DataTable, EventTime, EventType, Treatment, Intervention,
                                         TargetTime, TargetEvent, CVArg, Model, 
                                         MaxUpdateIter, OneStepEps, MinNuisance, 
                                         Verbose, GComp, ReturnModels, ID, RenameCovs)
    }
    
    with(ConcreteArgs, {
        # Miscellaneous Args ----
        checkBoolean(ArgList = list("Verbose" = Verbose, 
                                    "GComp" = GComp, 
                                    "ReturnModels" = ReturnModels, 
                                    "RenameCovs" = RenameCovs), 
                     Envir = ConcreteArgs)
        
        
        # Data Spec ----
        DataTable <- formatDataTable(DT = DataTable,
                                     EventTime = EventTime,
                                     EventType = EventType,
                                     Treatment = Treatment,
                                     ID = ID,
                                     LongTime = NULL,
                                     Verbose = Verbose,
                                     RenameCovs = RenameCovs) 
        
        
        # Interventions & Targets ----
        Regime <- getRegime(Intervention = Intervention, Data = DataTable)
        TargetEvent <- getTargetEvent(TargetEvent = TargetEvent, Data = DataTable)
        TargetTime <- getTargetTime(TargetTime = TargetTime, TargetEvent = TargetEvent, Data = DataTable)
        
        
        # Estimation Spec ----
        if (is.null(ConcreteArgs[["CVFolds"]]))
            CVFolds <- getCVFolds(CVArg = CVArg, Data = DataTable, CVSeed = sample(0:1e8, 1))
        Model <- getModel(Model = Model,
                          Data = DataTable,
                          Verbose = Verbose)
        
        # TMLE Update Spec ----
        MaxUpdateIter <- getMaxUpdateIter(MaxUpdateIter)
        OneStepEps <- checkOneStepEps(OneStepEps)
        MinNuisance <- getMinNuisance(MinNuisance)
    })
    
    print.ConcreteArgs(x = ConcreteArgs, Verbose = ConcreteArgs[["Verbose"]])
    return(ConcreteArgs)
}


makeConcreteArgs <- function(DataTable, EventTime, EventType, Treatment, Intervention,
                             TargetTime, TargetEvent, CVArg, Model, 
                             MaxUpdateIter, OneStepEps, MinNuisance, 
                             Verbose, GComp, ReturnModels, ID, RenameCovs) {
    ConcreteArgs <- new.env()
    with(ConcreteArgs, {
        DataTable <- DataTable
        EventTime <- EventTime
        EventType <- EventType
        Treatment <- Treatment
        Intervention <- Intervention
        TargetTime <- TargetTime
        TargetEvent <- TargetEvent
        CVArg <- CVArg
        Model <- Model
        MaxUpdateIter <- MaxUpdateIter
        OneStepEps <- OneStepEps
        MinNuisance <- MinNuisance
        Verbose <- Verbose
        GComp <- GComp
        ReturnModels <- ReturnModels
        ID <- ID
        RenameCovs <- RenameCovs
    })
    class(ConcreteArgs) <- union("ConcreteArgs", class(ConcreteArgs))
    return(ConcreteArgs)
}


checkBoolean <- function(ArgList, Envir) {
    lapply(seq_along(ArgList), function(i) {
        ArgOK <- try(all(is.logical(ArgList[[i]]), length(ArgList[[i]]) == 1, 
                         !is.na(ArgList[[i]])))
        if (any(inherits(ArgOK, "try-error"), !ArgOK, is.null(ArgOK))) {
            cat("Argument '", names(ArgList)[i], "' must be either TRUE or FALSE, ", 
                "so has been set to FALSE by default\n", sep = "")
            assign(names(ArgList)[i], FALSE, envir = Envir)
        } else {
            assign(names(ArgList)[i], ArgList[[i]], envir = Envir)
        }
    })
    invisible(ArgList)
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
    checkTreatment(Treatment = Treatment, EventType = EventType, DataTable = DT)
    IDArgs <- getID(ID = ID, DataTable = DT)
    ID <- IDArgs[["IDName"]]
    DT[[ID]] <- IDArgs[["IDVal"]]
    nEff <- length(unique(DT[[ID]]))
    LongTime <- NULL # LongTime <- getLongTime(LongTime = LongTime, DataTable = DT)
    
    OrigCovDT <- attr(DT, "OrigCovDT")
    SpecialCols <- c(ID, EventTime, EventType, Treatment, LongTime)
    CovNames <- setdiff(colnames(DT), SpecialCols)
    
    if (RenameCovs) {
        if (is.null(attr(DT, "CovNames"))) {
            CovDT <- getCovDataTable(DataTable = DT,
                                     EventTime = EventTime,
                                     EventType = EventType,
                                     Treatment = Treatment,
                                     ID = ID,
                                     LongTime = LongTime,
                                     Verbose = Verbose)
            DT <- cbind(DT[, .SD, .SDcols = SpecialCols], CovDT)
            attr(DT, "CovNames") <- attr(CovDT, "CovNames")
            OrigCovDT <- attr(CovDT, "OrigCovDT")
        }
    } else {
        if (is.null(attr(DT, "CovNames")))
            attr(DT, "CovNames") <- data.table(ColName = CovNames,
                                               CovName = CovNames,
                                               CovVal = rep_len(".", length(CovNames)))
    }
    DT <- structure(DT,
                    EventTime = EventTime,
                    EventType = EventType,
                    Treatment = Treatment,
                    LongTime = LongTime,
                    ID = ID, 
                    nEff = nEff,
                    RenameCovs = RenameCovs, 
                    OrigCovDT = OrigCovDT)
    setcolorder(DT, SpecialCols)
    return(DT)
}

checkEventTime <- function(EventTime, DataTable = NULL) {
    if (is.character(EventTime)) {
        tmp <- try(DataTable[[EventTime]])
        if (inherits(tmp, "try-error") | is.null(tmp))
            stop("No column named '", EventTime, "' was found in the supplied DataTable")
        if (any(!is.numeric(tmp), try(tmp <= 0), is.infinite(tmp),
                inherits(try(tmp <= 0), "try-error"), is.list(tmp)))
            stop("The 'EventTime' column must be finite, positive values without missingness.")
    } else
        stop("`EventTime` must be the name of the column containing the observed event times.")
    invisible(NULL)
}

checkEventType <- function(EventType, DataTable = NULL) {
    if (is.character(EventType)) {
        tmp <- try(DataTable[[EventType]])
        if (inherits(tmp, "try-error") | is.null(tmp))
            stop("No column named '", EventType, "' was found in the supplied DataTable.")
        if (any(!is.numeric(tmp), try(tmp < 0), inherits(try(tmp < 0), "try-error"), is.list(tmp)))
            stop("The 'EventType' column must be finite, non-negative values without missingness, ",
                 "with 0 indicating censoring")
    } else
        stop("`EventType` must be the name of the column containing the observed event types (",
             "with 0 indicating the onset of right censoring).")
    invisible(NULL)
}

checkTreatment <- function(Treatment, EventType, DataTable = NULL) {
    if (is.character(Treatment)) {
        tmp <- try(DataTable[, .SD , .SDcols = Treatment])
        if (inherits(tmp, "try-error") | is.null(tmp))
            stop("Column", ifelse(length(Treatment) > 1, "s", ""), " '", paste0(Treatment, collapse = ", "), 
                 "' ", ifelse(length(Treatment) > 1, "were", "was"), " not found in the supplied data.",  
                 "Check Treatment and DataTable argument inputs \n", attr(tmp, "condition"))
        if (any(tmp %in% DataTable[[EventType]]))
            stop("The name of the Treatment column(s) must be different from all EventType values. ", 
                 "Rename the treatment column(s) or recode EventType values to different integers.")
        attr(tmp, "var.name") <- Treatment
        apply(tmp, 2, function(trt) {
            if (any(!is.numeric(trt), is.nan(trt), is.infinite(trt), is.list(trt)))
                stop("Treatment must be a numeric vector with finite values. Encode binary variables ",
                     "as 0 or 1 and encode multinomial (factor) variables as integers.")
        })
    } else
        stop("The `Treatment` argument input must be a character vector containing the column name(s) corresponding to the intervention variable(s).")
    invisible(NULL)
}

getID <- function(ID, DataTable = NULL) {
    if (is.null(ID)) {
        ID <- "ID"
        IDVal <- 1:nrow(DataTable)
        # cat("No ID column specified. DataTable row numbers will be used as subject IDs, ",
        #     "which will not be appropriate for longitudinal data structures.\n", sep = "")
    } else if (is.character(ID)) {
        IDVal <- try(DataTable[[ID]])
        if (inherits(IDVal, "try-error") | is.null(IDVal))
            stop("No column named '", ID, "' was found in the supplied data. Check spelling ",
                 "or input argument into DataTable")
    }
    if (any(is.list(IDVal), is.null(IDVal), is.nan(IDVal), is.na(IDVal)))
        stop("ID column must not include missing values")
    return(list(IDVal = IDVal, IDName = ID))
}

# getLongTime <- function(LongTime, DataTable = NULL) {
#     if (!is.null(LongTime)) stop("Longitudinal data structures not yet supported.")
#     return(NULL)
# }

getCovDataTable <- function(DataTable, EventTime, EventType, Treatment, ID, LongTime, Verbose) {
    `(Intercept)` <- NULL
    CovNames <- setdiff(colnames(DataTable), c(EventTime, EventType, Treatment, ID, LongTime))
    CovDT <- DataTable[, .SD , .SDcols = CovNames]
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
        } 
        # else
        # cat("Categorical covariates detected: DataTable will be 1-hot encoded. New columns can ",
        #     "be linked with original columns through attr(.[['Data']], 'CovNames')\n", sep = "")
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
    attr(CovDT1Hot, "OrigCovDT") <- CovDT
    return(CovDT1Hot)
}

getRegime <- function(Intervention, Data) {
    TrtVal <- Data[, .SD, .SDcols = attr(Data, "Treatment")]
    CovDT <- attr(Data, "OrigCovDT")
    if (all(Intervention %in% c(0, 1), length(Intervention) %in% c(1, 2)))
        Intervention <- makeITT()[c("1", "0") %in% try(as.character(Intervention))]
    if (is.list(Intervention)) {
        Regimes <- lapply(seq_along(Intervention), function(i) {
            if (inherits(Intervention, "concreteIntervention") | 
                all(length(Intervention) %in% c(1, 2), sapply(Intervention, is.function))) {
                Intervention <- list(Intervention)
            }
            Regime <- Intervention[[i]]
            
            if (is.null(names(Intervention))) {
                RegName <- paste0("intervention", i)
            } else if (names(Intervention)[i] == "") {
                RegName <- paste0("intervention", i)
            } else
                RegName <- names(Intervention)[i]
            
            # intervention fn ----
            # Regime input as function, dataframe, vector of individual trt values, single value
            if (is.list(Regime)) {
                if (is.function(Regime[[1]]) & length(Regime) <= 2) {
                    # check intervention function
                    if (is.null(names(Regime)[1]))  
                        names(Regime)[1] <- "intervention"
                    RegimeVal <- try(do.call(Regime[["intervention"]], list(TrtVal, CovDT)))
                    if (inherits(RegimeVal, "try-error") & ncol(TrtVal) == 1)
                        RegimeVal <- try(do.call(Regime[["intervention"]], list(unlist(TrtVal), CovDT)))
                    if (any(inherits(RegimeVal, "try-error"), is.null(RegimeVal), 
                            dim(RegimeVal) != dim(TrtVal), !is.numeric(unlist(RegimeVal))))
                        stop("Intervention specification has failed to produce numeric RegimeVals ", 
                             "in a data.table with dimensions matching TrtVal") 
                    # g.star function checked in the last section of getRegime()
                    if (!inherits(RegimeVal, "data.table")) {
                        RegimeVal <- data.table::as.data.table(RegimeVal)
                        data.table::setnames(RegimeVal, new = colnames(TrtVal))
                    }
                } else {
                    stop("List inputs into the Intervention argument must be a list of length ",  
                         "<= 2, containing a treatment assignment function and an optional", 
                         "propensity score function. See concrete::makeITT() for an example. ")
                } 
            } else if (all(inherits(Regime, c("data.frame", "data.table", "matrix")), 
                           all(dim(Regime) == dim(TrtVal)))) {
                if (!all(apply(Regime, 2, typeof) == apply(TrtVal, 2, typeof))) {
                    stop("Regime value types mismatched with Treatment")
                } else {
                    RegimeVal <- Regime
                }
            } else if (isTRUE(length(Regime) == nrow(TrtVal) & is.numeric(Regime))) {
                RegimeVal <- data.table::copy(TrtVal)
                for (i in 1:ncol(RegimeVal)) {
                    RegimeVal[, i] <- Regime
                }
            } else if (isTRUE(length(Regime) == 1) & is.numeric(Regime)) {
                RegimeVal <- data.table::copy(TrtVal)
                for (i in 1:ncol(RegimeVal)) {
                    RegimeVal[, i] <- rep_len(Regime, nrow(RegimeVal))
                }
            } else
                stop("Something has gone wrong with Intervention specification. ", 
                     "Check 'Intervention=' argument input or debug concrete:::getRegime()")
            
            if (!all(dim(RegimeVal) == dim(TrtVal),
                     min(unlist(RegimeVal)) >= min(unlist(TrtVal)),
                     max(unlist(RegimeVal)) <= max(unlist(TrtVal)), na.rm = TRUE))
                stop("If providing an intervention function, the function output must return a ",
                     "numeric vector of desired treatment assignments (a*) with the same ",
                     "dimensions as the observed treatment and with values within the observed ",
                     "range. If providing intervention values, then the input must have the same ",
                     "dimensions as the observed treatment with values within the observed range.",
                     "Amend Intervention[[", RegName, "]] and try again")
            
            # check g.star ----
            if (!is.null(attr(Regime, "g.star"))) {
                attr(RegimeVal, "g.star") <- attr(Regime, "g.star")
            } else if (is.list(Regime)) {
                attr(RegimeVal, "g.star") <- Regime[["g.star"]]
            }
            if (is.null(attr(RegimeVal, "g.star"))) {
                attr(RegimeVal, "g.star") <- function(Treatment, Covariates, PropScore, Intervened) {
                    Probability <- data.table(1 * sapply(1:nrow(Treatment), function(i) 
                        all(Treatment[i, ] == Intervened[i, ])))
                    return(Probability)
                }
                cat("No g.star function specified, defaulting to the indicator that observed",
                    "treatment equals the desired treatment assignment, 1(A = a*).\n", sep = "")
            }
            PropScoreDummy <- data.table::copy(TrtVal)
            TrtNames <- colnames(TrtVal)
            setnames(PropScoreDummy, TrtNames)
            PropScoreDummy[, (TrtNames) := lapply(.SD, function(ps) rep_len(0.5, nrow(TrtVal))), 
                           .SDcols = TrtNames]
            GStarOK <- try(do.call(attr(RegimeVal, "g.star"), 
                                   list(TrtVal, CovDT, PropScoreDummy, RegimeVal)))
            
            if (any(inherits(GStarOK, "try-error"), 
                    inherits(try(as.numeric(unlist(GStarOK))), "try-error"), 
                    is.null(GStarOK))) {
                stop("Intervention must be a list of regimes specificed as list(intervention",
                     " = f(A, L), g.star = g(A, L)), and the g.star function f(A, L) must ",
                     "be a function of treatment and covariates that returns numeric ",
                     "probabilities bounded in [0, 1]. Amend Intervention[[", RegName,
                     "]] and try again")
            } else {
                if (min(as.numeric(unlist(GStarOK))) < 0 | max(as.numeric(unlist(GStarOK))) > 1)
                    stop("Intervention must be a list of regimes specificed as list(intervention",
                         " = f(A, L), g.star = g(A, L)), and the g.star function f(A, L) must ",
                         "be a function of treatment and covariates that returns numeric ",
                         "probabilities bounded in [0, 1]. Amend Intervention[[", RegName,
                         "]] and try again")
            }
            class(RegimeVal) <- union("concreteIntervention", class(RegimeVal))
            return(RegimeVal)
        })
        if (is.null(names(Intervention))) {
            names(Regimes) <- paste0("Regime", seq_along(Intervention))
        } else
            names(Regimes) <- names(Intervention)
    } else
        stop("Intervention must be integers corresponding to desired static treatment level(s), ",  
             "or a list of named regimes, each specificed as list(intervention = f(A, L), ",
             "g.star = g(A, L)). The intervention function f(A, L) must be a function of treatment ",
             "and covariates that returns numeric desired treatment assignments (a*) with the same ",
             "dimensions as the observed treatment. The g.star function f(A, L) must be a function ",
             "of treatment and covariates that returns numeric probabilities bounded in",
             "[0, 1]. See concrete::makeITT() for an example.")
    return(Regimes) 
}

getTargetEvent <- function(TargetEvent, Data) {
    UniqueEvents <- sort(unique(Data[[attr(Data, "EventType")]]))
    if (is.null(TargetEvent))
        TargetEvent <- UniqueEvents[UniqueEvents > 0]
    TargetEvent <- UniqueEvents[UniqueEvents != 0]
    if (any(!is.vector(TargetEvent), !is.numeric(TargetEvent), is.list(TargetEvent),
            length(setdiff(TargetEvent, UniqueEvents)) > 0))
        stop("TargetEvent must be a subset of the observed event types,",
             " DataTable[[\"", attr(Data, "EventType"), "\"]]): ", 
             paste0(UniqueEvents, collapse = ", "))
    return(TargetEvent)
}

getTargetTime <- function(TargetTime, TargetEvent, Data) {
    TypeVal <- TimeVal <- NULL
    Times <- data.table::data.table("TimeVal" = Data[[attr(Data, "EventTime")]], 
                                    "TypeVal" = Data[[attr(Data, "EventType")]])
    MaxTime <- Times[TypeVal %in% TargetEvent, ][, max(TimeVal)]
    MinTime <- Times[TypeVal %in% TargetEvent, ][, list(TimeVal = min(TimeVal)), by = "TypeVal"]
    MinTimeEvents <- MinTime[["TypeVal"]]
    MinTime <- MinTime[["TimeVal"]]
    
    if (!is.null(TargetTime)) {
        if (any(!is.vector(TargetTime), !is.numeric(TargetTime), is.list(TargetTime), try(TargetTime <= 0)))
            stop("TargetTime must be a positive numeric vector.")
        if (max(TargetTime) > MaxTime)
            stop("TargetTime must not target times after which all individuals are Censored, ", MaxTime)
        if (any(min(TargetTime) < MinTime))
            cat("TargetTime includes a time at which some events have not yet occurred - ",
                paste0(paste0("Event=", MinTimeEvents, ": ", MinTime), collapse = ", "), "\n", sep = "")
    } else{
        TargetTime <- MaxTime
        cat("No TargetTime provided; targeting the last observed event time by default, which may ",
            "result in estimates with high variance if most subjects have been censored by that time\n", sep = "")
    }
    return(TargetTime)
}

getCVFolds <- function(CVArg, Data, CVSeed = sample(0:1e8, 1)) {
    ## cross validation setup ----
    # stratified by event type to avoid regressions failing with rare events 
    # theory? but regressions may fail in practice with rare events otherwise 
    ### nice to do: make efficient CV representation
    nEff <- attr(Data, "nEff")
    V <- (nEff <= 30)*(nEff - 20) + (nEff <= 500)*10 + (nEff <= 5e3)*5 + (nEff <= 1e4)*2 + 3
    
    if (is.null(CVArg)) {
        CVArg <- list(n = nrow(Data), V = V, fold_fun = origami::folds_vfold, 
                      cluster_ids = Data[[attr(Data, "ID")]], 
                      strata_ids = Data[[attr(Data, "EventType")]])
    } else {
        if (is.list(CVArg)) {
            if (is.null(CVArg[["n"]]))
                CVArg[["n"]] <- nrow(Data)
            if (is.null(CVArg[["V"]]))
                CVArg[["V"]] <- V
            if (is.null(CVArg[["fold_fun"]]))
                CVArg[["fold_fun"]] <- origami::folds_vfold
            if (is.null(CVArg[["cluster_ids"]]))
                CVArg[["cluster_ids"]] <- Data[[attr(Data, "ID")]]
            if (is.null(CVArg[["strata_ids"]]))
                CVArg[["strata_ids"]] <- Data[[attr(Data, "EventType")]]
        } else {
            warning("CVArg input is not correctly formatted; default CVFolds have been generated.")
        }
    }
    set.seed(CVSeed)
    CVFolds <- try(do.call(origami::make_folds, CVArg))
    if (inherits(CVFolds, "try-error") | is.null(CVFolds))
        stop("CVArg must be a list of arguments to be passed into do.call(origami::make_folds, ...)")
    attr(CVFolds, "CVArg") <- CVArg
    attr(CVFolds, "CVSeed") <- CVSeed
    return(CVFolds)
}

getModel <- function(Model, Data, Verbose) {
    CovName <- NULL
    Treatment <- attr(Data, "Treatment")
    EventTime <- attr(Data, "EventTime")
    EventType <- attr(Data, "EventType")
    TypeVal <- Data[[EventType]]
    UniqueEvents <- sort(unique(TypeVal))
    CovNames <- attr(Data, "CovNames")
    RenameCovs <- attr(Data, "RenameCovs")
    
    Model <- makeModelList(Treatment = Treatment,
                           EventTime = EventTime, 
                           EventType = EventType, 
                           UniqueEvents = UniqueEvents, 
                           Model = Model, 
                           Verbose = Verbose)
    
    ## check hazard models
    CovNamesChanged <- FALSE
    for (FitVar in names(Model)) {
        if (!(FitVar %in% c(Treatment, UniqueEvents))) {
            cat("The Model[['", FitVar,"']] specification will be ignored. Check that model ",
                "specifications are named correspondingly to the treatment variable, or the ", 
                "numeric value representing a censoring or event type\n")
        } else {
            if (FitVar %in% UniqueEvents) {
                if (any(is.null(attr(Model[[FitVar]], "Backend")), 
                        attr(Model[[FitVar]], "Backend") != "coxph"))
                    attr(Model[[FitVar]], "Backend") <- "coxph"
                JBackend <- attr(Model[[FitVar]], "Backend")
                
                if (is.list(Model[[FitVar]])) {
                    if (is.null(names(Model[[FitVar]]))) {
                        names(Model[[FitVar]]) <- paste0("model", seq_along(Model[[FitVar]]))
                    } else if (any(names(Model[[FitVar]]) == "")) {
                        i <- which(names(Model[[FitVar]]) == "")
                        names(Model[[FitVar]])[i] <- paste0("model", i)
                    }
                } else {
                    Model[[FitVar]] <- list("model1" = Model[[FitVar]])
                }
                
                CoxLeft <- paste0("Surv(", EventTime, ", ", EventType, " == ", FitVar, ") ~ ")
                CoxLeftRegex <- paste0("^Surv\\(\\s*", EventTime, "\\s*,\\s*", EventType,
                                       "\\s*==\\s*", FitVar, "\\s*\\)\\s*~\\s*")
                for (j in seq_along(Model[[FitVar]])) {
                    # coxnet 
                    if (inherits(Model[[FitVar]][j], "SL.coxnet")) {
                        
                    } else if (inherits(Model[[FitVar]][j], "SL.cox")) {
                        # cox
                        
                    }
                    
                    Formula <- as.character(Model[[FitVar]][j])
                    CoxRight <- paste0(" ", sub("^.*~", "", Formula), " ")
                    
                    # rename covariates ----
                    
                    if (!isTRUE(attr(Model[[FitVar]][[j]], "NameChecked")) & RenameCovs) {
                        for (covar in unique(CovNames[["CovName"]])) {
                            OldColRegex <- paste0("([\\W\\D]|\\s){1}", covar, "([\\W\\D]|\\s){1}")
                            NewCol <- CovNames[CovName == covar, ][["ColName"]]
                            if (length(NewCol) > 1)
                                NewCol <- paste0(NewCol, collapse = "+")
                            CoxRight <- gsub(OldColRegex, paste0("\\1(", NewCol, ")\\2"), CoxRight)
                            CovNamesChanged <- TRUE
                        }
                    }
                    Model[[FitVar]][[j]] <- as.formula(paste0(CoxLeft, CoxRight))
                    attr(Model[[FitVar]][[j]], "NameChecked") <- TRUE
                }
                attr(Model[[FitVar]], "Backend") <- JBackend
            }
        }
    }
    attr(Model, "CovNamesChanged") <- CovNamesChanged
    return(Model)
}

makeModelList <- function(Treatment, EventTime, EventType, UniqueEvents, Model, Verbose) {
    if (is.null(Model))  
        Model <- list()
    
    # Propensity Score Estimators
    for (Trt in Treatment) {
        if (isTRUE(all(inherits(Model[[Trt]], "R6"), inherits(Model[[Trt]], "Lrnr_base")))) {
            attr(Model[[Trt]], "Backend") <- "sl3"
        } else if (isTRUE(inherits(Model[[Trt]], "SuperLearner"))) {
            attr(Model[[Trt]], "Backend") <- "SuperLearner"
        } else {
            SLSpecOK <- as.logical(all(sapply(Model[[Trt]], is.character)) * 
                                       (length(sapply(Model[[Trt]], is.character)) > 0))
            if (isTRUE(SLSpecOK)) {
                SLLrnrs <- suppressMessages(utils::capture.output(SuperLearner::listWrappers()))
                TrtLrnrs <- unlist(Model[[Trt]], recursive = TRUE)
                NonDefaultLrnrs <- TrtLrnrs[which(!(TrtLrnrs %in% SLLrnrs))]
                if (length(NonDefaultLrnrs > 0) & Verbose) {
                    cat("These candidate learners are not in the base SuperLearner package: ",
                        paste0(NonDefaultLrnrs, collapse = ","), "\n")
                }
            } else {
                if (!is.null(Model[[Trt]]))
                    cat("Model specification for \"", Trt, "\" is not recognized and will be ", 
                        "replaced with the default SuperLearner libraries.\n", sep = "")
                Model[[Trt]] <- c("SL.xgboost", "SL.glmnet")
            }
            attr(Model[[Trt]], "Backend") <- "SuperLearner"
        }
    }
    
    # Censoring and Event Hazard Estimators
    for (j in UniqueEvents) {
        if (is.null(Model[[as.character(j)]])) {
            Trt <- paste0(Treatment, collapse = " + ")
            HazModel <- list("TrtOnly" = as.formula(paste0("Surv(", EventTime, ", ",
                                                           EventType, " == ", j, ") ~ ", Trt)),
                             "MainTerms" = as.formula(paste0("Surv(", EventTime, ", ",
                                                             EventType, " == ", j, ") ~ .")))
            attr(HazModel[[1]], "NameChecked") <- TRUE
            attr(HazModel[[2]], "NameChecked") <- TRUE
            Model[[as.character(j)]] <- HazModel
            attr(Model[[as.character(j)]], "Backend") <- "CoxPH"
        }
        attr(Model[[as.character(j)]], "j") <- j
    }
    class(Model) <- union("ModelList", class(Model))
    return(Model)
}

getMaxUpdateIter <- function(MaxUpdateIter) {
    MaxUpdateIterOK <- try(all(is.numeric(MaxUpdateIter), length(MaxUpdateIter) == 1,
                               MaxUpdateIter > 0, !is.infinite(MaxUpdateIter)))
    if (any(inherits(MaxUpdateIterOK, "try-error"), is.null(MaxUpdateIterOK), !MaxUpdateIterOK)) {
        cat("MaxUpdateIter must a positive, finite whole number, so has been set to 100\n")
        MaxUpdateIter <- 100
    }
    return(ceiling(MaxUpdateIter))
}

checkOneStepEps <- function(OneStepEps) {
    OneStepEpsOK <- try(all(is.numeric(OneStepEps), length(OneStepEps) == 1,
                            OneStepEps > 0, OneStepEps <= 1))
    if (any(inherits(OneStepEpsOK, "try-error"), !OneStepEpsOK, is.null(OneStepEpsOK))) {
        cat("OneStepEps must a positive number between (0, 1], so has been set to 0.5\n")
        OneStepEps <- 0.5
    }
    return(OneStepEps)
}

getMinNuisance <- function(MinNuisance = 0.05) {
    MinNuisanceOK <- try(all(is.numeric(MinNuisance), length(MinNuisance) == 1,
                             MinNuisance > 0, MinNuisance <= 1))
    if (any(inherits(MinNuisanceOK, "try-error"), !MinNuisanceOK, is.null(MinNuisanceOK))) {
        cat("MinNuisance must a positive number between (0, 1], so has been set to 0.05\n")
        MinNuisance <- 0.05
    }
    return(MinNuisance)
}

#' @describeIn formatArguments makeITT ...
makeITT <- function(...) {
    Intrv <- list(...)
    if (length(Intrv) == 0) {
        ITT <- list("A=1" = list("intervention" = function(ObservedTrt, Covariates, PropScore) {
            TrtNames <- colnames(ObservedTrt)
            Intervened <- data.table::copy(ObservedTrt)
            Intervened[, (TrtNames) := lapply(.SD, function(a) 1), .SDcols = TrtNames]
            return(Intervened)
        },
        "g.star" = function(Treatment, Covariates, PropScore, Intervened) {
            Probability <- data.table(1 * sapply(1:nrow(Treatment), function(i) 
                all(Treatment[i, ] == Intervened[i, ])))
            return(Probability)
        }),
        "A=0" = list("intervention" = function(ObservedTrt, Covariates) {
            TrtNames <- colnames(ObservedTrt)
            Intervened <- data.table::copy(ObservedTrt)
            Intervened[, (TrtNames) := lapply(.SD, function(a) 0), .SDcols = TrtNames]
            return(Intervened)
        },
        "g.star" = function(Treatment, Covariates, PropScore, Intervened) {
            Probability <- data.table(1 * sapply(1:nrow(Treatment), function(i) 
                all(Treatment[i, ] == Intervened[i, ])))
            return(Probability)
        }))
    } else {
        if (!all(sapply(Intrv, function(a) inherits(a, c("data.frame", "data.table", "matrix"))))) 
            stop("Interventions must be appropriately named and sized data.frames or data.tables")
        ITT <- vector("list", length = length(Intrv))
        ITT <- lapply(Intrv, function(NewTrt) {
            list("intervention" = function(ObservedTrt, Covariates, PropScore) 
            {
                TrtNames <- colnames(ObservedTrt)
                Intervened <- data.table::copy(ObservedTrt)
                Intervened[, (TrtNames) := eval(NewTrt), .SDcols = TrtNames]
                return(Intervened)
            },
            "g.star" = function(Treatment, Covariates, PropScore, Intervened) 
            {
                Probability <- data.table(1 * sapply(1:nrow(Treatment), function(i) 
                    all(Treatment[i, ] == Intervened[i, ])))
                return(Probability)
            })
        })
    }
    for (i in seq_along(ITT)) {
        class(ITT[[i]]) <- union("concreteIntervention", class(ITT[[i]]))
    }
    return(ITT)
}


#' @describeIn formatArguments print.ConcreteArgs print method for "ConcreteArgs" class
#' @param x a ConcreteArgs object
#' @param ... additional arguments to be passed into print methods
#' @exportS3Method print ConcreteArgs

print.ConcreteArgs <- function(x, ...) {
    `.` <- `..a` <- NULL
    Args <- list(...)
    cat("\nObserved Data (", nrow(x$DataTable)," rows x ", ncol(x$DataTable), 
        " cols)\nUnique IDs: \"", attr(x$DataTable, "ID"), "\" (n=", 
        attr(x$DataTable, "nEff"), "),  Time-to-Event: \"", 
        attr(x$DataTable, "EventTime"), "\",  Event Type: \"", 
        attr(x$DataTable, "EventType"), "\",  Treatment", 
        ifelse(length(attr(x$DataTable, "Treatment")) > 1, "s", ""), ": \"", 
        paste0(attr(x$DataTable, "Treatment"), collapse = ", "), 
        "\"\n\n", sep = "")
    
    EventTypes <- x$DataTable[[attr(x$DataTable, "EventType")]]
    UniqueEvents <- sort(unique(EventTypes))
    Censoring <- ifelse(any(UniqueEvents <= 0), paste(UniqueEvents[UniqueEvents <= 0], sep = ","), "None")
    EventTimes <- x$DataTable[[attr(x$DataTable, "EventTime")]]
    Treatments <- x$DataTable[, .SD, .SDcols = attr(x$DataTable, "Treatment")]
    if (any(isTRUE(Args[["Verbose"]]), isTRUE(Args[["verbose"]]))) {
        cat("Events:\n")
        for (j in UniqueEvents) {
            cat(ifelse(j <= 0, paste0("Cens. ", j), paste0("Event ", j)), " : n=", sum(EventTypes == j), 
                " (", round(mean(EventTypes == j), 2), "),  [min,max] = [", 
                min(EventTimes[EventTypes == j]), ", ", max(EventTimes[EventTypes == j]), "]\n", sep = "")    
        }
        cat("\n")
        cat(ncol(Treatments), " Treatment Variable", ifelse(ncol(Treatments) > 1, "s", ""),":\n")
        if (ncol(Treatments) > 1) {
            print(head(Treatments[, .(count = .N), by = c(colnames(Treatments))], 10))
        } else {
            if (length(unique(unlist(Treatments))) <= 5) {
                cat(colnames(Treatments),": ")
                for (a in sort(unique(unlist(Treatments)))) {
                    cat(a, ": n=", sum(unlist(Treatments) == a), 
                        " (", round(mean(unlist(Treatments) == a), 2), ")   ", sep = "")
                }
                cat("\n")
            } else {
                cat("Treatment Quantiles:\n")
                print(stats::quantile(Treatments[, 1]), Args)
            }
        }
        
        cat("\n")
        cat(nrow(attr(x$DataTable, "CovNames")), "Baseline Covariates\n")
        print(head(attr(x$DataTable, "CovNames"), 4), Args)
        if (nrow(attr(x$DataTable, "CovNames")) > 4)
            cat("...", nrow(attr(x$DataTable, "CovNames")) - 4, "rows not shown")
        cat("\n")
    }
    
    cat(rep("-", 20), "\nEstimand Specification:\n")
    cat("Target Event", ifelse(length(x$TargetEvent) == 1, "", "s"), ": ", 
        paste0(sort(x$TargetEvent), collapse = ", "), "\n\n", sep = "")
    TargTimesCapped <- sort(x$TargetTime)
    if (length(TargTimesCapped) > 6) {
        TargTimesCapped <- paste0(
            paste0(sapply(head(TargTimesCapped, 3), function(tm) 
                paste0(tm, " (", sum(EventTimes > tm), "/", length(EventTimes), ")")), 
                collapse = ", "), ", ..., ",
            paste0(sapply(tail(TargTimesCapped, 3), function(tm) 
                paste0(tm, " (", sum(EventTimes > tm), "/", length(EventTimes), ")")), 
                collapse = ", "))
    } else {
        TargTimesCapped <- paste0(sapply(TargTimesCapped, function(tm) 
            paste0(tm, " (", sum(EventTimes > tm), "/", length(EventTimes), ")")), 
            collapse = ", ")
    }
    cat("Target Time", ifelse(length(x$TargetTime) == 1, "", "s"),  " (n at risk): ", 
        TargTimesCapped, "\n\n", sep = "")
    Regimes <- x$Regime
    cat('Intervention', ifelse(length(Regimes) > 1, "s", ""), "\n", sep = "")
    for (d in seq_along(Regimes)) {
        cat('  ', names(Regimes)[d], ': (', sep = "")
        for (a in 1:ncol(Regimes[[d]])) {
            cat("\"", colnames(Regimes[[d]])[a], "\" = [", 
                paste0(unlist(head(subset(Regimes[[d]], select = a), 10)), collapse = ','), ",...]", 
                ifelse(a == ncol(Regimes[[d]]), ")  -  ", ", "), sep = "")
        }
        cat("Observed Prevalence = ", 
            round(mean(sapply(1:nrow(Treatments), function(i) {
                all(Treatments[i, ] == Regimes[[d]][i, ])})), 2), "\n", sep = "")
    }
    cat("\n")
    cat(rep("-", 20), "\nEstimation Specification:\n")
    ## cross-validation
    cat(ifelse(is.null(attr(x$CVFolds, "CVArg")$strata_ids), "", "Stratified "), 
        length(x$CVFolds), "-Fold Cross Validation \n", sep = "")
    ## SL spec
    for (Trt in attr(x$DataTable, "Treatment")) {
        
        TrtMod <- x$Model[[Trt]]
        PSBackend <- attr(TrtMod, "Backend")
        if (PSBackend == "SuperLearner") {
            cat("\"", Trt,"\" Propensity Score Estimation (", 
                PSBackend, "): Default SL Selector, Default Loss Fn, ", 
                length(TrtMod), " candidate", ifelse(length(TrtMod) > 1, "s", ""), " - ", 
                paste0(head(TrtMod, 5), collapse = ", "), 
                ifelse(length(TrtMod) > 5, "...", ""), "\n", sep = "")
        } else {
            if (inherits(TrtMod, "Stack")) 
                lrnrs <- sapply(TrtMod$params$learners, 
                                function(x) sub("^Lrnr_([[:alpha:]]+)(.*)", "Lrnr_\\1", x$name))
            else
                lrnrs <- sub("^Lrnr_([[:alpha:]]+)(.*)", "Lrnr_\\1", TrtMod$name)
            cat("\"", Trt, "\" Propensity Score Estimation (", PSBackend, 
                "): Default SL Selector, Default Loss Fn, ", 
                length(lrnrs), " candidate", ifelse(length(lrnrs) > 1, "s", ""), " - ", 
                paste0(head(lrnrs, 5), collapse = ", "), 
                ifelse(length(lrnrs) > 5, "...", ""), "\n", sep = "")
        }
    }
    
    for (j in UniqueEvents) {
        JMod <- x$Model[[as.character(j)]]
        cat(ifelse(as.numeric(j) <= 0, "Cens. ", "Event "), j, 
            " Estimation (coxph): Discrete SL Selector, Log Partial-LL Loss, ", 
            length(JMod), " candidate", ifelse(length(JMod) > 1, "s", ""), sep = "")
        if (attr(JMod, "Backend") == "coxph"){
            cat(" - ", paste0(head(names(JMod), 5), collapse = ", "), 
                ifelse(length(JMod) > 5, ", ...", ""), sep = "")
        }
        cat("\n")
    }
    cat("\n")
    
    ## TMLE spec
    cat("One-step TMLE (finite sum approx.) simultaneously targeting all cause-specific Absolute Risks",
        "\ng nuisance bounds = [", signif(x$MinNuisance, 4), ", 1],  max update steps = ", 
        x$MaxUpdateIter, ",  starting one-step epsilon = ", x$OneStepEps, 
        "\n\n", sep = "")
    
    ## Misc
    if (isTRUE(attr(x$Model, "CovNamesChanged")) & x$Verbose) {
        cat("****\nCox model specifications have been renamed where necessary to reflect",
            " changed covariate names. Model specifications in .[['Model']] can be ",
            "checked against the covariate names in attr(.[['DataTable']], 'CovNames')",
            "\n****\n", sep = "")
    }
}

