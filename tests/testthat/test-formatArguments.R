set.seed(0)
data <- data.table::as.data.table(survival::pbc)[, c("time", "status", "trt", "id", "age", "sex")]
data[, trt := sample(0:1, length(trt), replace = TRUE)]

test_that("formatArguments works ", {
    concrete.args <- formatArguments(Data = data,
                                     EventTime = "time",
                                     EventType = "status",
                                     Treatment = "trt",
                                     ID = 'id',
                                     Intervention = 0:1,
                                     TargetTime = quantile(data[["time"]], probs = seq(.1, .9, .05)),
                                     TargetEvent = unique(data[["status"]])
    )
    expect_s3_class(concrete.args, class = "ConcreteArgs")
    expect_s3_class(formatArguments(concrete.args), class = "ConcreteArgs")
    expect_s3_class(formatArguments(ConcreteArgs = concrete.args), class = "ConcreteArgs")
    expect_error(formatArguments(ConcreteArgs = data))
})

test_that("Data with missingness or incorrect type throw errors", {
    require(data.table)
    DataWithMissing <- as.data.table(survival::pbc)[, c("time", "status", "trt", "id", "age", "sex")]
    expect_error(formatArguments(Data = DataWithMissing,
                                 EventTime = "time",
                                 EventType = "status",
                                 Treatment = "trt",
                                 ID = 'id',
                                 Intervention = 0:1,
                                 TargetTime = mean(data[["time"]]),
                                 TargetEvent = unique(data[["status"]])))
    expect_error(formatArguments(Data = as.data.frame(DataWithMissing),
                                 EventTime = "time",
                                 EventType = "status",
                                 Treatment = "trt",
                                 ID = 'id',
                                 Intervention = 0:1,
                                 TargetTime = mean(data[["time"]]),
                                 TargetEvent = unique(data[["status"]])))
    expect_error(formatArguments(Data = as.numeric(data$time),
                                 EventTime = "time",
                                 EventType = "status",
                                 Treatment = "trt",
                                 ID = 'id',
                                 Intervention = 0:1,
                                 TargetTime = mean(data[["time"]]),
                                 TargetEvent = unique(data[["status"]])))
    expect_error(formatArguments(Data = "foo",
                                 EventTime = "time",
                                 EventType = "status",
                                 Treatment = "trt",
                                 ID = 'id',
                                 Intervention = 0:1,
                                 TargetTime = mean(data[["time"]]),
                                 TargetEvent = unique(data[["status"]])))
})

test_that("EventTime is a positive, finite numeric vector", {
    test_vals <- list(NaN, NA, Inf, TRUE, "a", 0, -1)
    for (value in test_vals) {
        expect_error(concrete:::checkEventTime(value, data.frame("x" = value)))
        expect_error(concrete:::checkEventTime("x", data.frame("x" = value)))
    }
})

test_that("EventType is a non-negative numeric vector", {
    test_vals <- list(NaN, NA, TRUE, "a", -1)
    for (value in test_vals) {
        expect_error(concrete:::checkEventType(value, data.frame("x" = value)))
        expect_error(concrete:::checkEventType("x", data.frame("x" = value)))
    }
})

test_that("Treatment is a numeric vector", {
    test_vals <- list(NaN, NA, Inf, TRUE, "a")
    for (value in test_vals) {
        expect_error(concrete:::checkTreatment(value, data.frame("x" = value)))
        expect_error(concrete:::checkTreatment("x", data.frame("x" = value)))
    }
})

test_that("Intervention specifications", {
    test_vals <- list(NaN, NA, Inf, "a", matrix(1, 3, 3),
                      function(...) return(list(...)),
                      list(function(x) x, function(y) 1))
    for (value in test_vals) {
        expect_error(formatArguments(DataTable = data, EventTime = "time",
                                     EventType = "status", Treatment = "trt",
                                     ID = "id", Intervention = value))
    }
})

test_that("ID is a vector with non-\'null\'-type values", {
    require(data.table)
    data <- as.data.table(survival::pbc)[, c("time", "status", "trt", "id", "age", "sex")]

    set.seed(0)
    data[, trt := sample(0:1, length(trt), replace = TRUE)]

    expect_error(concrete:::getID(NULL, data), regexp = NA)

    test_vals <- list(NaN, NA)
    for (value in test_vals) {
        data <- data.frame("x" = value)
        expect_error(concrete:::getID(value, data))
        expect_error(concrete:::getID("x", data))
    }
})

test_that("Boolean cheecks for non-boolean values and resets values to FALSE", {
    require(data.table)
    data <- as.data.table(survival::pbc)[, c("time", "status", "trt", "id", "age", "sex")]

    set.seed(0)
    data[, trt := sample(0:1, length(trt), replace = TRUE)]
    concrete.args <- formatArguments(Data = data,
                                     EventTime = "time",
                                     EventType = "status",
                                     Treatment = "trt",
                                     ID = 'id',
                                     Intervention = 0:1,
                                     TargetTime = mean(data[["time"]]),
                                     TargetEvent = unique(data[["status"]]),
                                     Verbose = 2,
                                     GComp = NA,
                                     ReturnModels = "c",
                                     RenameCovs = Inf)
    for (bool in c("Verbose", "GComp", "ReturnModels", "RenameCovs")) {
        expect_equal(concrete.args[[bool]], FALSE)
    }
})

test_that("RenameCovs = FALSE gets processed correctly", {
    require(data.table)
    data <- as.data.table(survival::pbc)[, c("time", "status", "trt", "id", "age", "sex")]

    set.seed(0)
    data[, trt := sample(0:1, length(trt), replace = TRUE)]
    concrete.args <- formatArguments(Data = data,
                                     EventTime = "time",
                                     EventType = "status",
                                     Treatment = "trt",
                                     ID = 'id',
                                     Intervention = 0:1,
                                     TargetTime = mean(data[["time"]]),
                                     TargetEvent = unique(data[["status"]]),
                                     RenameCovs = FALSE)
    expect_equal(colnames(concrete.args$DataTable), colnames(data))
})

test_that("Treatment column NAME conflicting with EventType values throws error", {
    require(data.table)
    # Create data where treatment column NAME matches an EventType VALUE
    # This causes model formula conflicts since the package uses EventType values as list names
    data_conflict <- data.table(
        time = c(1, 2, 3, 4, 5),
        status = c(0, 1, 1, 0, 1),  # EventType values include 0 and 1
        trt = c(0, 1, 0, 1, 0),
        id = 1:5,
        age = c(50, 60, 55, 45, 65)
    )
    # Rename treatment column to "1" which matches EventType value 1
    setnames(data_conflict, "trt", "1")
    # This should error because treatment column NAME "1" matches EventType VALUE 1
    expect_error(
        formatArguments(DataTable = data_conflict,
                        EventTime = "time",
                        EventType = "status",
                        Treatment = "1",
                        ID = "id",
                        Intervention = makeITT()),
        regexp = "name of the Treatment"
    )
})

test_that("Empty data throws appropriate error", {
    require(data.table)
    empty_data <- data.table(
        time = numeric(0),
        status = numeric(0),
        trt = numeric(0),
        id = numeric(0),
        age = numeric(0)
    )
    # Empty data should error - either during MinNuisance calculation or data validation
    expect_error(
        formatArguments(DataTable = empty_data,
                        EventTime = "time",
                        EventType = "status",
                        Treatment = "trt",
                        ID = "id",
                        Intervention = makeITT(),
                        MinNuisance = 0.05)  # Provide MinNuisance to bypass nrow(0) issue
    )
})

test_that("Data with only one row throws appropriate error or warning", {
    require(data.table)
    single_row <- data.table(
        time = 100,
        status = 1,
        trt = 1,
        id = 1,
        age = 50
    )
    # Should error because can't do CV with one observation
    expect_error(
        formatArguments(DataTable = single_row,
                        EventTime = "time",
                        EventType = "status",
                        Treatment = "trt",
                        ID = "id",
                        Intervention = makeITT(),
                        MinNuisance = 0.05)  # Provide MinNuisance to ensure error is from CV
    )
})

test_that("getDefaultCVFolds returns expected values", {
    # Test the CV fold calculation helper function
    # Formula: (nEff <= 30)*(nEff - 20) + (nEff <= 500)*10 + (nEff <= 5e3)*5 + (nEff <= 1e4)*2 + 3
    expect_equal(concrete:::getDefaultCVFolds(20), 20L)    # (20-20) + 10 + 5 + 2 + 3 = 20
    expect_equal(concrete:::getDefaultCVFolds(25), 25L)    # (25-20) + 10 + 5 + 2 + 3 = 25
    expect_equal(concrete:::getDefaultCVFolds(30), 30L)    # (30-20) + 10 + 5 + 2 + 3 = 30
    expect_equal(concrete:::getDefaultCVFolds(100), 20L)   # 0 + 10 + 5 + 2 + 3 = 20
    expect_equal(concrete:::getDefaultCVFolds(500), 20L)   # 0 + 10 + 5 + 2 + 3 = 20
    expect_equal(concrete:::getDefaultCVFolds(1000), 10L)  # 0 + 0 + 5 + 2 + 3 = 10
    expect_equal(concrete:::getDefaultCVFolds(5000), 10L)  # 0 + 0 + 5 + 2 + 3 = 10
    expect_equal(concrete:::getDefaultCVFolds(7000), 5L)   # 0 + 0 + 0 + 2 + 3 = 5
    expect_equal(concrete:::getDefaultCVFolds(10000), 5L)  # 0 + 0 + 0 + 2 + 3 = 5
    expect_equal(concrete:::getDefaultCVFolds(20000), 3L)  # 0 + 0 + 0 + 0 + 3 = 3
})

test_that("Covariate with same name as reserved column throws error", {
    require(data.table)
    # Data where a covariate has same name as EventTime column
    data_dup <- data.table(
        time = c(1, 2, 3, 4, 5),
        status = c(0, 1, 1, 0, 1),
        trt = c(0, 1, 0, 1, 0),
        id = 1:5,
        age = c(50, 60, 55, 45, 65)
    )
    # Adding a second column named "time" (duplicate)
    # This tests internal handling - data.table will rename it
    expect_s3_class(
        formatArguments(DataTable = data_dup,
                        EventTime = "time",
                        EventType = "status",
                        Treatment = "trt",
                        ID = "id",
                        Intervention = makeITT(),
                        TargetTime = 3),
        "ConcreteArgs"
    )
})

test_that("Missing required columns throw informative errors", {
    require(data.table)
    data_missing <- data.table(
        time = c(1, 2, 3, 4, 5),
        status = c(0, 1, 1, 0, 1),
        id = 1:5
        # Missing trt column
    )
    expect_error(
        formatArguments(DataTable = data_missing,
                        EventTime = "time",
                        EventType = "status",
                        Treatment = "trt",
                        ID = "id",
                        Intervention = makeITT()),
        regexp = "trt|Treatment"
    )

    # Missing EventTime column
    expect_error(
        formatArguments(DataTable = data_missing,
                        EventTime = "nonexistent",
                        EventType = "status",
                        Treatment = "id",
                        ID = "id",
                        Intervention = makeITT()),
        regexp = "nonexistent|EventTime"
    )
})

# ==============================================================================
# Shared test fixtures
# ==============================================================================

# Minimal clean dataset: one numeric covariate, no NAs, times > 0
make_test_data <- function(n = 100, seed = 42) {
    set.seed(seed)
    data.table::data.table(
        id     = 1:n,
        time   = rexp(n, rate = 0.005) + 1,   # mean ~201, well above any target
        status = sample(0:1, n, replace = TRUE),
        trt    = sample(0:1, n, replace = TRUE),
        age    = rnorm(n, 50, 10)
    )
}

# Dataset with a factor covariate (triggers one-hot encoding in getCovDataTable)
make_factor_data <- function(n = 100, seed = 42) {
    set.seed(seed)
    data.table::data.table(
        id     = 1:n,
        time   = rexp(n, rate = 0.005) + 1,
        status = sample(0:1, n, replace = TRUE),
        trt    = sample(0:1, n, replace = TRUE),
        age    = rnorm(n, 50, 10),
        sex    = factor(sample(c("f", "m"), n, replace = TRUE))
    )
}

# Competing risks dataset: status in {0, 1, 2}
make_cr_data <- function(seed = 42) {
    set.seed(seed)
    n <- 99L  # divisible by 3 for exact equal groups
    data.table::data.table(
        id     = 1:n,
        time   = rexp(n, rate = 0.005) + 1,
        status = rep(0:2, n %/% 3L),  # exactly 33 of each
        trt    = sample(0:1, n, replace = TRUE),
        age    = rnorm(n, 50, 10)
    )
}

# ==============================================================================
# Block 1: ConcreteArgs output structure
# ==============================================================================

test_that("formatArguments output has ConcreteArgs class and all expected fields", {
    dt <- make_test_data()
    args <- suppressMessages(
        formatArguments(DataTable = dt, EventTime = "time", EventType = "status",
                        Treatment = "trt", ID = "id", Intervention = makeITT(),
                        TargetTime = 100, TargetEvent = 1, CVArg = list(V = 2),
                        MinNuisance = 0.05)
    )
    expect_s3_class(args, "ConcreteArgs")
    expected_fields <- c("DataTable", "Regime", "TargetTime", "TargetEvent",
                         "CVFolds", "Model", "MaxUpdateIter", "OneStepEps",
                         "MinNuisance", "Verbose", "GComp", "ReturnModels", "RenameCovs")
    for (field in expected_fields) {
        expect_false(is.null(args[[field]]), info = paste("Missing field:", field))
    }
})

test_that("formatted DataTable carries required attributes", {
    dt <- make_test_data()
    args <- suppressMessages(
        formatArguments(DataTable = dt, EventTime = "time", EventType = "status",
                        Treatment = "trt", ID = "id", Intervention = makeITT(),
                        TargetTime = 100, CVArg = list(V = 2), MinNuisance = 0.05)
    )
    Data <- args$DataTable
    expect_equal(attr(Data, "EventTime"),  "time")
    expect_equal(attr(Data, "EventType"),  "status")
    expect_equal(attr(Data, "Treatment"),  "trt")
    expect_equal(attr(Data, "ID"),         "id")
    expect_equal(attr(Data, "nEff"),       100L)
    expect_false(is.null(attr(Data, "CovNames")))
    expect_false(is.null(attr(Data, "RenameCovs")))
})

test_that("CovNames attribute is a data.table with ColName/CovName/CovVal columns", {
    dt <- make_test_data()
    args <- suppressMessages(
        formatArguments(DataTable = dt, EventTime = "time", EventType = "status",
                        Treatment = "trt", ID = "id", Intervention = makeITT(),
                        TargetTime = 100, CVArg = list(V = 2), MinNuisance = 0.05,
                        RenameCovs = FALSE)
    )
    CovNames <- attr(args$DataTable, "CovNames")
    expect_s3_class(CovNames, "data.table")
    expect_true(all(c("ColName", "CovName", "CovVal") %in% colnames(CovNames)))
    # With RenameCovs=FALSE and one covariate 'age', ColName == CovName == "age"
    expect_equal(CovNames$ColName, "age")
    expect_equal(CovNames$CovName, "age")
})

test_that("special columns are placed first in formatted DataTable", {
    dt <- make_test_data()
    args <- suppressMessages(
        formatArguments(DataTable = dt, EventTime = "time", EventType = "status",
                        Treatment = "trt", ID = "id", Intervention = makeITT(),
                        TargetTime = 100, CVArg = list(V = 2), MinNuisance = 0.05)
    )
    cols <- colnames(args$DataTable)
    # id, time, status, trt must be the first four columns
    expect_equal(cols[seq_along(c("id", "time", "status", "trt"))],
                 c("id", "time", "status", "trt"))
})

# ==============================================================================
# Block 2: RenameCovs = TRUE / getCovDataTable
# ==============================================================================

test_that("RenameCovs=TRUE renames numeric covariate to L1 and CovNames maps back correctly", {
    dt <- make_test_data()   # single numeric covariate: 'age'
    args <- suppressMessages(
        formatArguments(DataTable = dt, EventTime = "time", EventType = "status",
                        Treatment = "trt", ID = "id", Intervention = makeITT(),
                        TargetTime = 100, CVArg = list(V = 2), MinNuisance = 0.05,
                        RenameCovs = TRUE)
    )
    Data <- args$DataTable
    expect_true("L1" %in% colnames(Data))
    expect_false("age" %in% colnames(Data))
    CovNames <- attr(Data, "CovNames")
    expect_equal(CovNames$ColName, "L1")
    expect_equal(CovNames$CovName, "age")
})

test_that("RenameCovs=TRUE one-hot encodes a factor covariate after numeric covariates", {
    dt <- make_factor_data()   # age (numeric) + sex (factor)
    args <- suppressMessages(
        formatArguments(DataTable = dt, EventTime = "time", EventType = "status",
                        Treatment = "trt", ID = "id", Intervention = makeITT(),
                        TargetTime = 100, CVArg = list(V = 2), MinNuisance = 0.05,
                        RenameCovs = TRUE)
    )
    Data    <- args$DataTable
    CovNames <- attr(Data, "CovNames")

    # Numeric covariate processed first → L1
    expect_true("L1" %in% colnames(Data))
    expect_false("age" %in% colnames(Data))
    expect_equal(CovNames[ColName == "L1", CovName], "age")

    # Factor covariate (sex) one-hot encoded → L2
    expect_true("L2" %in% colnames(Data))
    expect_false("sex" %in% colnames(Data))
    expect_equal(CovNames[ColName == "L2", CovName], "sex")

    # L2 must be a binary indicator (0/1)
    expect_true(all(Data$L2 %in% c(0, 1)))
})

test_that("RenameCovs=TRUE: all CovNames$ColName values exist as actual columns (regression test)", {
    # Regression test for the named-list cbind bug (numeric./categorical. prefix) and
    # the setnames reference-aliasing bug (original_names overwritten before CovNames built)
    dt <- make_factor_data()
    args <- suppressMessages(
        formatArguments(DataTable = dt, EventTime = "time", EventType = "status",
                        Treatment = "trt", ID = "id", Intervention = makeITT(),
                        TargetTime = 100, CVArg = list(V = 2), MinNuisance = 0.05,
                        RenameCovs = TRUE)
    )
    Data     <- args$DataTable
    CovNames <- attr(Data, "CovNames")
    expect_true(all(CovNames$ColName %in% colnames(Data)),
                info = paste("Missing cols:", paste(setdiff(CovNames$ColName, colnames(Data)), collapse=", ")))
    expect_true(all(CovNames$CovName %in% c("age", "sex")))
})

test_that("RenameCovs=TRUE warns and produces intercept-only L1 when no covariates present", {
    dt <- data.table::data.table(
        id = 1:40, time = rexp(40, 0.005) + 1,
        status = sample(0:1, 40, TRUE), trt = sample(0:1, 40, TRUE)
    )
    expect_warning(
        args <- suppressMessages(
            formatArguments(DataTable = dt, EventTime = "time", EventType = "status",
                            Treatment = "trt", ID = "id", Intervention = makeITT(),
                            TargetTime = 100, CVArg = list(V = 2), MinNuisance = 0.05,
                            RenameCovs = TRUE)
        ),
        regexp = "No covariate|intercept"
    )
    CovNames <- attr(args$DataTable, "CovNames")
    expect_equal(nrow(CovNames), 1L)
    expect_equal(CovNames$ColName,  "L1")
    expect_equal(CovNames$CovName, "(Intercept)")
    expect_true("L1" %in% colnames(args$DataTable))
})

test_that("RenameCovs=FALSE preserves original covariate column names", {
    dt <- make_factor_data()
    args <- suppressMessages(
        formatArguments(DataTable = dt, EventTime = "time", EventType = "status",
                        Treatment = "trt", ID = "id", Intervention = makeITT(),
                        TargetTime = 100, CVArg = list(V = 2), MinNuisance = 0.05,
                        RenameCovs = FALSE)
    )
    Data <- args$DataTable
    expect_true("age" %in% colnames(Data))
    expect_true("sex" %in% colnames(Data))
    expect_false("L1"  %in% colnames(Data))
    expect_false("L2"  %in% colnames(Data))
})

# ==============================================================================
# Block 3: TargetTime
# ==============================================================================

test_that("TargetTime=NULL defaults to last observed event time and emits message", {
    dt      <- make_test_data()
    max_evt <- max(dt[status == 1, time])
    expect_message(
        args <- suppressWarnings(
            formatArguments(DataTable = dt, EventTime = "time", EventType = "status",
                            Treatment = "trt", ID = "id", Intervention = makeITT(),
                            TargetTime = NULL, CVArg = list(V = 2), MinNuisance = 0.05)
        ),
        regexp = "No TargetTime|targeting"
    )
    expect_equal(args$TargetTime, max_evt)
})

test_that("TargetTime scalar value is stored without modification", {
    dt   <- make_test_data()
    args <- suppressMessages(
        formatArguments(DataTable = dt, EventTime = "time", EventType = "status",
                        Treatment = "trt", ID = "id", Intervention = makeITT(),
                        TargetTime = 100, CVArg = list(V = 2), MinNuisance = 0.05)
    )
    expect_equal(args$TargetTime, 100)
})

test_that("TargetTime vector is stored without modification", {
    dt    <- make_test_data()
    times <- c(50, 100, 200)
    args  <- suppressMessages(
        formatArguments(DataTable = dt, EventTime = "time", EventType = "status",
                        Treatment = "trt", ID = "id", Intervention = makeITT(),
                        TargetTime = times, CVArg = list(V = 2), MinNuisance = 0.05)
    )
    expect_equal(args$TargetTime, times)
})

test_that("non-positive TargetTime values throw an error", {
    dt <- make_test_data()
    for (bad in c(0, -5)) {
        expect_error(
            suppressMessages(
                formatArguments(DataTable = dt, EventTime = "time", EventType = "status",
                                Treatment = "trt", ID = "id", Intervention = makeITT(),
                                TargetTime = bad, CVArg = list(V = 2), MinNuisance = 0.05)
            ),
            regexp = "positive"
        )
    }
})

test_that("TargetTime exceeding max event time throws an error", {
    dt <- make_test_data()
    expect_error(
        suppressMessages(
            formatArguments(DataTable = dt, EventTime = "time", EventType = "status",
                            Treatment = "trt", ID = "id", Intervention = makeITT(),
                            TargetTime = 1e9, CVArg = list(V = 2), MinNuisance = 0.05)
        ),
        regexp = "Censored"
    )
})

test_that("TargetTime before first event of a type issues a warning but succeeds", {
    # Construct data where all events (status==1) occur at time >= 50
    set.seed(7)
    n <- 60
    dt_late <- data.table::data.table(
        id     = 1:n,
        time   = c(rep(10, 30), seq(50, 200, length.out = 30)),
        status = c(rep(0, 30), rep(1, 30)),
        trt    = sample(0:1, n, TRUE),
        age    = rnorm(n)
    )
    # TargetTime = 25 < 50 (first event) but > 0 and < max event time
    expect_warning(
        suppressMessages(
            formatArguments(DataTable = dt_late, EventTime = "time", EventType = "status",
                            Treatment = "trt", ID = "id", Intervention = makeITT(),
                            TargetTime = 25, CVArg = list(V = 2), MinNuisance = 0.05)
        ),
        regexp = "not yet occurred"
    )
})

# ==============================================================================
# Block 4: TargetEvent
# ==============================================================================

test_that("TargetEvent=NULL targets all non-zero event types", {
    dt   <- make_test_data()  # status in {0, 1}
    args <- suppressMessages(
        formatArguments(DataTable = dt, EventTime = "time", EventType = "status",
                        Treatment = "trt", ID = "id", Intervention = makeITT(),
                        TargetTime = 100, TargetEvent = NULL, CVArg = list(V = 2),
                        MinNuisance = 0.05)
    )
    expect_equal(args$TargetEvent, 1)
    expect_false(0 %in% args$TargetEvent)
})

test_that("TargetEvent excludes event type 0 (censoring) regardless of data", {
    dt   <- make_cr_data()   # status in {0, 1, 2}
    args <- suppressMessages(
        formatArguments(DataTable = dt, EventTime = "time", EventType = "status",
                        Treatment = "trt", ID = "id", Intervention = makeITT(),
                        TargetTime = 100, CVArg = list(V = 2), MinNuisance = 0.05)
    )
    expect_false(0 %in% args$TargetEvent)
})

test_that("TargetEvent for competing risks data includes all non-zero event types", {
    dt   <- make_cr_data()   # events: 1 and 2
    args <- suppressMessages(
        formatArguments(DataTable = dt, EventTime = "time", EventType = "status",
                        Treatment = "trt", ID = "id", Intervention = makeITT(),
                        TargetTime = 100, CVArg = list(V = 2), MinNuisance = 0.05)
    )
    expect_equal(sort(args$TargetEvent), c(1L, 2L))
})

test_that("TargetEvent input is overridden: all non-zero events always returned (documents current behavior)", {
    # NOTE: getTargetEvent unconditionally overwrites the input with all non-zero events.
    # This test documents that behavior so a future change would be caught.
    dt   <- make_cr_data()
    args <- suppressMessages(
        formatArguments(DataTable = dt, EventTime = "time", EventType = "status",
                        Treatment = "trt", ID = "id", Intervention = makeITT(),
                        TargetTime = 100, TargetEvent = 1, CVArg = list(V = 2),
                        MinNuisance = 0.05)
    )
    expect_equal(sort(args$TargetEvent), c(1L, 2L))
})

# ==============================================================================
# Block 5: Intervention / getRegime
# ==============================================================================

test_that("makeITT() produces two regimes named 'A=1' and 'A=0'", {
    dt   <- make_test_data()
    args <- suppressMessages(
        formatArguments(DataTable = dt, EventTime = "time", EventType = "status",
                        Treatment = "trt", ID = "id", Intervention = makeITT(),
                        TargetTime = 100, CVArg = list(V = 2), MinNuisance = 0.05)
    )
    expect_equal(length(args$Regime), 2L)
    expect_equal(sort(names(args$Regime)), sort(c("A=1", "A=0")))
})

test_that("integer shorthand 0:1 produces two regimes", {
    dt   <- make_test_data()
    args <- suppressMessages(
        formatArguments(DataTable = dt, EventTime = "time", EventType = "status",
                        Treatment = "trt", ID = "id", Intervention = 0:1,
                        TargetTime = 100, CVArg = list(V = 2), MinNuisance = 0.05)
    )
    expect_equal(length(args$Regime), 2L)
})

test_that("single integer 1 produces only the 'A=1' regime", {
    dt   <- make_test_data()
    args <- suppressMessages(
        formatArguments(DataTable = dt, EventTime = "time", EventType = "status",
                        Treatment = "trt", ID = "id", Intervention = 1,
                        TargetTime = 100, CVArg = list(V = 2), MinNuisance = 0.05)
    )
    expect_equal(length(args$Regime), 1L)
    expect_equal(names(args$Regime), "A=1")
})

test_that("each regime is a data.table with class concreteIntervention and a g.star function", {
    dt   <- make_test_data()
    args <- suppressMessages(
        formatArguments(DataTable = dt, EventTime = "time", EventType = "status",
                        Treatment = "trt", ID = "id", Intervention = makeITT(),
                        TargetTime = 100, CVArg = list(V = 2), MinNuisance = 0.05)
    )
    for (regime in args$Regime) {
        expect_s3_class(regime, "concreteIntervention")
        expect_s3_class(regime, "data.table")
        expect_equal(nrow(regime), nrow(dt))
        expect_true(is.function(attr(regime, "g.star")))
    }
})

test_that("A=1 regime sets all trt values to 1 and A=0 regime sets all to 0", {
    dt   <- make_test_data()
    args <- suppressMessages(
        formatArguments(DataTable = dt, EventTime = "time", EventType = "status",
                        Treatment = "trt", ID = "id", Intervention = makeITT(),
                        TargetTime = 100, CVArg = list(V = 2), MinNuisance = 0.05)
    )
    expect_true(all(args$Regime[["A=1"]]$trt == 1))
    expect_true(all(args$Regime[["A=0"]]$trt == 0))
})

test_that("custom named intervention list preserves regime names", {
    dt       <- make_test_data()
    my_int   <- makeITT()
    names(my_int) <- c("treated", "control")
    args <- suppressMessages(
        formatArguments(DataTable = dt, EventTime = "time", EventType = "status",
                        Treatment = "trt", ID = "id", Intervention = my_int,
                        TargetTime = 100, CVArg = list(V = 2), MinNuisance = 0.05)
    )
    expect_true("treated" %in% names(args$Regime))
    expect_true("control" %in% names(args$Regime))
})

test_that("custom intervention with no g.star gets default indicator function and emits message", {
    dt     <- make_test_data()
    no_gstar <- list(
        "always_one" = list(
            intervention = function(trt, cov) {
                out <- data.table::copy(trt)
                out[, trt := 1L]
                out
            }
            # g.star deliberately omitted
        )
    )
    expect_message(
        formatArguments(DataTable = dt, EventTime = "time", EventType = "status",
                        Treatment = "trt", ID = "id", Intervention = no_gstar,
                        TargetTime = 100, CVArg = list(V = 2), MinNuisance = 0.05),
        regexp = "g.star|No g.star|default"
    )
    args <- suppressMessages(
        formatArguments(DataTable = dt, EventTime = "time", EventType = "status",
                        Treatment = "trt", ID = "id", Intervention = no_gstar,
                        TargetTime = 100, CVArg = list(V = 2), MinNuisance = 0.05)
    )
    expect_true(is.function(attr(args$Regime[["always_one"]], "g.star")))
})

test_that("intervention returning out-of-observed-range values throws an error", {
    dt      <- make_test_data()
    bad_int <- list(
        list(intervention = function(trt, cov) {
            out <- data.table::copy(trt)
            out[, trt := 999L]   # outside observed range [0, 1]
            out
        })
    )
    expect_error(
        suppressMessages(
            formatArguments(DataTable = dt, EventTime = "time", EventType = "status",
                            Treatment = "trt", ID = "id", Intervention = bad_int,
                            TargetTime = 100, CVArg = list(V = 2), MinNuisance = 0.05)
        )
    )
})

test_that("g.star returning values > 1 throws an error", {
    dt       <- make_test_data()
    bad_gstar <- list(
        list(
            intervention = function(trt, cov) {
                out <- data.table::copy(trt); out[, trt := 1L]; out
            },
            g.star = function(trt, cov, ps, intervened) {
                data.table::data.table(V1 = rep(2, nrow(trt)))  # probabilities > 1
            }
        )
    )
    expect_error(
        suppressMessages(
            formatArguments(DataTable = dt, EventTime = "time", EventType = "status",
                            Treatment = "trt", ID = "id", Intervention = bad_gstar,
                            TargetTime = 100, CVArg = list(V = 2), MinNuisance = 0.05)
        )
    )
})

# ==============================================================================
# Block 6: CVArg / getCVFolds
# ==============================================================================

test_that("CVArg=NULL generates folds with V from getDefaultCVFolds", {
    dt          <- make_test_data(n = 100)
    expected_V  <- concrete:::getDefaultCVFolds(100L)   # = 20
    args        <- suppressMessages(
        formatArguments(DataTable = dt, EventTime = "time", EventType = "status",
                        Treatment = "trt", ID = "id", Intervention = makeITT(),
                        TargetTime = 100, CVArg = NULL, MinNuisance = 0.05)
    )
    expect_equal(length(args$CVFolds), expected_V)
})

test_that("CVArg=list(V=2) produces exactly 2 folds", {
    dt   <- make_test_data()
    args <- suppressMessages(
        formatArguments(DataTable = dt, EventTime = "time", EventType = "status",
                        Treatment = "trt", ID = "id", Intervention = makeITT(),
                        TargetTime = 100, CVArg = list(V = 2), MinNuisance = 0.05)
    )
    expect_equal(length(args$CVFolds), 2L)
})

test_that("each fold has non-empty training_set and validation_set integer vectors", {
    dt   <- make_test_data()
    args <- suppressMessages(
        formatArguments(DataTable = dt, EventTime = "time", EventType = "status",
                        Treatment = "trt", ID = "id", Intervention = makeITT(),
                        TargetTime = 100, CVArg = list(V = 2), MinNuisance = 0.05)
    )
    for (fold in args$CVFolds) {
        expect_true(length(fold$training_set)   > 0L)
        expect_true(length(fold$validation_set) > 0L)
    }
})

test_that("fold training and validation sets partition 1:n with no overlap", {
    dt   <- make_test_data()
    args <- suppressMessages(
        formatArguments(DataTable = dt, EventTime = "time", EventType = "status",
                        Treatment = "trt", ID = "id", Intervention = makeITT(),
                        TargetTime = 100, CVArg = list(V = 2), MinNuisance = 0.05)
    )
    n <- nrow(dt)
    for (fold in args$CVFolds) {
        # No overlap between training and validation
        expect_equal(length(intersect(fold$training_set, fold$validation_set)), 0L)
        # Together they cover all rows
        expect_equal(sort(c(fold$training_set, fold$validation_set)), seq_len(n))
    }
})

test_that("non-list CVArg emits informative error", {
    dt <- make_test_data()
    expect_error(
        suppressMessages(
            formatArguments(DataTable = dt, EventTime = "time", EventType = "status",
                            Treatment = "trt", ID = "id", Intervention = makeITT(),
                            TargetTime = 100, CVArg = "not_a_list", MinNuisance = 0.05)
        ),
        regexp = "list|CVArg"
    )
})

test_that("CVFolds carries CVArg and CVSeed attributes", {
    dt   <- make_test_data()
    args <- suppressMessages(
        formatArguments(DataTable = dt, EventTime = "time", EventType = "status",
                        Treatment = "trt", ID = "id", Intervention = makeITT(),
                        TargetTime = 100, CVArg = list(V = 2), MinNuisance = 0.05)
    )
    expect_false(is.null(attr(args$CVFolds, "CVArg")))
    expect_false(is.null(attr(args$CVFolds, "CVSeed")))
    expect_equal(attr(args$CVFolds, "CVArg")[["V"]], 2L)
})

# ==============================================================================
# Block 7: Model / makeModelList
# ==============================================================================

test_that("Model=NULL generates default entries for treatment and all event types", {
    dt   <- make_test_data()
    args <- suppressMessages(
        formatArguments(DataTable = dt, EventTime = "time", EventType = "status",
                        Treatment = "trt", ID = "id", Intervention = makeITT(),
                        TargetTime = 100, CVArg = list(V = 2), Model = NULL,
                        MinNuisance = 0.05)
    )
    expect_true("trt" %in% names(args$Model))
    expect_true("0"   %in% names(args$Model))
    expect_true("1"   %in% names(args$Model))
})

test_that("default propensity model is a character vector containing SL.xgboost and SL.glmnet", {
    dt   <- make_test_data()
    args <- suppressMessages(
        formatArguments(DataTable = dt, EventTime = "time", EventType = "status",
                        Treatment = "trt", ID = "id", Intervention = makeITT(),
                        TargetTime = 100, CVArg = list(V = 2), Model = NULL,
                        MinNuisance = 0.05)
    )
    trt_model <- args$Model[["trt"]]
    expect_true(is.character(trt_model))
    expect_true("SL.xgboost" %in% trt_model)
    expect_true("SL.glmnet"  %in% trt_model)
})

test_that("default hazard model is a list of Lrnr.Cox formulas with correct LHS", {
    dt   <- make_test_data()
    args <- suppressMessages(
        formatArguments(DataTable = dt, EventTime = "time", EventType = "status",
                        Treatment = "trt", ID = "id", Intervention = makeITT(),
                        TargetTime = 100, CVArg = list(V = 2), Model = NULL,
                        MinNuisance = 0.05)
    )
    haz <- args$Model[["1"]]
    expect_true(is.list(haz))
    expect_true(length(haz) >= 1L)
    for (m in haz) {
        expect_true(inherits(m, "Lrnr.Cox"))
        expect_true(inherits(m, "formula"))
        # LHS should reference the event time and event type (deparse gives a single string)
        lhs <- deparse(m[[2]])
        expect_true(grepl("time",   lhs))
        expect_true(grepl("status", lhs))
    }
})

test_that("custom propensity character-vector model is accepted as-is", {
    dt   <- make_test_data()
    args <- suppressMessages(
        formatArguments(DataTable = dt, EventTime = "time", EventType = "status",
                        Treatment = "trt", ID = "id", Intervention = makeITT(),
                        TargetTime = 100, CVArg = list(V = 2), MinNuisance = 0.05,
                        Model = list("trt" = c("SL.glm", "SL.mean")))
    )
    # Strip the Backend attribute that makeModelList attaches before comparing
    expect_equal(as.character(args$Model[["trt"]]), c("SL.glm", "SL.mean"))
})

test_that("unrecognized propensity model spec is replaced with default and emits message", {
    dt <- make_test_data()
    expect_message(
        args <- formatArguments(DataTable = dt, EventTime = "time", EventType = "status",
                                Treatment = "trt", ID = "id", Intervention = makeITT(),
                                TargetTime = 100, CVArg = list(V = 2), MinNuisance = 0.05,
                                Model = list("trt" = 42)),  # not a character vector
        regexp = "not recognized|replaced|default"
    )
    expect_true(is.character(args$Model[["trt"]]))
})

test_that("custom hazard formula gets Lrnr.Cox class and correct LHS added", {
    dt   <- make_test_data()
    args <- suppressMessages(
        formatArguments(DataTable = dt, EventTime = "time", EventType = "status",
                        Treatment = "trt", ID = "id", Intervention = makeITT(),
                        TargetTime = 100, CVArg = list(V = 2), MinNuisance = 0.05,
                        RenameCovs = FALSE,
                        Model = list("trt" = c("SL.glm"),
                                     "0"   = list(Surv(time, status == 0) ~ trt),
                                     "1"   = list(Surv(time, status == 1) ~ trt + age)))
    )
    expect_true(inherits(args$Model[["1"]][[1]], "Lrnr.Cox"))
    expect_true(inherits(args$Model[["1"]][[1]], "formula"))
})

test_that("'coxnet' string in hazard model gets Lrnr.Coxnet class", {
    dt   <- make_test_data()
    args <- suppressMessages(
        formatArguments(DataTable = dt, EventTime = "time", EventType = "status",
                        Treatment = "trt", ID = "id", Intervention = makeITT(),
                        TargetTime = 100, CVArg = list(V = 2), MinNuisance = 0.05,
                        Model = list("trt" = c("SL.glm"),
                                     "0"   = list("coxnet"),
                                     "1"   = list("coxnet")))
    )
    expect_true(inherits(args$Model[["1"]][[1]], "Lrnr.Coxnet"))
})

test_that("unrecognized Model key emits message (ignored during formula processing)", {
    dt <- make_test_data()
    expect_message(
        formatArguments(DataTable = dt, EventTime = "time", EventType = "status",
                        Treatment = "trt", ID = "id", Intervention = makeITT(),
                        TargetTime = 100, CVArg = list(V = 2), MinNuisance = 0.05,
                        Model = list("trt"         = c("SL.glm"),
                                     "unknown_key" = list(Surv(time, status == 1) ~ .))),
        regexp = "ignored|unknown_key"
    )
})

test_that("RenameCovs=TRUE rewrites covariate name in custom hazard formula", {
    dt   <- make_test_data()   # covariate 'age' → renamed to 'L1'
    args <- suppressMessages(
        formatArguments(DataTable = dt, EventTime = "time", EventType = "status",
                        Treatment = "trt", ID = "id", Intervention = makeITT(),
                        TargetTime = 100, CVArg = list(V = 2), MinNuisance = 0.05,
                        RenameCovs = TRUE,
                        Model = list("trt" = c("SL.glm"),
                                     "0"   = list(Surv(time, status == 0) ~ trt + age),
                                     "1"   = list(Surv(time, status == 1) ~ trt + age)))
    )
    formula_rhs <- paste(as.character(args$Model[["1"]][[1]]), collapse = " ")
    expect_true( grepl("L1",  formula_rhs))
    expect_false(grepl("\\bage\\b", formula_rhs))
})

test_that("RenameCovs=FALSE leaves covariate names unchanged in custom hazard formula", {
    dt   <- make_test_data()
    args <- suppressMessages(
        formatArguments(DataTable = dt, EventTime = "time", EventType = "status",
                        Treatment = "trt", ID = "id", Intervention = makeITT(),
                        TargetTime = 100, CVArg = list(V = 2), MinNuisance = 0.05,
                        RenameCovs = FALSE,
                        Model = list("trt" = c("SL.glm"),
                                     "0"   = list(Surv(time, status == 0) ~ trt + age),
                                     "1"   = list(Surv(time, status == 1) ~ trt + age)))
    )
    formula_rhs <- paste(as.character(args$Model[["1"]][[1]]), collapse = " ")
    expect_true(grepl("age", formula_rhs))
    expect_false(grepl("L1", formula_rhs))
})

# ==============================================================================
# Block 8: Scalar parameter validation
# ==============================================================================

test_that("MaxUpdateIter valid integer is stored exactly", {
    dt   <- make_test_data()
    args <- suppressMessages(
        formatArguments(DataTable = dt, EventTime = "time", EventType = "status",
                        Treatment = "trt", ID = "id", Intervention = makeITT(),
                        TargetTime = 100, CVArg = list(V = 2), MinNuisance = 0.05,
                        MaxUpdateIter = 50)
    )
    expect_equal(args$MaxUpdateIter, 50L)
})

test_that("MaxUpdateIter fractional value is ceiling'd", {
    dt   <- make_test_data()
    args <- suppressMessages(
        formatArguments(DataTable = dt, EventTime = "time", EventType = "status",
                        Treatment = "trt", ID = "id", Intervention = makeITT(),
                        TargetTime = 100, CVArg = list(V = 2), MinNuisance = 0.05,
                        MaxUpdateIter = 2.7)
    )
    expect_equal(args$MaxUpdateIter, 3L)
})

test_that("invalid MaxUpdateIter values reset to 100 with a message", {
    dt <- make_test_data()
    for (bad in list(-1, 0, Inf, "a", NA, NULL)) {
        expect_message(
            args <- formatArguments(DataTable = dt, EventTime = "time", EventType = "status",
                                    Treatment = "trt", ID = "id", Intervention = makeITT(),
                                    TargetTime = 100, CVArg = list(V = 2), MinNuisance = 0.05,
                                    MaxUpdateIter = bad),
            regexp = "MaxUpdateIter|100",
            info   = paste("bad value:", deparse(bad))
        )
        expect_equal(args$MaxUpdateIter, 100L, info = paste("bad value:", deparse(bad)))
    }
})

test_that("OneStepEps valid values in (0, 1] are stored exactly", {
    dt <- make_test_data()
    for (val in c(0.1, 0.5, 1.0)) {
        args <- suppressMessages(
            formatArguments(DataTable = dt, EventTime = "time", EventType = "status",
                            Treatment = "trt", ID = "id", Intervention = makeITT(),
                            TargetTime = 100, CVArg = list(V = 2), MinNuisance = 0.05,
                            OneStepEps = val)
        )
        expect_equal(args$OneStepEps, val, info = paste("val:", val))
    }
})

test_that("invalid OneStepEps values reset to 0.5 with a message", {
    dt <- make_test_data()
    for (bad in list(0, -0.1, 1.5, "a", NA)) {
        expect_message(
            args <- formatArguments(DataTable = dt, EventTime = "time", EventType = "status",
                                    Treatment = "trt", ID = "id", Intervention = makeITT(),
                                    TargetTime = 100, CVArg = list(V = 2), MinNuisance = 0.05,
                                    OneStepEps = bad),
            regexp = "OneStepEps|0.5",
            info   = paste("bad value:", deparse(bad))
        )
        expect_equal(args$OneStepEps, 0.5, info = paste("bad value:", deparse(bad)))
    }
})

test_that("MinNuisance=NULL is auto-computed as 5/sqrt(n)/log(n)", {
    dt       <- make_test_data(n = 100)
    args     <- suppressMessages(
        formatArguments(DataTable = dt, EventTime = "time", EventType = "status",
                        Treatment = "trt", ID = "id", Intervention = makeITT(),
                        TargetTime = 100, CVArg = list(V = 2), MinNuisance = NULL)
    )
    expected <- 5 / sqrt(100) / log(100)
    expect_equal(args$MinNuisance, expected)
})

test_that("auto-computed MinNuisance decreases with increasing sample size", {
    dt_s  <- make_test_data(n = 100)
    dt_l  <- make_test_data(n = 500, seed = 99)
    args_s <- suppressMessages(
        formatArguments(DataTable = dt_s, EventTime = "time", EventType = "status",
                        Treatment = "trt", ID = "id", Intervention = makeITT(),
                        TargetTime = 100, CVArg = list(V = 2), MinNuisance = NULL)
    )
    args_l <- suppressMessages(
        formatArguments(DataTable = dt_l, EventTime = "time", EventType = "status",
                        Treatment = "trt", ID = "id", Intervention = makeITT(),
                        TargetTime = 100, CVArg = list(V = 2), MinNuisance = NULL)
    )
    expect_gt(args_s$MinNuisance, args_l$MinNuisance)
})

test_that("invalid MinNuisance values reset to 0.05 with a message", {
    dt <- make_test_data()
    for (bad in list(-0.1, 0, 1.5, "a", NA)) {
        expect_message(
            args <- formatArguments(DataTable = dt, EventTime = "time", EventType = "status",
                                    Treatment = "trt", ID = "id", Intervention = makeITT(),
                                    TargetTime = 100, CVArg = list(V = 2),
                                    MinNuisance = bad),
            regexp = "MinNuisance|0.05",
            info   = paste("bad value:", deparse(bad))
        )
        expect_equal(args$MinNuisance, 0.05, info = paste("bad value:", deparse(bad)))
    }
})

# ==============================================================================
# Block 9: Competing risks
# ==============================================================================

test_that("competing risks data generates Model entries for all three event types", {
    dt   <- make_cr_data()   # events: 0, 1, 2
    args <- suppressMessages(
        formatArguments(DataTable = dt, EventTime = "time", EventType = "status",
                        Treatment = "trt", ID = "id", Intervention = makeITT(),
                        TargetTime = 100, CVArg = list(V = 2), MinNuisance = 0.05)
    )
    expect_true("0" %in% names(args$Model))
    expect_true("1" %in% names(args$Model))
    expect_true("2" %in% names(args$Model))
})

test_that("competing risks: both non-censoring event types are targeted by default", {
    dt   <- make_cr_data()
    args <- suppressMessages(
        formatArguments(DataTable = dt, EventTime = "time", EventType = "status",
                        Treatment = "trt", ID = "id", Intervention = makeITT(),
                        TargetTime = 100, CVArg = list(V = 2), MinNuisance = 0.05)
    )
    expect_equal(sort(args$TargetEvent), c(1L, 2L))
    expect_false(0L %in% args$TargetEvent)
})

test_that("competing risks: hazard models for events 1 and 2 use correct event indicator", {
    dt   <- make_cr_data()
    args <- suppressMessages(
        formatArguments(DataTable = dt, EventTime = "time", EventType = "status",
                        Treatment = "trt", ID = "id", Intervention = makeITT(),
                        TargetTime = 100, CVArg = list(V = 2), MinNuisance = 0.05)
    )
    # as.character() on a call returns a vector; use deparse() for a single string
    lhs_1 <- deparse(args$Model[["1"]][[1]][[2]])
    lhs_2 <- deparse(args$Model[["2"]][[1]][[2]])
    expect_true(grepl("== 1", lhs_1))
    expect_true(grepl("== 2", lhs_2))
})

# ==============================================================================
# Block 10: Data input flexibility and re-validation
# ==============================================================================

test_that("data.frame input is coerced to data.table and succeeds", {
    df   <- as.data.frame(make_test_data())
    args <- suppressMessages(
        formatArguments(DataTable = df, EventTime = "time", EventType = "status",
                        Treatment = "trt", ID = "id", Intervention = makeITT(),
                        TargetTime = 100, CVArg = list(V = 2), MinNuisance = 0.05)
    )
    expect_s3_class(args,            "ConcreteArgs")
    expect_s3_class(args$DataTable,  "data.table")
})

test_that("partial argument name 'Data' matches parameter 'DataTable'", {
    dt   <- make_test_data()
    args <- suppressMessages(
        formatArguments(Data = dt, EventTime = "time", EventType = "status",
                        Treatment = "trt", ID = "id", Intervention = makeITT(),
                        TargetTime = 100, CVArg = list(V = 2), MinNuisance = 0.05)
    )
    expect_s3_class(args, "ConcreteArgs")
})

test_that("ID=NULL auto-generates a sequential integer 'ID' column", {
    dt_noid <- make_test_data()[, !"id"]
    args    <- suppressMessages(
        formatArguments(DataTable = dt_noid, EventTime = "time", EventType = "status",
                        Treatment = "trt", ID = NULL, Intervention = makeITT(),
                        TargetTime = 100, CVArg = list(V = 2), MinNuisance = 0.05)
    )
    expect_equal(attr(args$DataTable, "ID"), "ID")
    expect_equal(args$DataTable$ID, seq_len(nrow(dt_noid)))
})

test_that("ConcreteArgs re-validation via positional and named argument both succeed", {
    dt    <- make_test_data()
    args1 <- suppressMessages(
        formatArguments(DataTable = dt, EventTime = "time", EventType = "status",
                        Treatment = "trt", ID = "id", Intervention = makeITT(),
                        TargetTime = 100, CVArg = list(V = 2), MinNuisance = 0.05)
    )
    # Positional
    args2 <- suppressMessages(formatArguments(args1))
    expect_s3_class(args2, "ConcreteArgs")
    # Named
    args3 <- suppressMessages(formatArguments(ConcreteArgs = args1))
    expect_s3_class(args3, "ConcreteArgs")
})

test_that("passing a non-ConcreteArgs object as ConcreteArgs throws an error", {
    expect_error(
        formatArguments(ConcreteArgs = list(a = 1, b = 2)),
        regexp = "ConcreteArgs"
    )
})
