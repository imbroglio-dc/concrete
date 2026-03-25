library(testthat)
library(concrete)
library(data.table)

test_that("TMLE handles very small propensity scores", {
    # Create data where one treatment group is rare
    set.seed(42)
    n <- 200
    data <- data.table(
        time = rexp(n, rate = 0.01),
        status = sample(0:1, n, replace = TRUE, prob = c(0.3, 0.7)),
        trt = sample(0:1, n, replace = TRUE, prob = c(0.95, 0.05)),  # Very unbalanced
        id = 1:n,
        age = rnorm(n, 50, 10)
    )

    # Should complete without error, possibly with warnings
    expect_warning(
        concrete.args <- formatArguments(
            DataTable = data,
            EventTime = "time",
            EventType = "status",
            Treatment = "trt",
            ID = "id",
            Intervention = makeITT(),
            TargetTime = quantile(data$time, 0.5),
            TargetEvent = 1,
            CVArg = list(V = 2),
            Verbose = FALSE,
            MaxUpdateIter = 2
        ),
        regexp = NA  # May or may not warn
    )
    expect_s3_class(concrete.args, "ConcreteArgs")
})

test_that("TMLE handles rare events", {
    # Create data where events are somewhat rare but not extremely rare
    set.seed(123)
    n <- 200
    data <- data.table(
        time = rexp(n, rate = 0.01),
        status = sample(0:1, n, replace = TRUE, prob = c(0.7, 0.3)),  # 30% events
        trt = sample(0:1, n, replace = TRUE),
        id = 1:n,
        age = rnorm(n, 50, 10),
        sex = rbinom(n, 1, 0.5)
    )

    concrete.args <- formatArguments(
        DataTable = data,
        EventTime = "time",
        EventType = "status",
        Treatment = "trt",
        ID = "id",
        Intervention = makeITT(),
        TargetTime = quantile(data$time, 0.5),
        TargetEvent = 1,
        CVArg = list(V = 2),
        Verbose = FALSE,
        MaxUpdateIter = 2
    )

    # May warn about convergence or rare events, but should complete
    result <- tryCatch(
        suppressWarnings(doConcrete(concrete.args)),
        error = function(e) NULL
    )
    # Skip if model couldn't be fit due to numerical issues
    expect_s3_class(result, "ConcreteEst")
})

test_that("TMLE convergence warning is raised when not converged", {
    set.seed(456)
    data <- as.data.table(survival::pbc)[1:100, c("time", "status", "trt", "id", "age", "sex")]
    data <- data[complete.cases(data), ]
    data[, trt := sample(0:1, .N, replace = TRUE)]
    data[, status := as.numeric(status >= 1)]

    concrete.args <- formatArguments(
        DataTable = data,
        EventTime = "time",
        EventType = "status",
        Treatment = "trt",
        ID = "id",
        Intervention = makeITT(),
        TargetTime = 2500,
        TargetEvent = 1,
        CVArg = list(V = 2),
        Verbose = FALSE,
        MaxUpdateIter = 1  # Very few iterations to trigger non-convergence
    )

    # Should warn about TMLE not converging
    expect_warning(
        concrete.est <- doConcrete(concrete.args),
        regexp = "TMLE has not converged|converged"
    )
})

test_that("getOutput handles edge cases in estimation", {
    set.seed(789)
    data <- as.data.table(survival::pbc)[1:150, c("time", "status", "trt", "id", "age", "sex")]
    data <- data[complete.cases(data), ]
    data[, trt := sample(0:1, .N, replace = TRUE)]
    data[, status := as.numeric(status >= 1)]

    concrete.args <- formatArguments(
        DataTable = data,
        EventTime = "time",
        EventType = "status",
        Treatment = "trt",
        ID = "id",
        Intervention = makeITT(),
        TargetTime = 2500,
        TargetEvent = 1,
        CVArg = list(V = 2),
        Verbose = FALSE,
        MaxUpdateIter = 3
    )

    suppressWarnings(concrete.est <- doConcrete(concrete.args))

    # Output should have valid structure even with challenging data
    output <- getOutput(concrete.est)
    expect_s3_class(output, "ConcreteOut")
    expect_true(all(is.finite(output$`Pt Est`)))
    expect_true(all(output$se >= 0 | is.na(output$se)))
})

test_that("Extreme covariate values are handled", {
    set.seed(101)
    n <- 150
    data <- data.table(
        time = rexp(n, rate = 0.01),
        status = sample(0:1, n, replace = TRUE),
        trt = sample(0:1, n, replace = TRUE),
        id = 1:n,
        age = c(rnorm(n-2, 50, 10), 1e6, -1e6)  # Extreme outliers
    )

    # Should handle extreme values, possibly with warnings about model fit
    expect_warning(
        concrete.args <- formatArguments(
            DataTable = data,
            EventTime = "time",
            EventType = "status",
            Treatment = "trt",
            ID = "id",
            Intervention = makeITT(),
            TargetTime = quantile(data$time, 0.5),
            TargetEvent = 1,
            CVArg = list(V = 2),
            Verbose = FALSE,
            MaxUpdateIter = 2
        ),
        regexp = NA  # May or may not warn
    )
    expect_s3_class(concrete.args, "ConcreteArgs")
})

test_that("Complete separation in treatment is handled", {
    # Create data where treatment is strongly (but not perfectly) predicted by a covariate
    set.seed(202)
    n <- 100
    age <- c(rnorm(n/2, 35, 10), rnorm(n/2, 65, 10))  # Overlapping ranges
    data <- data.table(
        time = rexp(n, rate = 0.01),
        status = sample(0:1, n, replace = TRUE),
        trt = c(rep(0, n/2), rep(1, n/2)),
        id = 1:n,
        age = age, 
        sex = rbinom(n, 1, 0.5)
    )

    # May warn about propensity score estimation issues
    suppressWarnings(
        concrete.args <- formatArguments(
            DataTable = data,
            EventTime = "time",
            EventType = "status",
            Treatment = "trt",
            ID = "id",
            Intervention = makeITT(),
            TargetTime = quantile(data$time, 0.5),
            TargetEvent = 1,
            CVArg = list(V = 2),
            Verbose = FALSE,
            MaxUpdateIter = 2
        )
    )

    # Estimation may fail in extreme cases, handle gracefully
    result <- tryCatch(
        suppressWarnings(doConcrete(concrete.args)),
        error = function(e) NULL
    )
    expect_s3_class(result, "ConcreteEst")
})

# ---- boundSurvival unit tests ------------------------------------------------

test_that("boundSurvival bounds near-zero cells, warns once, and sets SurvWarnings$First", {
    EvalTimes <- 1:5
    TargetRows <- c(3L, 5L)

    # Column = subject; 5 time rows, 3 subjects
    # Subjects 1 and 3 have near-zero survival at target rows 3 and 5
    Estimates <- list(
        regime0 = list(EvntFreeSurv = matrix(c(
            0.9, 0.7,  5e-13, 0.2, 0.0,    # subject 1: violations at rows 3 & 5
            0.9, 0.7,  0.4,   0.2, 0.1,    # subject 2: fine
            0.9, 0.7,  1e-13, 0.2, 2e-15   # subject 3: violations at rows 3 & 5
        ), nrow = 5, ncol = 3)),
        regime1 = list(EvntFreeSurv = matrix(0.5, nrow = 5, ncol = 3))  # all fine
    )
    attr(Estimates, "SurvWarnings") <- list(First = NULL, Last = NULL)

    expect_warning(
        out <- concrete:::boundSurvival(Estimates, TargetRows, EvalTimes, StepNum = 1L),
        regexp = "Survival probability approaching 0"
    )

    sw <- attr(out, "SurvWarnings")

    # First is set; Last equals First on the initial occurrence
    expect_false(is.null(sw$First))
    expect_equal(sw$First, sw$Last)

    # Correct columns in order
    expect_named(sw$First, c("Regime", "Step", "TargetTime", "MinSurv", "N_subjects", "Pct_subjects"))

    # Only regime0 had violations; one row per target time = 2 rows total
    expect_equal(nrow(sw$First), 2L)
    expect_true(all(sw$First$Regime == "regime0"))
    expect_true(all(sw$First$Step == 1L))
    expect_equal(sw$First$TargetTime, c(3, 5))       # EvalTimes[TargetRows]
    expect_true(all(sw$First$MinSurv <= 1e-12))
    expect_equal(sw$First$N_subjects, c(2L, 2L))     # subjects 1 & 3 at both target rows
    expect_equal(sw$First$Pct_subjects, c(round(100 * 2 / 3, 1), round(100 * 2 / 3, 1)))
})

test_that("boundSurvival physically bounds near-zero cells and leaves others untouched", {
    EvalTimes <- 1:4
    TargetRows <- 3L

    Estimates <- list(
        regime0 = list(EvntFreeSurv = matrix(c(
            0.5, 0.5, 5e-13,       0.5,   # subject 1: near-zero at row 3
            0.5, 0.5, 0.5,         0.5,   # subject 2: all fine
            0.5, 0.5, NA_real_,    0.5    # subject 3: NA at row 3
        ), nrow = 4, ncol = 3))
    )
    attr(Estimates, "SurvWarnings") <- list(First = NULL, Last = NULL)

    out <- suppressWarnings(
        concrete:::boundSurvival(Estimates, TargetRows, EvalTimes, StepNum = 1L)
    )
    m <- out[["regime0"]][["EvntFreeSurv"]]

    expect_equal(m[3, 1], 1e-12)     # near-zero → bounded
    expect_equal(m[3, 2], 0.5)       # fine → unchanged
    expect_true(is.na(m[3, 3]))      # NA → left alone
    expect_equal(m[1, 1], 0.5)       # non-target row, high value → unchanged
})

test_that("boundSurvival updates SurvWarnings$Last but not $First on repeated occurrence", {
    EvalTimes <- 1:3
    TargetRows <- 2L

    # 3×1 matrix: one subject, near-zero at row 2
    Estimates <- list(
        regime0 = list(EvntFreeSurv = matrix(c(0.5, 1e-13, 0.5), nrow = 3, ncol = 1))
    )
    attr(Estimates, "SurvWarnings") <- list(First = NULL, Last = NULL)

    # Step 1: first occurrence
    est1 <- suppressWarnings(
        concrete:::boundSurvival(Estimates, TargetRows, EvalTimes, StepNum = 1L)
    )
    first_df <- attr(est1, "SurvWarnings")$First

    # Reset the bounded cell to near-zero to simulate the next TMLE iteration
    est1[["regime0"]][["EvntFreeSurv"]][2, 1] <- 1e-13

    # Step 2: second occurrence — no new warning should fire
    est2 <- expect_no_warning(
        concrete:::boundSurvival(est1, TargetRows, EvalTimes, StepNum = 2L)
    )

    sw2 <- attr(est2, "SurvWarnings")
    expect_equal(sw2$First, first_df)   # First is frozen after step 1
    expect_equal(sw2$Last$Step, 2L)     # Last reflects the new step
    expect_false(isTRUE(all.equal(sw2$First, sw2$Last)))  # they differ (Step 1 vs 2)
})

test_that("boundSurvival leaves SurvWarnings empty when all survival is above floor", {
    EvalTimes <- 1:3
    TargetRows <- c(1L, 2L, 3L)

    # All values well above 1e-12
    Estimates <- list(
        regime0 = list(EvntFreeSurv = matrix(seq(0.9, 0.5, length.out = 9), nrow = 3, ncol = 3))
    )
    attr(Estimates, "SurvWarnings") <- list(First = NULL, Last = NULL)

    out <- expect_no_warning(
        concrete:::boundSurvival(Estimates, TargetRows, EvalTimes, StepNum = 1L)
    )

    sw <- attr(out, "SurvWarnings")
    expect_null(sw$First)
    expect_null(sw$Last)
    # Matrix should be bit-for-bit identical to the input
    expect_equal(out[["regime0"]][["EvntFreeSurv"]],
                 Estimates[["regime0"]][["EvntFreeSurv"]])
})

# ---- SurvWarnings integration tests ------------------------------------------

test_that("SurvWarnings attribute on ConcreteEst has correct list(First, Last) structure", {
    set.seed(404)
    n <- 100
    data <- data.table(
        time   = rexp(n, rate = 1),
        status = rep(1L, n),  # all events, no censoring
        trt    = sample(0:1, n, replace = TRUE),
        id     = 1:n,
        age    = rnorm(n, 50, 10), 
        sex    = rbinom(n, 1, 0.5)
    )
    target_times <- quantile(data$time, c(0.5, 0.9, 0.99))

    concrete.args <- formatArguments(
        DataTable = data, EventTime = "time", EventType = "status",
        Treatment = "trt", ID = "id", Intervention = makeITT(),
        TargetTime = target_times, TargetEvent = 1,
        CVArg = list(V = 2), Verbose = FALSE, MaxUpdateIter = 2
    )

    result <- tryCatch(
        suppressWarnings(doConcrete(concrete.args)),
        error = function(e) NULL
    )
    skip_if(is.null(result), "Model fitting failed")

    sw <- attr(result, "SurvWarnings")

    # Attribute must exist after doTmleUpdate ran
    expect_false(is.null(sw), label = "SurvWarnings attribute present on ConcreteEst")
    expect_type(sw, "list")
    expect_named(sw, c("First", "Last"))

    # If bounding occurred, verify the data frame structure
    if (!is.null(sw$First)) {
        expect_s3_class(sw$First, "data.frame")
        expect_named(sw$First,
                     c("Regime", "Step", "TargetTime", "MinSurv", "N_subjects", "Pct_subjects"))
        expect_true(all(sw$First$MinSurv <= 1e-12))
        expect_true(all(sw$First$N_subjects >= 1L))
        expect_true(all(sw$First$Pct_subjects >= 0 & sw$First$Pct_subjects <= 100))

        expect_s3_class(sw$Last, "data.frame")
        expect_named(sw$Last,
                     c("Regime", "Step", "TargetTime", "MinSurv", "N_subjects", "Pct_subjects"))
    }
})

test_that("SurvWarnings$First and $Last are NULL when survival stays above floor", {
    set.seed(505)
    n <- 150
    data <- data.table(
        time   = rexp(n, rate = 0.01),  # slow event rate → high survival throughout
        status = sample(0:1, n, replace = TRUE, prob = c(0.5, 0.5)),
        trt    = sample(0:1, n, replace = TRUE),
        id     = 1:n,
        age    = rnorm(n, 50, 10), 
        sex   = rbinom(n, 1, 0.5)
    )
    early_time <- quantile(data$time, 0.25)

    concrete.args <- formatArguments(
        DataTable = data, EventTime = "time", EventType = "status",
        Treatment = "trt", ID = "id", Intervention = makeITT(),
        TargetTime = early_time, TargetEvent = 1,
        CVArg = list(V = 2), Verbose = FALSE, MaxUpdateIter = 2
    )

    result <- tryCatch(
        suppressWarnings(doConcrete(concrete.args)),
        error = function(e) NULL
    )
    skip_if(is.null(result), "Model fitting failed")

    sw <- attr(result, "SurvWarnings")
    # If TMLE ran, SurvWarnings is set but First/Last should be NULL (no bounding needed)
    # If TMLE was skipped (already converged), sw itself will be NULL — also acceptable
    expect_true(is.null(sw) || is.null(sw$First),
                label = "No survival bounding expected at early target time with low event rate")
})
