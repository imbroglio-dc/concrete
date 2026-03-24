library(testthat)
library(concrete)
library(data.table)

# Skip tests if survivalSL is not available
skip_if_not_installed("survivalSL")

test_that("makeSurvivalSL creates valid specification", {
    spec <- makeSurvivalSL(
        methods = c("LIB_AFTgamma", "LIB_PHgompertz"),
        metric = "ibs",
        cv = 5,
        seed = 123
    )

    expect_s3_class(spec, "Lrnr.SurvivalSL")
    expect_true(spec$survivalSL)
    expect_equal(spec$methods, c("LIB_AFTgamma", "LIB_PHgompertz"))
    expect_equal(spec$metric, "ibs")
    expect_equal(spec$cv, 5)
    expect_equal(spec$seed, 123)
})

test_that("makeSurvivalSL validates methods", {
    # Should warn about invalid methods
    expect_warning(
        spec <- makeSurvivalSL(
            methods = c("LIB_AFTgamma", "INVALID_METHOD", "LIB_PHgompertz"),
            metric = "ibs"
        ),
        "Unknown survivalSL methods"
    )

    # After warning, only valid methods should remain
    expect_equal(spec$methods, c("LIB_AFTgamma", "LIB_PHgompertz"))
})

test_that("makeSurvivalSL requires at least 2 methods", {
    expect_error(
        makeSurvivalSL(methods = c("LIB_AFTgamma")),
        "at least 2 learner methods"
    )
})

test_that("makeSurvivalSL validates metric", {
    expect_error(
        makeSurvivalSL(
            methods = c("LIB_AFTgamma", "LIB_PHgompertz"),
            metric = "invalid_metric"
        ),
        "Invalid metric"
    )
})

test_that("isSurvivalSLSpec correctly identifies specs", {
    spec <- makeSurvivalSL()
    expect_true(concrete:::isSurvivalSLSpec(spec))

    # List with survivalSL = TRUE
    list_spec <- list(survivalSL = TRUE, methods = c("LIB_AFTgamma", "LIB_PHgompertz"))
    expect_true(concrete:::isSurvivalSLSpec(list_spec))

    # Regular formula should not be identified as survivalSL
    formula_spec <- Surv(time, status) ~ age + sex
    expect_false(concrete:::isSurvivalSLSpec(formula_spec))

    # Regular list without survivalSL key
    regular_list <- list(a = 1, b = 2)
    expect_false(concrete:::isSurvivalSLSpec(regular_list))
})

test_that("formatArguments accepts survivalSL model specs", {
    data <- as.data.table(survival::pbc)
    data <- data[1:100, .SD, .SDcols = c("id", "time", "status", "trt", "age", "sex")]
    data <- data[complete.cases(data), ]
    data[, trt := as.numeric(trt) - 1]
    data[, status := ifelse(status == 0, 0, 1)]

    intervention <- makeITT()
    sl_spec <- makeSurvivalSL(
        methods = c("LIB_AFTgamma", "LIB_PHgompertz"),
        metric = "ibs",
        cv = 3
    )

    model <- list(
        "trt" = c("SL.glm"),
        "0" = sl_spec,
        "1" = sl_spec
    )

    # Should not error
    expect_message(
        concrete_args <- formatArguments(
            DataTable = data,
            EventTime = "time",
            EventType = "status",
            Treatment = "trt",
            ID = "id",
            TargetTime = 2000,
            TargetEvent = 1,
            Intervention = intervention,
            CVArg = list(V = 2),
            Model = model,
            Verbose = TRUE
        ),
        "survivalSL"
    )

    # Check that model specs are preserved
    expect_s3_class(concrete_args$Model[["0"]], "Lrnr.SurvivalSL")
    expect_s3_class(concrete_args$Model[["1"]], "Lrnr.SurvivalSL")
})

test_that("survivalSL integration works in full pipeline", {
    skip_on_cran()

    data <- as.data.table(survival::pbc)
    data <- data[1:150, .SD, .SDcols = c("id", "time", "status", "trt", "age", "sex")]
    data <- data[complete.cases(data), ]
    data[, trt := as.numeric(trt) - 1]
    data[, status := ifelse(status == 0, 0, 1)]

    intervention <- makeITT()
    sl_spec <- makeSurvivalSL(
        methods = c("LIB_AFTgamma", "LIB_PHgompertz"),
        metric = "ibs",
        cv = 3,
        seed = 42
    )

    model <- list(
        "trt" = c("SL.glm"),
        "0" = sl_spec,
        "1" = sl_spec
    )

    concrete_args <- formatArguments(
        DataTable = data,
        EventTime = "time",
        EventType = "status",
        Treatment = "trt",
        ID = "id",
        TargetTime = 2000,
        TargetEvent = 1,
        Intervention = intervention,
        CVArg = list(V = 2),
        Model = model,
        Verbose = FALSE,
        MaxUpdateIter = 2
    )

    # Run estimation - may warn about convergence
    expect_warning(
        concrete_est <- doConcrete(concrete_args),
        regexp = NULL
    )

    # Check output structure
    expect_s3_class(concrete_est, "ConcreteEst")

    # Should be able to get output
    concrete_out <- getOutput(concrete_est)
    expect_s3_class(concrete_out, "ConcreteOut")
    expect_true(nrow(concrete_out) > 0)
})

test_that("survivalSL methods list is complete", {
    validMethods <- c(
        "LIB_AFTgamma", "LIB_AFTggamma", "LIB_AFTllogis",
        "LIB_AFTweibull", "LIB_COXaic", "LIB_COXall",
        "LIB_COXen", "LIB_COXlasso", "LIB_COXridge",
        "LIB_PHexponential", "LIB_PHgompertz", "LIB_PHspline",
        "LIB_PLANN", "LIB_RSF"
    )

    # All valid methods should be accepted without warning
    for (method in validMethods) {
        # Need two methods minimum, use AFTgamma as second
        if (method != "LIB_AFTgamma") {
            expect_silent(
                makeSurvivalSL(methods = c(method, "LIB_AFTgamma"))
            )
        }
    }
})
