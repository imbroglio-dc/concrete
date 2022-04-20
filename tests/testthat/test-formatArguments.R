
test_that("EventTimes is a positive, finite numeric vector", {
    test_vals <- list(NULL, NaN, NA, Inf, TRUE, "a", 0, matrix(1, ncol = 3, nrow = 3), as.list(1:3))
    for (value in test_vals) {
        expect_error(checkEventTimes(value))
    }
})

test_that("EventTypes is a non-negative numeric vector", {
    test_vals <- list(NULL, NaN, NA, TRUE, "a", -1, matrix(1, 3, 3), as.list(1:3))
    for (value in test_vals) {
        expect_error(checkEventTypes(value))
    }
})

test_that("Treatment is a numeric vector", {
    test_vals <- list(NULL, NaN, NA, Inf, TRUE, "a", matrix(1, 3, 3), as.list(1:3))
    for (value in test_vals) {
        expect_error(checkTreatment(value))
    }
})

test_that("Intervention", {
    test_vals <- list(NULL, NaN, NA, Inf, TRUE, "a", 1, matrix(1, 3, 3),
                      list(function(x) x, function(y) y))
    for (value in test_vals) {
        expect_error(checkIntervetion(value))
    }
})

test_that("CovDataTable is a matrix, data.frame, or data.table", {
    test_vals <- list(NULL, NaN, NA, Inf, TRUE, "a", 1, as.list(1:3))
    for (value in test_vals) {
        expect_error(checkCovDataTable(value))
    }
})

test_that("ID is a vector with non-\'null\'-type values", {
    test_vals <- list(NULL, NaN, NA, Inf, matrix(1, 3, 3), as.list(1:3))
    for (value in test_vals) {
        expect_error(checkID(value))
    }
})

