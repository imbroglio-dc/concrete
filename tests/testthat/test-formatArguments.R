test_that("formatArguments works for a few sample analyses", {
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
                                     TargetTime = quantile(data[["time"]], probs = seq(.1, .9, .05)), 
                                     TargetEvent = unique(data[["status"]])
    )
    expect_s3_class(concrete.args, class = "ConcreteArgs")
    expect_s3_class(formatArguments(concrete.args), class = "ConcreteArgs")
    expect_error(formatArguments(ConcreteArgs = data))
    # data[, status := as.numeric(status >= 1)] # to make simple right-censored survival
    # data[status == 0, status := sample(1:2, sum(status == 0), replace = TRUE)]
})

test_that("Data with missingness or incorrect type throw errors", {
    require(data.table)
    data <- as.data.table(survival::pbc)[, c("time", "status", "trt", "id", "age", "sex")]
    expect_error(formatArguments(Data = data, 
                                 EventTime = "time", 
                                 EventType = "status", 
                                 Treatment = "trt", 
                                 ID = 'id', 
                                 Intervention = 0:1, 
                                 TargetTime = mean(data[["time"]]), 
                                 TargetEvent = unique(data[["status"]])))
    expect_error(formatArguments(Data = as.data.frame(data), 
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
    expect_error(formatArguments(Data = "blah", 
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
        data <- data.frame("x" = value)
        expect_error(checkEventTime(value, data))
        expect_error(checkEventTime("x", data))
    }
})

test_that("EventType is a non-negative numeric vector", {
    test_vals <- list(NaN, NA, TRUE, "a", -1)
    for (value in test_vals) {
        data <- data.frame("x" = value)
        expect_error(checkEventType(value, data))
        expect_error(checkEventType("x", data))
    }
})

test_that("Treatment is a numeric vector", {
    test_vals <- list(NaN, NA, Inf, TRUE, "a")
    for (value in test_vals) {
        data <- data.frame("x" = value)
        expect_error(checkTreatment(value, data))
        expect_error(checkTreatment("x", data))
    }
})

test_that("Intervention specifications", {
    require(data.table)
    data <- as.data.table(survival::pbc)[, c("time", "status", "trt", "id", "age", "sex")]
    
    set.seed(0)
    data[, trt := sample(0:1, length(trt), replace = TRUE)]
    test_vals <- list(NaN, NA, Inf, TRUE, "a", 1, matrix(1, 3, 3), function(...) return(list(...)),
                      list(function(x) x, function(y) y))
    for (value in test_vals) {
        expect_error(getRegime(value, data = data))
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

