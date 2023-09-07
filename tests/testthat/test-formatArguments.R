# set.seed(0)
# data <- data.table::as.data.table(survival::pbc)[, c("time", "status", "trt", "id", "age", "sex")]
# data[, trt := sample(0:1, length(trt), replace = TRUE)]
# 
# test_that("formatArguments works ", {
#     concrete.args <- formatArguments(Data = data, 
#                                      EventTime = "time", 
#                                      EventType = "status", 
#                                      Treatment = "trt", 
#                                      ID = 'id', 
#                                      Intervention = 0:1, 
#                                      TargetTime = quantile(data[["time"]], probs = seq(.1, .9, .05)), 
#                                      TargetEvent = unique(data[["status"]])
#     )
#     expect_s3_class(concrete.args, class = "ConcreteArgs")
#     expect_s3_class(formatArguments(concrete.args), class = "ConcreteArgs")
#     expect_s3_class(formatArguments(ConcreteArgs = concrete.args), class = "ConcreteArgs")
#     expect_error(formatArguments(ConcreteArgs = data))
# })
# 
# test_that("Data with missingness or incorrect type throw errors", {
#     require(data.table)
#     DataWithMissing <- as.data.table(survival::pbc)[, c("time", "status", "trt", "id", "age", "sex")]
#     expect_error(formatArguments(Data = DataWithMissing, 
#                                  EventTime = "time", 
#                                  EventType = "status", 
#                                  Treatment = "trt", 
#                                  ID = 'id', 
#                                  Intervention = 0:1, 
#                                  TargetTime = mean(data[["time"]]), 
#                                  TargetEvent = unique(data[["status"]])))
#     expect_error(formatArguments(Data = as.data.frame(DataWithMissing), 
#                                  EventTime = "time", 
#                                  EventType = "status", 
#                                  Treatment = "trt", 
#                                  ID = 'id', 
#                                  Intervention = 0:1, 
#                                  TargetTime = mean(data[["time"]]), 
#                                  TargetEvent = unique(data[["status"]])))
#     expect_error(formatArguments(Data = as.numeric(data$time), 
#                                  EventTime = "time", 
#                                  EventType = "status", 
#                                  Treatment = "trt", 
#                                  ID = 'id', 
#                                  Intervention = 0:1, 
#                                  TargetTime = mean(data[["time"]]), 
#                                  TargetEvent = unique(data[["status"]])))
#     expect_error(formatArguments(Data = "foo", 
#                                  EventTime = "time", 
#                                  EventType = "status", 
#                                  Treatment = "trt", 
#                                  ID = 'id', 
#                                  Intervention = 0:1, 
#                                  TargetTime = mean(data[["time"]]), 
#                                  TargetEvent = unique(data[["status"]])))
# })
# 
# test_that("EventTime is a positive, finite numeric vector", {
#     test_vals <- list(NaN, NA, Inf, TRUE, "a", 0, -1)
#     for (value in test_vals) {
#         expect_error(concrete:::checkEventTime(value, data.frame("x" = value)))
#         expect_error(concrete:::checkEventTime("x", data.frame("x" = value)))
#     }
# })
# 
# test_that("EventType is a non-negative numeric vector", {
#     test_vals <- list(NaN, NA, TRUE, "a", -1)
#     for (value in test_vals) {
#         expect_error(concrete:::checkEventType(value, data.frame("x" = value)))
#         expect_error(concrete:::checkEventType("x", data.frame("x" = value)))
#     }
# })
# 
# test_that("Treatment is a numeric vector", {
#     test_vals <- list(NaN, NA, Inf, TRUE, "a")
#     for (value in test_vals) {
#         expect_error(concrete:::checkTreatment(value, data.frame("x" = value)))
#         expect_error(concrete:::checkTreatment("x", data.frame("x" = value)))
#     }
# })
# 
# test_that("Intervention specifications", {
#     test_vals <- list(NaN, NA, Inf, "a", matrix(1, 3, 3), 
#                       function(...) return(list(...)),
#                       list(function(x) x, function(y) 1))
#     for (value in test_vals) {
#         expect_error(formatArguments(DataTable = data, EventTime = "time", 
#                                      EventType = "status", Treatment = "trt", 
#                                      ID = "id", Intervention = value))
#     }
# })
# 
# test_that("ID is a vector with non-\'null\'-type values", {
#     require(data.table)
#     data <- as.data.table(survival::pbc)[, c("time", "status", "trt", "id", "age", "sex")]
#     
#     set.seed(0)
#     data[, trt := sample(0:1, length(trt), replace = TRUE)]
#     
#     expect_error(concrete:::getID(NULL, data), regexp = NA)
#     
#     test_vals <- list(NaN, NA)
#     for (value in test_vals) {
#         data <- data.frame("x" = value)
#         expect_error(concrete:::getID(value, data))
#         expect_error(concrete:::getID("x", data))
#     }
# })
# 
# test_that("Boolean cheecks for non-boolean values and resets values to FALSE", {
#     require(data.table)
#     data <- as.data.table(survival::pbc)[, c("time", "status", "trt", "id", "age", "sex")]
#     
#     set.seed(0)
#     data[, trt := sample(0:1, length(trt), replace = TRUE)]
#     concrete.args <- formatArguments(Data = data, 
#                                      EventTime = "time", 
#                                      EventType = "status", 
#                                      Treatment = "trt", 
#                                      ID = 'id', 
#                                      Intervention = 0:1, 
#                                      TargetTime = mean(data[["time"]]), 
#                                      TargetEvent = unique(data[["status"]]),
#                                      Verbose = 2, 
#                                      GComp = NA, 
#                                      ReturnModels = "c", 
#                                      RenameCovs = Inf)
#     for (bool in c("Verbose", "GComp", "ReturnModels", "RenameCovs")) {
#         expect_equal(concrete.args[[bool]], FALSE)
#     }
# })
# 
# test_that("RenameCovs = FALSE gets processed correctly", {
#     require(data.table)
#     data <- as.data.table(survival::pbc)[, c("time", "status", "trt", "id", "age", "sex")]
#     
#     set.seed(0)
#     data[, trt := sample(0:1, length(trt), replace = TRUE)]
#     concrete.args <- formatArguments(Data = data, 
#                                      EventTime = "time", 
#                                      EventType = "status", 
#                                      Treatment = "trt", 
#                                      ID = 'id', 
#                                      Intervention = 0:1, 
#                                      TargetTime = mean(data[["time"]]), 
#                                      TargetEvent = unique(data[["status"]]),
#                                      RenameCovs = FALSE)
#     expect_equal(colnames(concrete.args$DataTable), colnames(data))
# })
