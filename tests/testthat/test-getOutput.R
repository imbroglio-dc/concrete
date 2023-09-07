# test_that("getOutput() throws expected errors", 
#           code = {
#               # competing risk
#               data <- as.data.table(survival::pbc)[, c("time", "status", "trt", "id", "age", "sex")]
#               set.seed(0)
#               data[, trt := sample(0:1, length(trt), replace = TRUE)]
#               
#               data[, status := as.numeric(status >= 1)]
#               
#               concrete.args.SL <- formatArguments(Data = data, EventTime = "time", EventType = "status",
#                                                   Treatment = "trt", ID = "id", Intervention = makeITT(),
#                                                   TargetTime = 2500, TargetEvent = NULL,
#                                                   Model = NULL, Verbose = TRUE, ReturnModels = TRUE)
#               concrete.est <- doConcrete(ConcreteArgs = concrete.args.SL)
#               
#               expect_error(object = {getOutput(concrete.args.SL)})
#               expect_error(object = {getOutput(concrete.est, Estimand = "a")})
#               expect_error(object = {getOutput(concrete.est, Intervention = "a")})
#               expect_error(object = {getOutput(concrete.est, Estimand = "RD", Intervention = 1)})
#               expect_error(object = {getOutput(concrete.est, Signif = NULL)})
#               expect_error(object = {getOutput(concrete.est, GComp = 5)})
#           }
# )
