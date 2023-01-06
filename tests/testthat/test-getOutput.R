test_that("getOutput() returns with competing risks", 
          code = {
              # competing risk
              expect_error(object = {
                  data <- as.data.table(survival::pbc)[, c("time", "status", "trt", "id", "age", "sex")]
                  set.seed(0)
                  data[, trt := sample(0:1, length(trt), replace = TRUE)]
                  
                  concrete.args.SL <- try(
                      formatArguments(Data = data, EventTime = "time", EventType = "status",
                                      Treatment = "trt", ID = "id", Intervention = makeITT(),
                                      TargetTime = 2500, TargetEvent = NULL,
                                      Model = NULL, Verbose = TRUE, ReturnModels = TRUE)
                  )
                  
                  concrete.ests <- doConcrete(ConcreteArgs = concrete.args.SL)
                  concrete.out <- getOutput(concrete.ests)
              }, regexp = NA)
          }
)

test_that("doConcrete() runs right-censored survival", 
          code = {
              # competing risk
              expect_error(object = {
                  data <- as.data.table(survival::pbc)[, c("time", "status", "trt", "id", "age", "sex")]
                  set.seed(0)
                  data[, trt := sample(0:1, length(trt), replace = TRUE)]
                  
                  data[, status := as.numeric(status >= 1)]
                  
                  concrete.args.SL <- try(
                      formatArguments(Data = data, EventTime = "time", EventType = "status",
                                      Treatment = "trt", ID = "id", Intervention = makeITT(),
                                      TargetTime = 2500, TargetEvent = NULL,
                                      Model = NULL, Verbose = TRUE, ReturnModels = TRUE)
                  )
                  
                  concrete.ests <- doConcrete(ConcreteArgs = concrete.args.SL)
                  concrete.out <- getOutput(concrete.ests)
              }, regexp = NA)
          }
)