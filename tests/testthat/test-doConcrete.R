test_that("doConcrete() runs with competing risks", 
          code = {
              # competing risk
              expect_error(object = {
                  data <- as.data.table(survival::pbc)[, c("time", "status", "trt", "id", "age", "sex")]
                  set.seed(0)
                  data[, trt := sample(0:1, length(trt), replace = TRUE)]
                  
                  # data[, status := as.numeric(status >= 1)] # competing risk or not
                  # data[status == 0, status := sample(1:2, sum(status == 0), replace = TRUE)]
                  
                  library(sl3)
                  a_lrnrs <- make_learner(Stack, Lrnr_glm$new(), Lrnr_glmnet$new())
                  
                  concrete.args.SL <- try(
                      formatArguments(Data = data, EventTime = "time", EventType = "status",
                                      Treatment = "trt", ID = "id", Intervention = makeITT(),
                                      TargetTime = 2500, TargetEvent = NULL,
                                      Model = NULL, Verbose = TRUE, ReturnModels = TRUE)
                  )
                  
                  concrete.args.sl3 <- concrete.args.SL
                  concrete.args.sl3[["Model"]][["trt"]] <- a_lrnrs
                  concrete.args.sl3[["PropScoreBackend"]] <- "sl3"
                  
                  concrete.ests <- list("SL" = try(doConcrete(ConcreteArgs = concrete.args.SL)), 
                                        "sl3" = try(doConcrete(ConcreteArgs = concrete.args.sl3)))
              }, regexp = NA)
          }
)


test_that("doConcrete() runs right-censored survival", 
          code = {
              # competing risk
              expect_snapshot_value(x = {
                  data <- as.data.table(survival::pbc)[, c("time", "status", "trt", "id", "age", "sex")]
                  set.seed(0)
                  data[, trt := sample(0:1, length(trt), replace = TRUE)]
                  
                  data[, status := as.numeric(status >= 1)] # competing risk or not
                  # data[status == 0, status := sample(1:2, sum(status == 0), replace = TRUE)]
                  
                  library(sl3)
                  a_lrnrs <- make_learner(Stack, Lrnr_glm$new(), Lrnr_glmnet$new())
                  
                  concrete.args.SL <- try(
                      formatArguments(Data = data, EventTime = "time", EventType = "status",
                                      Treatment = "trt", ID = "id", Intervention = makeITT(),
                                      TargetTime = 2500, TargetEvent = NULL,
                                      Model = NULL, Verbose = TRUE, ReturnModels = TRUE)
                  )
                  
                  concrete.args.sl3 <- concrete.args.SL
                  concrete.args.sl3[["Model"]][["trt"]] <- a_lrnrs
                  concrete.args.sl3[["PropScoreBackend"]] <- "sl3"
                  
                  concrete.ests <- list("SL" = try(doConcrete(ConcreteArgs = concrete.args.SL)), 
                                        "sl3" = try(doConcrete(ConcreteArgs = concrete.args.sl3)))
              })
          }
)
