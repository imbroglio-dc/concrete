test_that("doConcrete() runs with competing risks", 
          code = {
              # competing risk
              expect_error(object = {
                  require(data.table)
                  data <- as.data.table(survival::pbc)[, c("time", "status", "trt", "id", "age", "sex")]
                  set.seed(0)
                  data[, trt := sample(0:1, length(trt), replace = TRUE)]
                  
                  # data[, status := as.numeric(status >= 1)] # competing risk or not
                  # data[status == 0, status := sample(1:2, sum(status == 0), replace = TRUE)]
                  
                  concrete.args.SL <- formatArguments(Data = data, EventTime = "time", EventType = "status",
                                                      Treatment = "trt", ID = "id", Intervention = 0:1,
                                                      TargetTime = 2500, TargetEvent = NULL,
                                                      Model = NULL, Verbose = TRUE, ReturnModels = TRUE)
                  
                  concrete.est <- doConcrete(ConcreteArgs = concrete.args.SL)
                  print(concrete.est, Verbose = FALSE)
                  plot(concrete.est, convergence = FALSE, propscores = TRUE, ask = FALSE)
                  
                  if(require(sl3) & require(Rsolnp)){
                      a_lrnrs <- make_learner(Stack, Lrnr_glm$new(), Lrnr_glmnet$new())
                      concrete.args.sl3 <- formatArguments(Data = data, EventTime = "time", EventType = "status",
                                                           Treatment = "trt", ID = "id", Intervention = 0:1,
                                                           TargetTime = 2500, TargetEvent = NULL,
                                                           Model = NULL, Verbose = TRUE, ReturnModels = TRUE)
                      concrete.args.sl3$Model$trt <- a_lrnrs
                      concrete.args.sl3$PropScoreBackend <- "sl3"
                      concrete.args.sl3 <- formatArguments(concrete.args.sl3)
                      concrete.est <- doConcrete(ConcreteArgs = concrete.args.sl3)
                      print(concrete.est)
                      plot(concrete.est, convergence = FALSE, propscores = TRUE, ask = FALSE)
                  }
                  
              }, regexp = NA)
          }
)


test_that("doConcrete() runs right-censored survival", 
          code = {
              # competing risk
              expect_error(object = {
                  require(data.table)
                  data <- as.data.table(survival::pbc)[, c("time", "status", "trt", "id", "age", "sex")]
                  set.seed(0)
                  data[, trt := sample(0:1, length(trt), replace = TRUE)]
                  
                  data[, status := as.numeric(status >= 1)] # competing risk or not
                  # data[status == 0, status := sample(1:2, sum(status == 0), replace = TRUE)]
                  
                  
                  concrete.args.SL <- formatArguments(Data = data, EventTime = "time", EventType = "status",
                                                      Treatment = "trt", ID = "id", Intervention = 0:1,
                                                      TargetTime = 2500, TargetEvent = NULL,
                                                      Model = NULL, Verbose = TRUE, ReturnModels = TRUE)
                  concrete.est <- doConcrete(ConcreteArgs = concrete.args.SL)
                  print(concrete.est, Verbose = FALSE)
                  plot(concrete.est, convergence = TRUE, propscores = FALSE, ask = FALSE)
                  
                  if (require(sl3) & require(Rsolnp)) {
                      a_lrnrs <- make_learner(Stack, Lrnr_glm$new(), Lrnr_glmnet$new())
                      
                      concrete.args.sl3 <- formatArguments(Data = data, EventTime = "time", EventType = "status",
                                                           Treatment = "trt", ID = "id", Intervention = 0:1,
                                                           TargetTime = 2500, TargetEvent = NULL,
                                                           Model = NULL, Verbose = TRUE, ReturnModels = TRUE)
                      concrete.args.sl3$Model$trt <- a_lrnrs
                      concrete.args.sl3$PropScoreBackend <- "sl3"
                      concrete.args.sl3 <- formatArguments(concrete.args.sl3)
                      concrete.est <- doConcrete(ConcreteArgs = concrete.args.sl3)
                      print(concrete.est)
                      plot(concrete.est, convergence = TRUE, propscores = FALSE, ask = FALSE)
                  }
              }, regexp = NA)
          }
)

test_that("doConcrete() throws an error if input is not a ConcreteArgs object", 
          code = {
              # competing risk
              expect_error(object = {
                  require(data.table)
                  data <- as.data.table(survival::pbc)[, c("time", "status", "trt", "id", "age", "sex")]
                  set.seed(0)
                  data[, trt := sample(0:1, length(trt), replace = TRUE)]
                  
                  # data[, status := as.numeric(status >= 1)] # competing risk or not
                  # data[status == 0, status := sample(1:2, sum(status == 0), replace = TRUE)]
                  
                  
                  concrete.args.SL <- formatArguments(Data = data, EventTime = "time", EventType = "status",
                                                      Treatment = "trt", ID = "id", Intervention = 0:1,
                                                      TargetTime = 2500, TargetEvent = NULL,
                                                      Model = NULL, Verbose = TRUE, ReturnModels = TRUE)
                  concrete.est <- doConcrete(ConcreteArgs = data)
              })
          }
)
