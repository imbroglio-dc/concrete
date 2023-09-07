# test_that("doConcrete() runs with competing risks",
#           code = {
#               # competing risk
#               PnEIC <- NULL
#               data <- as.data.table(survival::pbc)[, c("time", "status", "trt", "id", "age", "sex")]
#               set.seed(0)
#               data[, trt := sample(0:1, length(trt), replace = TRUE)]
#               data <- data[1:200, ]
# 
#               # using superlearner
#               expect_error(object = {
#                   concrete.args.SL <- formatArguments(Data = data, EventTime = "time", EventType = "status",
#                                                       Treatment = "trt", ID = "id", Intervention = 0:1,
#                                                       TargetTime = 2500, TargetEvent = NULL, CVArg = list(V=2),
#                                                       MaxUpdateIter = 2, Model = NULL)
#                   concrete.args.SL$Model$trt <- c("SL.glmnet")
#                   formatArguments(concrete.args.SL)
#                   concrete.est <- doConcrete(ConcreteArgs = concrete.args.SL)
#                   print(concrete.est, Verbose = FALSE)
#                   plot(concrete.est, convergence = FALSE, propscores = TRUE, ask = FALSE)
#                   }, regexp = NA)
# 
#               # getOutput
#               expect_error(object = {
#                   concrete.out <- getOutput(concrete.est)
#                   plot(concrete.out, "RR", NullLine = TRUE, GComp = TRUE, ask = FALSE)
#                   plot(concrete.out, "RR", NullLine = TRUE, GComp = FALSE, ask = FALSE)
#                   plot(concrete.out, "RR", NullLine = FALSE, GComp = TRUE, ask = FALSE)
#                   plot(concrete.out, "RR", NullLine = FALSE, GComp = FALSE, ask = FALSE)
# 
#                   plot(concrete.out, "RD", NullLine = TRUE, GComp = TRUE, ask = FALSE)
#                   plot(concrete.out, "RD", NullLine = TRUE, GComp = FALSE, ask = FALSE)
#                   plot(concrete.out, "RD", NullLine = FALSE, GComp = TRUE, ask = FALSE)
#                   plot(concrete.out, "RD", NullLine = FALSE, GComp = FALSE, ask = FALSE)
#               }, regexp = NA)
# 
#               # not-converged tmle
#               expect_warning(object = {
#                   concrete.args.SL$MaxUpdateIter <- 1
#                   formatArguments(concrete.args.SL)
#                   doConcrete(concrete.args.SL)
#               })
# 
#               # getNormPnEIC using sigma - needs stricter tests to actually be used
#               expect_error(object = {
#                   concrete:::getNormPnEIC(concrete.est$`A=1`$SummEIC[, PnEIC], NaN)
#               }, regexp = NA)
# 
#               # using sl3
#               # expect_error(object = {
#               #     if (requireNamespace("sl3", quietly = TRUE) & requireNamespace("Rsolnp", quietly = TRUE)) {
#               #         a_lrnrs <- make_learner(Stack, Lrnr_glm$new(), Lrnr_glmnet$new())
#               #         concrete.args.sl3 <- formatArguments(Data = data, EventTime = "time", EventType = "status",
#               #                                              Treatment = "trt", ID = "id", Intervention = 0:1,
#               #                                              TargetTime = 2500, TargetEvent = NULL,
#               #                                              Model = NULL, Verbose = TRUE, ReturnModels = TRUE)
#               #         concrete.args.sl3$Model$trt <- a_lrnrs
#               #         concrete.args.sl3 <- formatArguments(concrete.args.sl3)
#               #         concrete.est <- doConcrete(ConcreteArgs = concrete.args.sl3)
#               #         print(concrete.est)
#               #         plot(concrete.est, convergence = FALSE, propscores = TRUE, ask = FALSE)
#               #     }
#               # }, regexp = NA)
#           }
# )
# 
# 
# test_that("doConcrete() runs right-censored survival",
#           code = {
#               # competing risk
#               expect_error(object = {
#                   data <- as.data.table(survival::pbc)[, c("time", "status", "trt", "id", "age", "sex")]
#                   set.seed(0)
#                   data[, trt := sample(0:1, length(trt), replace = TRUE)]
#                   data[, status := as.numeric(status >= 1)]
#                   data <- data[1:200, ]
# 
#                   concrete.args.SL <- formatArguments(Data = data, EventTime = "time", EventType = "status",
#                                                       Treatment = "trt", ID = "id", Intervention = 0:1,
#                                                       TargetTime = 2500, TargetEvent = NULL, CVArg = list(V=2),
#                                                       MaxUpdateIter = 2, Model = NULL)
#                   concrete.est <- doConcrete(ConcreteArgs = concrete.args.SL)
#                   print(concrete.est, Verbose = FALSE)
#                   plot(concrete.est, convergence = TRUE, propscores = FALSE, ask = FALSE)
# 
#                   # using sl3
#                   # if (requireNamespace("sl3", quietly = TRUE) & requireNamespace("Rsolnp", quietly = TRUE)) {
#                   #     a_lrnrs <- make_learner(Stack, Lrnr_glm$new(), Lrnr_glmnet$new())
#                   # 
#                   #     concrete.args.sl3 <- formatArguments(Data = data, EventTime = "time", EventType = "status",
#                   #                                          Treatment = "trt", ID = "id", Intervention = 0:1,
#                   #                                          TargetTime = 2500, TargetEvent = NULL, MaxUpdateIter = 2,
#                   #                                          Model = NULL, Verbose = TRUE, ReturnModels = TRUE)
#                   #     concrete.args.sl3$Model$trt <- a_lrnrs
#                   #     concrete.args.sl3 <- formatArguments(concrete.args.sl3)
#                   #     concrete.est <- doConcrete(ConcreteArgs = concrete.args.sl3)
#                   #     print(concrete.est)
#                   #     plot(concrete.est, convergence = TRUE, propscores = FALSE, ask = FALSE)
#                   # }
#               }, regexp = NA)
#           }
# )