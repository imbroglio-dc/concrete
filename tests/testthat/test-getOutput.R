test_that("getOutput() throws expected errors",
          code = {
              # competing risk
              data <- as.data.table(survival::pbc)[, c("time", "status", "trt", "id", "age", "sex")]
              set.seed(0)
              data[, trt := sample(0:1, length(trt), replace = TRUE)]

              data[, status := as.numeric(status >= 1)]

              concrete.args.SL <- formatArguments(Data = data, EventTime = "time", EventType = "status",
                                                  Treatment = "trt", ID = "id", Intervention = makeITT(),
                                                  TargetTime = 2500, TargetEvent = NULL,
                                                  Model = NULL, Verbose = TRUE, ReturnModels = TRUE)
              concrete.est <- doConcrete(ConcreteArgs = concrete.args.SL)

              expect_error(object = {getOutput(concrete.args.SL)})
              expect_error(object = {getOutput(concrete.est, Estimand = "a")})
              expect_error(object = {getOutput(concrete.est, Intervention = "a")})
              expect_error(object = {getOutput(concrete.est, Estimand = "RD", Intervention = 1)})
              expect_error(object = {getOutput(concrete.est, Signif = NULL)})
              expect_error(object = {getOutput(concrete.est, GComp = 5)})
          }
)

test_that("getOutput() returns ConcreteOut class with correct structure",
          code = {
              data <- as.data.table(survival::pbc)[, c("time", "status", "trt", "id", "age", "sex")]
              set.seed(123)
              data[, trt := sample(0:1, length(trt), replace = TRUE)]
              data[, status := as.numeric(status >= 1)]

              concrete.args <- formatArguments(Data = data, EventTime = "time", EventType = "status",
                                               Treatment = "trt", ID = "id", Intervention = makeITT(),
                                               TargetTime = 2500, TargetEvent = NULL,
                                               Model = NULL, Verbose = FALSE, MaxUpdateIter = 2)
              suppressWarnings(concrete.est <- doConcrete(concrete.args))

              output <- getOutput(concrete.est)

              expect_s3_class(output, "ConcreteOut")
              expect_s3_class(output, "data.table")

              # Check essential columns exist
              expect_true("Intervention" %in% colnames(output))
              expect_true("Time" %in% colnames(output))
              expect_true("Event" %in% colnames(output))
              expect_true("Estimand" %in% colnames(output))
              expect_true("Pt Est" %in% colnames(output))
          }
)

test_that("getOutput() Estimand parameter works correctly",
          code = {
              data <- as.data.table(survival::pbc)[, c("time", "status", "trt", "id", "age", "sex")]
              set.seed(123)
              data[, trt := sample(0:1, length(trt), replace = TRUE)]
              data[, status := as.numeric(status >= 1)]

              concrete.args <- formatArguments(Data = data, EventTime = "time", EventType = "status",
                                               Treatment = "trt", ID = "id", Intervention = makeITT(),
                                               TargetTime = 2500, TargetEvent = NULL,
                                               Model = NULL, Verbose = FALSE, MaxUpdateIter = 2)
              suppressWarnings(concrete.est <- doConcrete(concrete.args))

              # Test Risk estimand filter - note: stored as "Abs Risk" in output
              output_risk <- getOutput(concrete.est, Estimand = "Risk")
              expect_true(all(output_risk$Estimand == "Abs Risk"))

              # Test multiple estimands (RD, RR require 2 interventions)
              output_risk2 <- getOutput(concrete.est, Estimand = c("Risk"))
              expect_true(length(unique(output_risk2$Estimand)) >= 1)
          }
)

test_that("getOutput() GComp parameter works correctly",
          code = {
              data <- as.data.table(survival::pbc)[, c("time", "status", "trt", "id", "age", "sex")]
              set.seed(123)
              data[, trt := sample(0:1, length(trt), replace = TRUE)]
              data[, status := as.numeric(status >= 1)]

              concrete.args <- formatArguments(Data = data, EventTime = "time", EventType = "status",
                                               Treatment = "trt", ID = "id", Intervention = makeITT(),
                                               TargetTime = 2500, TargetEvent = NULL,
                                               Model = NULL, Verbose = FALSE, MaxUpdateIter = 2,
                                               GComp = TRUE)
              suppressWarnings(concrete.est <- doConcrete(concrete.args))

              # With GComp = TRUE - gcomp estimates should be included via Estimator column
              output_gcomp <- getOutput(concrete.est, GComp = TRUE)
              expect_true(attr(output_gcomp, "GComp"))
              # GComp adds rows with Estimator = "gcomp"
              expect_true("gcomp" %in% output_gcomp$Estimator)

              # With GComp = FALSE - only tmle estimates
              output_no_gcomp <- getOutput(concrete.est, GComp = FALSE)
              expect_false(attr(output_no_gcomp, "GComp"))
              expect_false("gcomp" %in% output_no_gcomp$Estimator)
          }
)
