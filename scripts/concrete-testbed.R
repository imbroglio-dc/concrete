devtools::load_all(".")
library(data.table)
data <- as.data.table(survival::pbc)
set.seed(12345)
data[is.na(trt), trt := sample(1:2, sum(is.na(trt)), replace = TRUE)][, trt := trt - 1]

# data[, status := as.numeric(status >= 1)] # competing risk or not
intervention <- list("A == 1" = list("intervention" = function(a, L) {rep_len(1, length(a))},
                                     "g.star" = function(a, L) {as.numeric(a == 1)}),
                     "A == 0" = list("intervention" = function(a, L) {rep_len(0, length(a))},
                                     "g.star" = function(a, L) {as.numeric(a == 0)}))
target.time <- 500 * (2:4)
target.event <- sort(unique(data[status > 0, status]))
# a_lrnrs <- make_learner(Stack, Lrnr_glm$new(), Lrnr_glmnet$new())
# model <- list("trt" = a_lrnrs,
#               "0" = list(mod1 = Surv(time, status == 0) ~ trt + age + sex),
#               "1" = list(mod1 = Surv(time, status == 1) ~ trt + age + sex))

concrete.args <- formatArguments(DataTable = data[, c("time", "status", "trt", "id", "age", "sex")],
                                 EventTime = "time", EventType = "status",
                                 Treatment = "trt", ID = "id", Intervention = intervention,
                                 TargetTime = target.time, TargetEvent = target.event,
                                 Model = NULL, PropScoreBackend = "SuperLearner", Verbose = TRUE)
formatArguments(ConcreteArgs = concrete.args)
concrete.est <- doConcrete(ConcreteArgs = concrete.args)

concrete.ate <- getOutput(Estimate = concrete.est, Estimand = c("rd"), TargetTime = target.time,
                          TargetEvent = target.event, GComp = TRUE)
concrete.ate$RD[order(Time, Event, Estimator)]