devtools::load_all(".")
data <- as.data.table(survival::pbc)
set.seed(12345)
data[is.na(trt), trt := sample(1:2, sum(is.na(trt)), replace = TRUE)][, trt := trt - 1]

# data[, status := as.numeric(status >= 1)] # competing risk or not
intervention <- makeITT()
target.time <- 500 * (5:7)
target.event <- sort(unique(data[status > 0, status]))
# a_lrnrs <- make_learner(Stack, Lrnr_glm$new(), Lrnr_glmnet$new())
# model <- list("trt" = a_lrnrs,
#               "0" = list(mod1 = Surv(time, status == 0) ~ trt + age + sex),
#               "1" = list(mod1 = Surv(time, status == 1) ~ trt + age + sex))

concrete.args <- formatArguments(DataTable = data[, c("time", "status", "trt", "id", "age", "sex")],
                                 EventTime = "time", EventType = "status",
                                 Treatment = "trt", ID = "id", Intervention = intervention,
                                 TargetTime = target.time, TargetEvent = target.event,
                                 Model = NULL, Verbose = TRUE, ReturnModels = TRUE)
concrete.args <- formatArguments(ConcreteArgs = concrete.args)
concrete.est <- doConcrete(ConcreteArgs = concrete.args)

concrete.out<- getOutput(Estimate = concrete.est, TargetTime = target.time,
                          TargetEvent = target.event, GComp = TRUE)
str(concrete.out)
concrete.rd <- concrete.out$RD[order(Estimator, Time, Event)]
library(ggplot2)
ggplot(data = concrete.rd, aes(x = as.factor(Time), y = RD, colour = Estimator, group = Estimator)) + facet_wrap(~Event, nrow = 2) + 
    geom_errorbar(aes(ymin = RD - 1.96*se, ymax = RD + 1.96*se), width = 0.8, position = position_dodge(width=0.3)) +
    geom_point(size = 2, position = position_dodge(width=0.3)) + theme_minimal()
