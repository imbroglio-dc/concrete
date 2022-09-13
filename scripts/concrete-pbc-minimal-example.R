devtools::install_github("imbroglio-dc/concrete")
library(data.table)
library(concrete)
set.seed(12345)

## formatting the data ----
data <- as.data.table(survival::pbc)
data[, trt := sample(0:1, nrow(data), replace = TRUE)]
data[is.na(stage), stage := sample(1:4, sum(is.na(stage)), replace = TRUE)][, stage := factor(stage)]
# data[, status := as.numeric(status >= 1)] # to eliminate competing risks
data <- data[, c("id", "time", "status", "trt", "age", "sex", "albumin", "stage")]

## specifying the tmle problem ----
concrete_args <- formatArguments(DataTable = data,
                                 EventTime = "time", 
                                 EventType = "status",
                                 Treatment = "trt", 
                                 ID = "id", 
                                 Intervention = makeITT(),
                                 TargetTime = 500 * (2:5), 
                                 TargetEvent = 1:2,
                                 Model = NULL, 
                                 Verbose = TRUE, 
                                 ReturnModels = TRUE)

## example of modifying 'ConcreteArgs' objects & rechecking them with formatArguments() ----
concrete_args$Model$`1`$model2 <- "~ trt + sex + stage"
concrete_args <- formatArguments(concrete_args)

## running Helene's one-step continuous-time survival TMLE ----
concrete_est <- doConcrete(ConcreteArgs = concrete_args)

# getOutput needs more work, i.e. pretty print, summary, and plot methods
concrete_out <- getOutput(Estimate = concrete_est, GComp = TRUE)
concrete_rd <- concrete_out$RD
concrete_rr <- concrete_out$RR
concrete_risks <- concrete_out$Risk

# example of future plotting output
library(ggplot2)
concrete_rd %>% mutate(Time = as.factor(Time)) %>% 
ggplot(aes(x = Time, y = RD, colour = Estimator, group = Estimator)) + 
    facet_wrap(~Event, nrow = 2) + 
    geom_errorbar(aes(ymin = RD - 1.96*se, ymax = RD + 1.96*se), 
                  width = 0.8, position = position_dodge(width = 0.3)) +
    geom_point(size = 2, position = position_dodge(width = 0.3)) + 
    theme_minimal()
