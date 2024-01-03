library(concrete)
data <- survival::pbc[, c("time", "status", "trt", "age", "sex", "albumin")]
data <- subset(data, subset = !is.na(data$trt))
data$trt <- data$trt - 1

# Specify Analysis
ConcreteArgs <- formatArguments(
  DataTable = data, 
  EventTime = "time",
  EventType = "status", 
  Treatment = "trt",
  Intervention = 0:1, 
  TargetTime = 365.25/2 * (6:12), 
  TargetEvent = 1:2, 
  CVArg = list(V = 10), 
  Verbose = FALSE
)

# Compute
ConcreteEst <- doConcrete(ConcreteArgs)

# Return Output
ConcreteOut <- getOutput(ConcreteEst, Estimand = "RD", Simultaneous = TRUE)
plot(ConcreteOut, NullLine = TRUE, ask = FALSE)

ggplot2::ggsave(filename = "PBC-RD.png")


ConcreteArgs <- formatArguments(
  DataTable = data,               # data.frame or data.table
  EventTime = "time",             # name of event time variable
  EventType = "status",           # name of event status variable
  Treatment = "trt",              # name of treatment variable
  ID = NULL,                      # name of the ID variable if present in input data 
  Intervention = 0:1,             # 2 static interventions
  TargetTime = 365.25/2 * (6:12), # 7 target times: 3-6 years biannually
  TargetEvent = 1:2,              # 2 competing risks
  CVArg = list(V = 10),           # 10-Fold Cross-Validation
  Model = NULL,                   # using default Super Learner libraries
  Verbose = FALSE                 # less verbose warnings and progress messages
)

attr(ConcreteArgs[["DataTable"]], "CovNames")

TreatOver60 <- list(
  "intervention" = function(ObservedTrt, Covariates, PropScore) {
    TrtNames <- colnames(ObservedTrt)
    Intervened <- data.table::copy(ObservedTrt)
    Intervened[, (TrtNames) := lapply(.SD, function(a) {
      as.numeric(Covariates[["age"]] > 60)
    }), .SDcols = TrtNames]
    return(Intervened)
  }, 
  # g.star function copied from from makeITT()
  "g.star" = function(Treatment, Covariates, PropScore, Intervened) {
    Probability <- data.table::data.table(1 * sapply(1:nrow(Treatment), function(i) 
      all(Treatment[i, ] == Intervened[i, ])))
    return(Probability)
  }
)
TreatUnder60 <- list(
  "intervention" = function(ObservedTrt, Covariates, PropScore) {
    TrtNames <- colnames(ObservedTrt)
    Intervened <- data.table::copy(ObservedTrt)
    Intervened[, (TrtNames) := lapply(.SD, function(a) {
      as.numeric(Covariates[["age"]] <= 60)
    }), .SDcols = TrtNames]
    return(Intervened)
  }
  # if g.star function left empty, one will be copied from from makeITT()
)

ConcreteArgs <- formatArguments(
  DataTable = data, 
  EventTime = "time", 
  EventType = "status", 
  Treatment = "trt",
  Intervention = list("Treat>60" = TreatOver60, 
                      "Treat<=60" = TreatUnder60), 
  TargetTime = 365.25/2 * (6:12), 
  TargetEvent = 1:2, 
  CVArg = list(V = 10), 
  RenameCovs = FALSE, ## turn off covariate pre-processing ##
  Verbose = FALSE
)

# specify regression models
ConcreteArgs$Model <- list(
  "trt" = c("SL.glmnet", "SL.bayesglm", "SL.xgboost", "SL.glm", "SL.ranger"),
  "0" = NULL, # will use the default library
  "1" = list(Surv(time, status == 1) ~ trt, Surv(time, status == 1) ~ .),
  "2" = list("Surv(time, status == 2) ~ trt", "Surv(time, status == 2) ~ .")
)



# decrease the maximum tmle update number to 50
ConcreteArgs$MaxUpdateIter <- 50

# add a candidate regression with treatment interactions
ConcreteArgs[["Model"]][["2"]][[3]] <- "Surv(time, status == 2) ~ trt*."

# validate new estimation specification
ConcreteArgs <- formatArguments(ConcreteArgs)


print(ConcreteArgs, Verbose = FALSE)


ConcreteEst <- doConcrete(ConcreteArgs)


print(ConcreteEst, Verbose = FALSE)

plot(ConcreteEst)
ggplot2::ggsave("ConcreteEst-plot.png", width = 5, height = 3, units = "in")

ConcreteOut <- getOutput(ConcreteEst = ConcreteEst, 
                         Estimand = "RD",
                         Intervention = 1:2, 
                         GComp = TRUE, 
                         Simultaneous = TRUE, 
                         Signif = 0.05)
head(ConcreteOut, 12)
