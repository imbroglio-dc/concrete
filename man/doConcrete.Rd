% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/doConCRTmle.R
\name{doConcrete}
\alias{doConcrete}
\title{doConcrete}
\usage{
doConcrete(ConcreteArgs)
}
\arguments{
\item{ConcreteArgs}{: output from concrete::formatArguments}
}
\value{
tbd
}
\description{
doConcrete
}
\examples{
library(data.table)
library(survival)
library(concrete)
data <- as.data.table(survival::pbc)
data[, trt := sample(0:1, nrow(data), TRUE)]
cols <- c("id", "time", "status", "trt",
          "age", "albumin", "sex", "bili")
data <- data[, .SD, .SDcols = cols]

intervention <- concrete:::ITT
target.time <- 2500
target.event <- 1:2
model <- list("trt" = c("SL.glm", "SL.glmnet"),
              "0" = list(Surv(time, status == 0) ~ .),
              "1" = list(Surv(time, status == 1) ~ .),
              "2" = list(Surv(time, status == 2) ~ .))

# formatArguments() returns correctly formatted arguments for doConcrete()
concrete.args <- formatArguments(DataTable = data,
                                 EventTime = "time",
                                 EventType = "status",
                                 Treatment = "trt",
                                 ID = "id",
                                 Intervention = intervention,
                                 TargetTime = target.time,
                                 TargetEvent = target.event,
                                 Model = model, Verbose = FALSE)

# doConcrete() returns tmle (and g-comp plug-in) estimates of targeted risks
concrete.est <- doConcrete(concrete.args)

# getOutput returns risk difference, relative risk, and treatment-specific risks
concrete.out <- getOutput(Estimate = concrete.est, Estimand = c("rd", "rr", "risk"), 
                          TargetTime = target.time, TargetEvent = target.event, GComp = TRUE)
concrete.out$RD
concrete.out$RR
concrete.out$Risk
}
