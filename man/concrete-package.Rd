\name{concrete-package}
\alias{concrete-package}
\alias{concrete}
\docType{package}
\title{
  One-step continuous-time Targeted Minimum Loss-Based Estimator (TMLE) for outcome-specific absolute risk estimands in right-censored survival settings with or without competing risks
}
\description{
Implements the methodology described in Rytgaard et al. (2023) <doi:10.1111/biom.13856> and Rytgaard and van der Laan (2023) <doi:10.1007/s10985-022-09576-2>. Currently can be used to estimate the effects of static or dynamic interventions on binary treatments given at baseline, cross-validated initial estimation of treatment propensity is done using the 'SuperLearner' package, and initial estimation of conditional hazards is done using ensembles of Cox regressions from the 'survival' package or Coxnet from the 'glmnet' package.
}
\details{
formatArguments()
many check...(), format...() functions
getInitialEstimates()
getPropScores()
getHazEstimates()
getEIC()
getIC()
doTMLEUpdate()
getOutput()
}
\author{
David Chen, <david.chen49@berkeley.edu>
Maintainer: David Chen <david.chen49@berkeley.edu>
}
\references{
  Rytgaard et al. (2023) <doi:10.1111/biom.13856>  
  Rytgaard and van der Laan (2023) <doi:10.1007/s10985-022-09576-2>
}
\keyword{ package }
\seealso{
  \link[SuperLearner]{SuperLearner}  
  \link[survival]{coxph}
  \link[glmnet]{glmnet}
}
\examples{
library(concrete)
library(data.table)
set.seed(12345)
data <- as.data.table(survival::pbc)
data <- data[!is.na(trt), ][, trt := trt - 1]
data <- data[, c("time", "status", "trt", "age", "sex", "albumin")]
\donttest{
ConcreteArgs <- formatArguments(DataTable = data,
                                EventTime = "time",
                                EventType = "status",
                                Treatment = "trt",
                                Intervention = 0:1,
                                TargetTime = 1500,
                                TargetEvent = 1:2,
                                MaxUpdateIter = 250,
                                CVArg = list(V = 10),
                                Verbose = FALSE)
ConcreteEst <- doConcrete(ConcreteArgs)
ConcreteOut <- getOutput(ConcreteEst)
}

## Joint Intervention
data <- data[, trt2 := sample(0:1, .N, replace = TRUE, prob = c(0.3, .7))]
Intervention <- makeITT("A1" = data.frame(trt = rep_len(1, nrow(data)), 
                                          trt2 = rep_len(1, nrow(data))), 
                        "A0" = data.frame(trt = rep_len(0, nrow(data)), 
                                          trt2 = rep_len(0, nrow(data))))
\donttest{
ConcreteArgs <- formatArguments(DataTable = data,
                                EventTime = "time",
                                EventType = "status",
                                Treatment = c("trt", "trt2"),
                                Intervention = Intervention,
                                TargetTime = 2000,
                                TargetEvent = 1:2,
                                MaxUpdateIter = 250,
                                CVArg = list(V = 10),
                                Verbose = FALSE)
ConcreteEst <- doConcrete(ConcreteArgs)
ConcreteOut <- getOutput(ConcreteEst)
}
}
