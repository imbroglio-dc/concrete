library(survival); library(zoo); library(prodlim); library(nleqslv)
# library(foreach); library(doRNG); library(doParallel)
library(tidyverse); library(data.table)
library(sl3); library(origami)
try(setwd(dir = "/Shared/Projects/continuousTMLE/"))
lapply(paste0("R/", list.files("R/")), source)
try(setwd(dir = "/Shared/Projects/ConCR-TMLE/"))
lapply(paste0("R/", list.files("R/")), source)

data <- as.data.table(survival::pbc)
set.seed(12345)
data[is.na(trt), trt := sample(1:2, sum(is.na(trt)), replace = TRUE)][, trt := trt - 1]
data[, status := as.numeric(status >= 1)]

intervention <- list("A == 1" = list("intervention" = function(a, L) {rep_len(1, length(a))},
                                     "g.star" = function(a, L) {as.numeric(a == 1)}),
                     "A == 0" = list("intervention" = function(a, L) {rep_len(0, length(a))},
                                     "g.star" = function(a, L) {as.numeric(a == 0)}))
target.time <- 500 * (2:4)
target.event <- sort(unique(data[status > 0, status]))

logreg <- make_learner(Lrnr_glm)
lasso <- make_learner(Lrnr_glmnet) # alpha default is 1
ridge <- Lrnr_glmnet$new(alpha = 0)
e_net <- make_learner(Lrnr_glmnet, alpha = 0.5)
a_lrnrs <- logreg # make_learner(Stack, logreg, lasso, ridge, e_net)

model <- list("Trt" = a_lrnrs,
              "0" = list(mod1 = Surv(time, status == 0) ~ trt + age + sex),
              "1" = list(mod1 = Surv(time, status == 1) ~ trt + age + sex))

estimation <- list("cause1" = list(fit = "cox",
                                   model = Surv(time, status == 1) ~ trt + age + sex),
                   "cens" = list(fit = "cox",
                                 model = Surv(time, status == 0) ~ trt + age + sex))
concreteArgs <- formatArguments(DataTable = data[, c("time", "status", "trt", "id", "age", "sex")],
                                EventTime = "time", EventType = "status",
                                Treatment = "trt", ID = "id", Intervention = intervention,
                                TargetTime = target.time, TargetEvent = target.event,
                                Model = model, Verbose = TRUE)
output <- with(concreteArgs, doConCRTmle(Data = Data,
                                      EventTime = EventTime,
                                      EventType = EventType,
                                      Treatment = Treatment,
                                      CovDataTable = CovDataTable,
                                      LongTime = LongTime,
                                      ID = ID,
                                      Events = Events,
                                      Censored = Censored,
                                      TargetTime = TargetTime,
                                      TargetEvent = TargetEvent,
                                      Regime = Regime,
                                      CVArg = CVArg,
                                      Model = Model,
                                      PropScoreBackend = PropScoreBackend,
                                      MaxUpdateIter = MaxUpdateIter,
                                      OneStepEps = OneStepEps,
                                      MinNuisance = MinNuisance,
                                      Verbose = Verbose,
                                      GComp = GComp))

tmp <- lapply(output, function(out.a) {
    do.call(rbind, lapply(sort(unique(data[status > 0, status])), function(j) {
        risks <- apply(out.a[["Hazards"]][[as.character(j)]] * out.a[["EvntFreeSurv"]], 2, cumsum)
        Psi <- cbind("tau" = target.time, "tmle.est" = rowMeans(risks[attr(output, "times") %in% target.time, ]))
        tmle.se <- subset(out.a$SummEIC[Event == j, ], select = c("Time","seEIC"))
        Psi <- merge(Psi, tmle.se[, list("tau" = Time, "tmle.se" = seEIC / sqrt(ncol(risks)))], by = "tau")
        return(cbind("J" = j, Psi))
    }))
})

# ate
tmp <- as.data.table(merge(tmp[[1]], tmp[[2]], by = c("J", "tau")))[order(tau)]
tmp[, list(J = J, tau = tau, tmle.est = tmle.est.x - tmle.est.y, tmle.se = sqrt(tmle.se.x^2 + tmle.se.y^2))]


contmle_output <- contmle(dt = data,
                          target = target.event,
                          iterative = FALSE,
                          treat.effect = "ate",
                          tau = target.time,
                          estimation = estimation,
                          treat.model = trt ~ sex + age,
                          sl.models = list(mod1 = list(Surv(time, status == 1) ~ trt + age + sex))
)
contmle_output$tmle
