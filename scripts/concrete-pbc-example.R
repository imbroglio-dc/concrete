
# load libraries ------------------------------------------------------------------------------

library(survival)
library(data.table)
library(zoo)
library(prodlim)
library(nleqslv)
library(sl3)
library(origami)
# devtools::install_github("imbroglio-dc/concrete")
library(concrete)

# prepare dataset -----------------------------------------------------------------------------

data <- as.data.table(survival::pbc)
set.seed(12345)
data[is.na(trt), trt := sample(1:2, sum(is.na(trt)), replace = TRUE)][, trt := trt - 1]
# data[, status := as.numeric(status >= 1)]

# data$time = event/censoring time
# data$status = event type (0 for censored)
# data$trt = intervention variable (binary point treatment)

for (j in seq_along(unique(data$status)[unique(data$status) > 0])) {
    plot(prodlim(Hist(time, status) ~ strata(trt), data = data),
         cause = unique(data$status)[unique(data$status) > 0][j], add = j > 1, lty = j, legend = FALSE)
}


# specify targets -----------------------------------------------------------------------------

target.event <- sort(unique(data[status > 0, status]))
target.time <- quantile(data$time, 1:3*0.25)

intervention <- list("A == 1" = list("intervention" = function(a, L) {rep_len(1, length(a))},
                                     "g.star" = function(a, L) {as.numeric(a == 1)}),
                     "A == 0" = list("intervention" = function(a, L) {rep_len(0, length(a))},
                                     "g.star" = function(a, L) {as.numeric(a == 0)}))
# concrete:::IntentToTreat # template intervention for binary intervention

a_lrnrs <- make_learner(Lrnr_glm) # sl3 learner object: see tlverse sl3 documentation
cv.arg <- list(n = nrow(data), fold_fun = folds_vfold, # list of arguments to be passed into origami::make_folds()
               cluster_ids = NULL, strata_ids = NULL) # see ?origami::make_folds for details

model <- list("Trt" = a_lrnrs,
              "0" = list(mod1 = Surv(time, status == 0) ~ trt + age + sex),
              "1" = list(mod1 = Surv(time, status == 1) ~ trt + age + sex),
              "2" = list(mod1 = Surv(time, status == 2) ~ trt + age + sex))

concrete.args <- formatArguments(DataTable = data[, c("time", "status", "trt", "id", "age", "sex")],
                                 EventTime = "time", EventType = "status",
                                 Treatment = "trt", ID = "id", Intervention = intervention,
                                 TargetTime = target.time, TargetEvent = target.event,
                                 CVArg = cv.arg, Model = model, Verbose = TRUE, )
concrete.est <- doConcrete(ConcreteArgs = concrete.args)

concrete.outputs <- getOutput(Estimate = concrete.est, Estimand = c("RR", "RD", "Risk"), TargetTime = target.time,
                              TargetEvent = target.event, GComp = TRUE)

concrete.outputs$RD
concrete.outputs$RR
attr(concrete.outputs$RR, "regime")

x <- lapply(concrete.outputs$Risk, function(risks) {
    risks <- setDT(risks)[Estimator == "tmle", ]
    y <- lapply(unique(risks$Event), function(j) {
        risk.j <- risks[Event == j, ]
        lines(risk.j$Time, risk.j$Risk, type = "b")
    })
})

lines(x = concrete.outputs$Risk$`A == 1`$Time, y = concrete.outputs$Risk$`A == 1`$Risk)

