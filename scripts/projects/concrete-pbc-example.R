


# prepare dataset -----------------------------------------------------------------------------


# data$time = event/censoring time
# data$status = event type (0 for censored)
# data$trt = intervention variable (binary point treatment)

# plot data
for (j in seq_along(unique(data$status)[unique(data$status) > 0])) {
    plot(prodlim(Hist(time, status) ~ strata(trt), data = data),
         cause = unique(data$status)[unique(data$status) > 0][j],
         add = j > 1, lty = j, confint = FALSE, legend = FALSE)
}


# specify targets -----------------------------------------------------------------------------
target.event <- sort(unique(data[status > 0, status]))
target.time <- quantile(data$time, 2:4*0.2)

intervention <- concrete:::IntentToTreat
# intervention <- list("A == 1" = list("intervention" = function(a, L) {rep_len(1, length(a))},
#                                      "g.star" = function(a, L) {as.numeric(a == 1)}),
#                      "A == 0" = list("intervention" = function(a, L) {rep_len(0, length(a))},
#                                      "g.star" = function(a, L) {as.numeric(a == 0)}))

# specify intervention(s) ------------------------------------------------------------------------
# list of arguments to be passed into origami::make_folds(), see ?origami::make_folds for details


# specify learners ----------------------------------------------------------------------------

cv.arg <- list(n = nrow(data), fold_fun = folds_vfold, cluster_ids = NULL, strata_ids = NULL)

# propensity score model
a_lrnrs <- make_learner(Lrnr_glm) # sl3 learner object: see tlverse sl3 documentation
cens_lib <- list(model.null = Surv(time, status == 0) ~ 1,
                 Surv(time, status == 0) ~ trt + age + sex)
evnt1_lib <- list(mod.main = Surv(time, status == 1) ~ trt + age + sex,
                  mod.int = Surv(time, status == 1) ~ trt:age + trt:sex)
evnt2_lib <- list(mod1 = Surv(time, status == 2) ~ trt + age + sex)

model <- list("Trt" = a_lrnrs,
              "0" = cens_lib,
              "1" = evnt1_lib,
              "2" = evnt2_lib)


# plug data and arguments in to checking function ---------------------------------------------
concrete.args <- formatArguments(DataTable = data[, c("time", "status", "trt", "id", "age", "sex")],
                                 EventTime = "time", EventType = "status",
                                 Treatment = "trt", ID = "id", Intervention = intervention,
                                 TargetTime = target.time, TargetEvent = target.event,
                                 CVArg = cv.arg, Model = model, Verbose = TRUE)


# run concrete --------------------------------------------------------------------------------
concrete.est <- doConcrete(ConcreteArgs = concrete.args)

# get estimates of target params (RR, RD, or Risks) -------------------------------------------
concrete.outputs <- getOutput(Estimate = concrete.est, Estimand = c("RR", "RD", "Risk"),
                              TargetTime = target.time, TargetEvent = target.event, GComp = TRUE)

concrete.outputs$RD
concrete.outputs$RR
attr(concrete.outputs$RR, "regime")


# plot estimates ------------------------------------------------------------------------------
ax <- lapply(concrete.outputs$Risk, function(risks) {
    risks <- setDT(risks)[Estimator == "tmle", ]
    y <- lapply(unique(risks$Event), function(j) {
        risk.j <- risks[Event == j, ]
        points(risk.j$Time, risk.j$Risk, type = "b")
    })
})


