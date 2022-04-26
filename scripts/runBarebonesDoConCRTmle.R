# helene's repo
contmle.dir <- c("/Shared/Projects/continuousTMLE/R", "~/research/SoftWare/continuousTMLE/R")
x <- lapply(contmle.dir, function(dir) lapply(list.files(dir, full.names = TRUE), source))

# obsolete concrete with hard coded variable names
obs.concrete.dir <- c("/Shared/Projects/ConCR-TMLE/obsolete-concrete/R",
                      "~/research/SoftWare/devel-tmle-survival/ConCR-TMLE/obsolete-concrete/R")
x <- lapply(obs.concrete.dir, function(dir) lapply(list.files(dir, full.names = TRUE), source))

# concrete
try(setwd(dir = "/Shared/Projects/ConCR-TMLE/"), silent = TRUE)
try(setwd("~/research/SoftWare/devel-tmle-survival/ConCR-TMLE/"), silent = TRUE)
x <- lapply(list.files("R/*",path = "R",full.names = TRUE), source)
source("scripts/packages.R")
source("scripts/prepare-pbc.R")
intervention <- list("A == 1" = list("intervention" = function(a, L) {rep_len(1, length(a))},
                                     "g.star" = function(a, L) {as.numeric(a == 1)}),
                     "A == 0" = list("intervention" = function(a, L) {rep_len(0, length(a))},
                                     "g.star" = function(a, L) {as.numeric(a == 0)}))
target.time <- 500 * (2:4)
target.event <- sort(unique(data[status > 0, status]))
a_lrnrs <- make_learner(Lrnr_glm)
model <- list("Trt" = a_lrnrs,
              "0" = list(mod1 = Surv(time, status == 0) ~ trt + age + sex),
              "1" = list(mod1 = Surv(time, status == 1) ~ trt + age + sex))

# plot(prodlim(Hist(time,status)~trt,data = data))

concrete.args <- formatArguments(DataTable = data[, c("time", "status", "trt", "id", "age", "sex")],
                                 EventTime = "time", EventType = "status",
                                 Treatment = "trt", ID = "id", Intervention = intervention,
                                 TargetTime = target.time, TargetEvent = target.event,
                                 Model = model, Verbose = TRUE)
output <- doConcrete(ConcreteArgs = concrete.args)

concrete.ate <- lapply(output, function(out.a) {
    do.call(rbind, lapply(sort(unique(data[status > 0, status])), function(j) {
        risks <- apply(out.a[["Hazards"]][[as.character(j)]] * out.a[["EvntFreeSurv"]], 2, cumsum)
        Psi <- cbind("tau" = target.time, "tmle.est" = rowMeans(risks[attr(output, "times") %in% target.time, ]))
        tmle.se <- subset(out.a$SummEIC[Event == j, ], select = c("Time","seEIC"))
        Psi <- merge(Psi, tmle.se[, list("tau" = Time, "tmle.se" = seEIC / sqrt(ncol(risks)))], by = "tau")
        return(cbind("J" = j, Psi))
    }))
})

# ate
concrete.ate <- as.data.table(merge(concrete.ate[[1]], concrete.ate[[2]], by = c("J", "tau")))[order(tau)]
concrete.ate[, list(J = J, tau = tau, tmle.est = tmle.est.x - tmle.est.y, tmle.se = sqrt(tmle.se.x^2 + tmle.se.y^2))]

#=================================================================================================================
# obsolete
#=================================================================================================================
# Intervention <- list(
#     "A == 1" = function(a, L) {
#         regime <- rep_len(1, length(a))
#         attr(regime, "g.star") <- function(a) {as.numeric(a == 1)}
#         return(regime)
#     },
#     "A == 0" = function(a, L) {
#         regime <- rep_len(0, length(a))
#         attr(regime, "g.star") <- function(a) {as.numeric(a == 0)}
#         return(regime)
#     })
# Models <- list("Trt" = a_lrnrs,
#                "0" = list(mod1 = Surv(Time, Event == 0) ~ Trt + age + sex),
#                "1" = list(mod1 = Surv(Time, Event == 1) ~ Trt + age + sex))
# OBSoutput <- OBSdoConCRTmle(EventTime = data$time,
#                             EventType = data$status,
#                             Treatment = data$trt,
#                             Intervention = Intervention,
#                             CovDataTable = data[, c("age", "sex")],
#                             ID = data$id,
#                             TargetTimes = 500*2:4,
#                             TargetEvents = sort(unique(data[status > 0, status])),
#                             Models = Models,
#                             CVArgs = NULL,
#                             NumUpdateSteps = 25,
#                             OneStepEps = 0.1,
#                             MinNuisance = 0.05,
#                             PropScoreBackend = "sl3",
#                             Verbose = TRUE,
#                             GComp = TRUE)
#




#=================================================================================================================
# contmle
#=================================================================================================================
estimation <- list("cause1" = list(fit = "cox",
                                   model = Surv(time, status == 1) ~ trt + age + sex),
                   "cens" = list(fit = "cox",
                                 model = Surv(time, status == 0) ~ trt + age + sex))

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
