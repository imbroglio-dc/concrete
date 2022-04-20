# helene's repo
try(setwd(dir = "/Shared/Projects/continuousTMLE/"),silent = TRUE)
try(setwd(dir = "~/research/SoftWare/continuousTMLE/"),silent = TRUE)
x <- lapply(paste0("R/", list.files("R/")), source)
# concrete 
try(setwd(dir = "/Shared/Projects/ConCR-TMLE/"),silent = TRUE)
try(setwd("~/research/SoftWare/devel-tmle-survival/ConCR-TMLE/"),silent = TRUE)
source("scripts/packages.R")
x <- lapply(list.files("R/*",path = "R",full.names = TRUE), source)
source("scripts/prepare-pbc.R")
# obsolete concrete with hard coded variable names



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
                                      ## MaxUpdateIter = MaxUpdateIter,
                                      MaxUpdateIter = 2,
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
