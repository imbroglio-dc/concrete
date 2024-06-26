loadPackages <- function() {
  library(survtmle)
  # library(MOSS); library(SuperLearner)
  # surv_tmle.dir <- c("/Shared/Projects/Roadmap_CVOT/R/functions/")
  # x <- lapply(surv_tmle.dir, function(dir) lapply(list.files(dir, full.names = TRUE),
  #                                                 function(x) try(source(x), silent = TRUE)))
  # unloadNamespace("mice")
  
  library(survival); library(zoo)
  contmle.dir <- c("/Shared/Projects/continuousTMLE/R", "~/research/SoftWare/continuousTMLE/R")
  x <- lapply(contmle.dir, function(dir) lapply(list.files(dir, full.names = TRUE), source))
  
  library(data.table); library(concrete)
  concrete.dir <- c("/Shared/Projects/concrete/")
  x <- lapply(concrete.dir, function(dir) try(source(paste0(dir, "scripts/sim-data/sim_functions.R"))))
  library(tidyverse)
  invisible(NULL)
}
loadPackages()
data.table::setDTthreads(threads = parallel::detectCores(), restore_after_fork = TRUE)
concrete.dir <- "/Shared/Projects/concrete/"
TargetTime <- seq(730, 1460, length.out = 5)

# Simulation Params
B <- 1e3
OutputPath <- paste0(concrete.dir, "scripts/sim-data/sim-out/")

# true risks --------------------------------------------------------------

if (file.exists(paste0(concrete.dir, "scripts/sim-data/TrueRisks.csv"))) {
  risks <- read.csv(paste0(concrete.dir, "scripts/sim-data/TrueRisks.csv"))[, -1]
} else {
  risks1 <- cbind("Intervention" = "A=1",
                  getTrueRisks(n = 5e6, target_time = sort(union(1:2000, TargetTime)), 
                               assign_A = function(W, n) return(rep_len(1, n)), 
                               # ltfu_coefs = c(7.5e-5, 1, 4, .3),
                               parallel = FALSE, rep = 1)); gc()
  risks0 <- cbind("Intervention" = "A=0",
                  getTrueRisks(n = 5e6, target_time = sort(union(1:2000, TargetTime)), 
                               assign_A = function(W, n) return(rep_len(0, n)), 
                               # ltfu_coefs = c(7.5e-5, 1, 4, .3),
                               parallel = FALSE, rep = 1)); gc()
  risks <- melt(rbind(risks1, risks0), id.vars = c("Intervention", "Time"), 
                variable.name = "Event", value.name = "True")
  write.csv(risks, paste0(concrete.dir, "scripts/sim-data/TrueRisks.csv"))
  rm(list = c("risks0", "risks1")); gc()
}
risks %>% mutate(Event = factor(Event)) %>%
  ggplot(aes(x = Time, y = True, colour = Event, linetype = Intervention)) + 
  geom_line(linewidth = .8) + theme_minimal() + labs(title = "True Competing Risks")


# simulation --------------------------------------------------------------

if (length(list.files(OutputPath, full.names = TRUE)) < B) {
  library(foreach)
  library(doParallel)
  n_cores <- min(parallel::detectCores(), 12)
  j <- 0
  set.seed(0)
  seeds <- sample(0:12345678, size = B)
  Start <- Sys.time()
  while(j < B) {
    cl <- makeCluster(n_cores, type = "FORK")
    registerDoParallel(cl)
    sim.result <- try({
      foreach(i = seq_along(seeds),
              .errorhandling = "remove",
              .packages = "readr") %dopar% {
                loadPackages()
                DiagPath <- paste0(concrete.dir, "scripts/sim-data/diagnostic.txt")
                estimators <- c(
                  "concrete",
                  # "contmle",
                  "survtmle-6mo",
                  # "survtmle-4mo",
                  "survtmle-3mo", 
                  # "survtmle-2mo", 
                  "survtmle-1mo",
                  "Aalen-Johansen")
                
                discretizations <- c("survtmle-6mo" = 2, "survtmle-4mo" = 3,
                                     "survtmle-3mo" = 4, "survtmle-2mo" = 6,
                                     "survtmle-1mo" = 12)
                
                if (file.exists(paste0(OutputPath, i, ".RDS"))) {
                  out <- read_rds(paste0(OutputPath, i, ".RDS"))
                  out$meta$time <- left_join(data.table("Estimator" = estimators), 
                                             out$meta$time)
                } else {
                  out <- list()
                  out$estimates <- data.table()
                  out$meta <- list("time" = data.table("Estimator" = estimators, 
                                                       "Time" = NaN),
                                   "error" = list(), 
                                   "notconverge" = list())
                }
                
                MaxIter <- 500
                Data <- simConCR(n = 2.5e3, random_seed = seeds[i])
                
                # helper functions --------------------------------------------------------
                simConcrete <- function(TargetTime, Data, MaxIter, DiagPath = NULL, j = 1) {
                  require(concrete); require(data.table)
                  if (!is.null(DiagPath)) {
                    readr::write_lines(paste0("Run ", j + i - 1, " ; concrete ; start"),
                                       file = DiagPath, append = TRUE)
                  }
                  concreteArgs <- formatArguments(Data = Data,
                                                  EventTime = "TIME",
                                                  EventType = "EVENT",
                                                  Treatment = "ARM",
                                                  ID = "id",
                                                  Intervention = 0:1,
                                                  TargetTime = TargetTime,
                                                  MaxUpdateIter = MaxIter,
                                                  Verbose = FALSE)
                  concreteArgs$Model$ARM <- c("SL.glm", "SL.glmnet")
                  for (VarName in names(concreteArgs$Model)) {
                    if (VarName %in% unique(Data[["EVENT"]])) {
                      concreteArgs[["Model"]][[as.character(VarName)]][["Coxnet"]] <- "coxnet"
                    }
                  }
                  concreteArgs <- formatArguments(concreteArgs)
                  Start <- Sys.time()
                  concreteEst <- doConcrete(concreteArgs)
                  concreteOut <- getOutput(concreteEst, Estimand = "Risk")
                  Stop <- Sys.time()
                  estimates <- cbind(Package = "concrete", concreteOut)
                  attr(estimates, "converged") <- attr(concreteEst, "TmleConverged")$converged
                  attr(estimates, "time") <- difftime(Stop, Start, units = "mins")
                  if (!is.null(DiagPath)) {
                    readr::write_lines(paste0("Run ", j + i - 1, " ; concrete ; ",
                                              format(unclass(attr(estimates, "time")), digits = 3),
                                              " ", attr(attr(estimates, "time"), "units")),
                                       file = DiagPath, append = TRUE)
                  }
                  return(estimates)
                }
                
                simContmle <- function(TargetTime, Data, MaxIter, DiagPath = NULL, j = 1) {
                  required <- c("data.table", "zoo", "glmnet", "survival", "stringr",
                                "nleqslv","prodlim", "Matrix", "coefplot", "hdnom")
                  z <- sapply(required, function(p) try(library(p, character.only = TRUE),
                                                        silent = TRUE))
                  if (!is.null(DiagPath)) {
                    readr::write_lines(paste0("Run ", j + i - 1, " ; contmle ; start"),
                                       file = DiagPath, append = TRUE)
                  }
                  Start <- Sys.time()
                  contmleRun <- contmle(
                    Data, #-- Data
                    target = sort(setdiff(unique(Data$EVENT), 0)), #-- target events
                    iterative = FALSE, #-- use one-step tmle to target F1 and F2 simultaneously
                    treat.effect = "ate", #-- target the ate directly
                    tau = TargetTime, #-- time-point of interest
                    estimation = list(
                      "cens" = list(fit = "cox",
                                    model = Surv(TIME, EVENT == 0) ~ ARM+GEOGR1+SEX+AGE+MIFL+
                                      EGFMDRBC+CREATBL+BMIBL+HBA1CBL+ETHNIC+SMOKER+STROKSFL),
                      "cause1" = list(fit = "cox",
                                      model = Surv(TIME, EVENT == 1) ~ ARM+GEOGR1+SEX+AGE+MIFL+
                                        EGFMDRBC+CREATBL+BMIBL+HBA1CBL+ETHNIC+SMOKER+STROKSFL),
                      "cause2" = list(fit = "cox",
                                      model = Surv(TIME, EVENT == 2) ~ ARM+GEOGR1+SEX+AGE+MIFL+
                                        EGFMDRBC+CREATBL+BMIBL+HBA1CBL+ETHNIC+SMOKER+STROKSFL) #,
                      # "cause3" = list(fit = "cox",
                      #                 model = Surv(TIME, EVENT == 3) ~ ARM+GEOGR1+SEX+AGE+MIFL+
                      #                     EGFMDRBC+CREATBL+BMIBL+HBA1CBL+ETHNIC+SMOKER+STROKSFL)
                    ),
                    treat.model = ARM ~ GEOGR1+SEX+AGE+MIFL+EGFMDRBC+CREATBL+BMIBL+
                      HBA1CBL+ETHNIC+SMOKER+STROKSFL,
                    verbose = FALSE,
                    sl.models = NULL,
                    no.small.steps = MaxIter
                  )
                  Stop <- Sys.time()
                  attr(contmleRun, "time") <- difftime(Stop, Start, units = "mins")
                  if (!is.null(DiagPath)) {
                    readr::write_lines(paste0("Run ", j + i - 1, " ; contmle ; ",
                                              format(unclass(attr(contmleRun, "time")), digits = 3),
                                              " ", attr(attr(contmleRun, "time"), "units")),
                                       file = DiagPath, append = TRUE)
                  }
                  return(contmleRun)
                }
                
                formatContmle <- function(contmleOutput) {
                  tmleOutput <- data.table(
                    "J" = rep(names(contmleOutput$tmle), each = 2),
                    "val" = c("ATE", "se"),
                    do.call(rbind, contmleOutput$tmle)
                  )
                  tmleOutput <- melt(tmleOutput, id.vars = c("J", "val"))
                  setnames(tmleOutput, "variable", "time")
                  tmleOutput[, J := gsub("F", "", J)]
                  tmleOutput[, time := as.numeric(gsub("tau=", "", time))]
                  tmleOutput <- dcast(tmleOutput, J + time ~ val, value.var = "value")
                  setnames(tmleOutput, c("J", "time", 'ATE'), c("Event", "Time", "RD"))
                  tmleOutput <- cbind("Package" = "contmle", "Estimator" = "tmle", tmleOutput)
                  
                  gcompOutput <-
                    data.table(
                      "J" = rep(names(contmleOutput$init), each = 2),
                      "val" = c("ATE", "se"),
                      do.call(rbind, contmleOutput$init)
                    )
                  gcompOutput <- melt(gcompOutput, id.vars = c("J", "val"))
                  setnames(gcompOutput, "variable", "time")
                  gcompOutput[, J := gsub("F", "", J)]
                  gcompOutput[, time := as.numeric(gsub("tau=", "", time))]
                  gcompOutput <- dcast(gcompOutput, J + time ~ val, value.var = "value")
                  setnames(gcompOutput, c("J", "time", 'ATE'), c("Event", "Time", "RD"))
                  gcompOutput <- cbind("Package" = "contmle", "Estimator" = "gcomp", gcompOutput)
                  
                  return(as.data.frame(rbind(tmleOutput, gcompOutput)))
                }
                
                simSurvtmle <- function(TargetTime, Bins, Data, MaxIter, DiagPath = NULL, j = 1) {
                  require(SuperLearner); require(survtmle)
                  survMonths <- paste0("survtmle", as.integer(12 / Bins), "mo")
                  if (!is.null(DiagPath)) {
                    readr::write_lines(paste0("Run ", j + i - 1, " ; ", survMonths," ; start"),
                                       file = DiagPath, append = TRUE)
                  }
                  sl_lib_failure <- c("SL.glmnet", "SL.glm")
                  sl_lib_censor <- c("SL.glmnet", "SL.xgboost")
                  sl_lib_g <- c("SL.glmnet", "SL.xgboost")
                  TargetTime <- ceiling(TargetTime / 365 * Bins)
                  TargetEvent <- sort(setdiff(unique(Data$EVENT), 0))
                  W <- as.data.frame(Data[, !c("id", "TIME", "EVENT", "ARM")])
                  
                  Start <- Sys.time()
                  tmleFit <- survtmle(ftime = ceiling(Data$TIME / 365 * Bins),
                                      ftype = Data$EVENT,
                                      trt = Data$ARM,
                                      adjustVars = W,
                                      t0 = max(TargetTime),
                                      SL.ftime = sl_lib_failure,
                                      SL.ctime = sl_lib_censor,
                                      SL.trt = sl_lib_g,
                                      returnIC = TRUE,
                                      returnModels = TRUE,
                                      method = "hazard",
                                      verbose = FALSE,
                                      maxIter = 100,
                                      Gcomp = FALSE)
                  survtmleOut <- getSurvtmleTbl(survtmleOut = print(timepoints(object = tmleFit, times = TargetTime)), 
                                                TargetEvent = TargetEvent, TargetTime = TargetTime, Bins = Bins)
                  
                  Stop <- Sys.time()
                  attr(survtmleOut, "time") <- difftime(Stop, Start, units = "mins")
                  if (!is.null(DiagPath)) {
                    readr::write_lines(paste0("Run ", j + i - 1, " ; ", survMonths," ; ",
                                              format(unclass(attr(survtmleOut, "time")), digits = 3), " ",
                                              attr(attr(survtmleOut, "time"), "units")),
                                       file = DiagPath, append = TRUE)
                  }
                  return(survtmleOut)
                }
                
                getSurvtmleTbl <- function(survtmleOut, TargetEvent, TargetTime, Bins) {
                  require(dplyr)
                  x <- c("Risk", "se")
                  survtmle.list <- lapply(seq_along(x), function(i) {
                    survtmle.est <- cbind("Estimator" = paste0("survtmle-", 12 / Bins, "mo"),
                                          "Intervention" = rep(c("A=0", "A=1"),
                                                               times = length(TargetEvent)),
                                          "Event" = rep(TargetEvent, each = 2),
                                          survtmleOut[[i]])
                    survtmle.est <- melt(setDT(survtmle.est),
                                         id.vars = c("Estimator", "Intervention", "Event"),
                                         value.name = x[i], variable.name = "Time")
                    survtmle.est[, Time := rep(TargetTime, each = length(TargetEvent) * 2) / Bins * 365]
                    return(survtmle.est)
                  })
                  out <- do.call(dplyr::full_join, survtmle.list)
                  setDT(out)
                  setnames(out, "Risk", "Pt Est")
                }
                
                simAJ <- function(TargetTime, Data, DiagPath = NULL, j = 1) {
                  require(survival); require(tidyverse); requireNamespace("stringr")
                  if (!is.null(DiagPath)) {
                    readr::write_lines(paste0("Run ", j + i - 1, " ; Aalen-Johansen ; start"),
                                       file = DiagPath, append = TRUE)
                  }
                  TargetEvent <- sort(setdiff(unique(Data$EVENT), 0))
                  
                  Start <- Sys.time()
                  
                  AJ <- survfit(Surv(time = TIME, event = as.factor(EVENT)) ~ ARM, 
                                data = Data, ctype = 1)
                  AJ <- cbind(as.data.table(summary(AJ, times = TargetTime)$pstate[, -1]),
                              as.data.table(summary(AJ, times = TargetTime)$std.err[, -1]))
                  setnames(AJ, c(paste0("Risk", TargetEvent), paste0("se", TargetEvent)))
                  AJ <- cbind("Time" = TargetTime, 
                              "Intervention" = rep(0:1, each = length(TargetTime)), 
                              AJ) %>% 
                    melt(., id.vars = c("Time", "Intervention")) %>% 
                    .[, "Event" := stringr::str_extract(string = variable, pattern = "[0-9]+")] %>%
                    .[, variable := stringr::str_extract(string = variable, pattern = "[^(0-9)]+")] %>%
                    dcast(., Event + Intervention + Time ~ variable) %>% 
                    cbind(., "Package" = "survival", "Estimator" = "Aalen-Johansen", 
                          "Estimand" = "Abs Risk", "CI Low" = NA, "CI Hi" = NA, 
                          "SimCI Low" = NA, "SimCI Hi" = NA) %>% 
                    .[, Intervention := paste0("A=", Intervention)]
                  setnames(AJ, "Risk", "Pt Est")
                  
                  Stop <- Sys.time()
                  attr(AJ, "time") <- difftime(Stop, Start, units = "mins")
                  if (!is.null(DiagPath)) {
                    readr::write_lines(paste0("Run ", j + i - 1, " ; Aalen-Johansen ; ",
                                              format(unclass(attr(AJ, "time")), digits = 3), " ",
                                              attr(attr(AJ, "time"), "units")),
                                       file = DiagPath, append = TRUE)
                  }
                  return(AJ)
                } 
                
                # concrete --------------------------------------------------------------------
                if ("concrete" %in% estimators & !("concrete" %in% out$estimates[["Package"]])) {
                  concreteOut <- try(simConcrete(TargetTime, Data, MaxIter, DiagPath),
                                     silent = TRUE)
                  if (inherits(concreteOut, "try-error")) {
                    out$meta$error <- c(out$meta$error, "concrete")
                  } else {
                    out$estimates <- concreteOut
                    out$meta$time[Estimator == "concrete", Time := attr(concreteOut, "time")]
                    if (!attr(concreteOut, "converge"))
                      out$meta$notconverge <- c(out$meta$notconverge, "concrete")
                  }
                  rm(concreteOut); gc()
                }
                
                # contmle ------------------------------------------------------------------------------
                if ("contmle" %in% estimators & !("contmle" %in% out$estimates[["Package"]])) {
                  contmleOut <- try(simContmle(TargetTime, Data, MaxIter, DiagPath),
                                    silent = TRUE)
                  
                  if (inherits(contmleOut, "try-error")) {
                    out$contmle <- data.table()
                    out$meta$error <- c(out$meta$error, "contmle")
                  } else {
                    out$contmle <- formatContmle(contmleOut)
                    out$meta$time[Estimator == "contmle", Time := attr(contmleOut, "time")]
                    if (contmleOut$convergenced.at.step >= MaxIter)
                      out$meta$notconverge <- c(out$meta$notconverge, "contmle")
                  }
                  rm(contmleOut); gc()
                }
                
                
                # survtmle ----------------------------------------------------------------
                
                for (estimator in names(discretizations)) {
                  if (estimator %in% estimators & !(estimator %in% out$estimates[["Estimator"]])) {
                    Bins <- discretizations[[estimator]]
                    survtmle_out <- try(simSurvtmle(TargetTime, Bins, Data, MaxIter, DiagPath), 
                                        silent = TRUE)
                    
                    if (inherits(survtmle_out, "try-error")) {
                      out$meta$error <- c(out$meta$error, estimator)
                    } else {
                      out$estimates <- rbind(out$estimates,
                                             cbind(Package = estimator,
                                                   survtmle_out, 
                                                   Estimand = "Abs Risk", 
                                                   "CI Low" = NA, 
                                                   "CI Hi"  = NA, 
                                                   "SimCI Low" = NA, 
                                                   "SimCI Hi"  = NA))
                      out$meta$time[Estimator == estimator,
                                    Time := attr(survtmle_out, "time")]
                    }
                  }
                  rm(survtmle_out); gc()
                }
                
                
                # Aalen-Johansen ----------------------------------------------------------
                if ("Aalen-Johansen" %in% estimators & !("Aalen-Johansen" %in% out$estimates[["Estimator"]])) {
                  
                  AJ <- try(simAJ(TargetTime = TargetTime, Data = Data, DiagPath = DiagPath))
                  
                  if (inherits(AJ, "try-error")) {
                    out$meta$error <- c(out$meta$error, "Aalen-Johansen")
                  } else {
                    out$estimates <- rbind(out$estimates, AJ)
                    out$meta$time[Estimator == "Aalen-Johansen",
                                  Time := attr(AJ, "time")]
                  }
                  rm(AJ); gc()
                }
                
                # return ------------------------------------------------------------------
                saveRDS(object = out, file = paste0(OutputPath, i, ".RDS"))
                rm(out); gc()
                return(NULL)
              }
    }, silent = TRUE)
    j <- length(list.files(OutputPath))
    registerDoSEQ()
    stopCluster(cl)
    rm(cl); gc()
    if (!inherits(sim.result, "try-error")) {
      Stop <- Sys.time()
      print(difftime(Stop, Start, units = "hours"))
    }
  }
}

# Processing --------------------------------------------------------------
Output <- lapply(list.files(OutputPath, full.names = TRUE),
                 function(run) readRDS(file = run))

estimates <- do.call(rbind, lapply(seq_along(Output), function(iter) {
  if (inherits(Output[[iter]]$estimates, "data.frame"))
    return(cbind("iter" = iter, Output[[iter]]$estimates))
  NULL
}))
estimates[grepl("survtmle", Package), se := sqrt(se)]

meta <- data.table("Estimator" = setdiff(unique(estimates[["Estimator"]]), "gcomp"), 
                   "Time" = 0, 
                   "Error" = 0, 
                   "NotConverge" = 0)
for (i in seq_along(Output)) {
  meta[, Time := Time + Output[[i]]$meta$time$Time]
  meta[, Error := Error + Estimator %in% Output[[i]]$meta$error]
  meta[, NotConverge := NotConverge + Estimator %in% Output[[i]]$meta$notconverge]
}
meta[, 2:4] <- meta[, 2:4] / length(Output)
meta

# True Risks --------------------------------------------------------------

TargetTimes <- sort(unique(estimates[["Time"]]))
setDT(risks)
risks <- rbind(risks, 
               risks[, list("Intervention" = "RD",
                            "True" = True[Intervention == "A=1"] - True[Intervention == "A=0"]),
                     by = c("Event", "Time")],
               risks[, list("Intervention" = "RR",
                            "True" = True[Intervention == "A=1"] / True[Intervention == "A=0"]),
                     by = c("Event", "Time")])

estimates <- rbind(estimates,
                   estimates[, list("Intervention" = "RD",
                                    "Pt Est" = `Pt Est`[Intervention == "A=1"] - `Pt Est`[Intervention == "A=0"],
                                    "se" = sqrt(se[Intervention == "A=1"]^2 + se[Intervention == "A=0"]^2)),
                             by = c("iter", "Package", "Estimator", "Event", "Time")],
                   estimates[, list("Intervention" = "RR",
                                    "Pt Est" = `Pt Est`[Intervention == "A=1"] / `Pt Est`[Intervention == "A=0"],
                                    "se" = sqrt((se[Intervention == "A=1"] / `Pt Est`[Intervention == "A=0"])^2 +
                                                  (se[Intervention == "A=0"] * `Pt Est`[Intervention == "A=1"] /
                                                     `Pt Est`[Intervention == "A=0"]^2)^2)),
                             by = c("iter", "Package", "Estimator", "Event", "Time")],
                   fill = TRUE)

# # for including contmle
# estimates <- rbind(estimates,
#                    do.call(rbind, lapply(seq_along(Output), function(iter) {
#                        if (inherits(Output[[iter]]$contmle, "data.frame")) {
#                            out <- setDT(Output[[iter]]$contmle)
#                            out[, Intervention := "RD"]
#                            setnames(out, "RD", "Pt Est")
#                            return(cbind("iter" = iter, out))
#                        }
#                        NULL
#                    })))



perf <- left_join(estimates[, Event := as.factor(Event)], 
                  mutate(risks, Event = factor(Event)), 
                  by = c("Intervention", "Event", "Time")) %>% 
  .[, list(l = quantile(`Pt Est`, 0.025), u = quantile(`Pt Est`, 0.975), 
           True = mean(True), 
           Bias = mean(`Pt Est` - True), 
           MSE = mean((`Pt Est` - True)^2), 
           "95% Cov" = mean(`Pt Est` + 1.96*se > True & `Pt Est` - 1.96*se < True), 
           `Pt Est` = mean(`Pt Est`), se = mean(se)),
    by = c("Package", "Intervention", "Estimator", "Event", "Time")] %>% 
  group_by(Package, Intervention, Estimator, Event, Time) %>% 
  dplyr::mutate(l.se = `Pt Est` - 1.96*se, u.se = `Pt Est` + 1.96*se) %>% 
  setDT()

# plots -------------------------------------------------------------------

melt(perf[, Time := as.numeric(Time)][, Bias := abs(Bias)], 
     id.vars = c("Package", "Intervention", "Estimator", "Event", "Time"), 
     measure.vars = c("Bias", "MSE", "95% Cov"), 
     variable.name = "Metric", value.name = "Value") %>% 
  mutate(Estimator = factor(Estimator, 
                            levels = c("Aalen-Johansen", "tmle", "survtmle-1mo", 
                                       "survtmle-3mo", "survtmle-6mo"), 
                            labels = c("Aalen-Johansen", "contmle", "survtmle-1mo", 
                                       "survtmle-3mo", "survtmle-6mo"))) %>% 
  dplyr::filter(Estimator != "gcomp") %>% 
  ggplot(aes(x = Time, y = Value, colour = Estimator)) +
  facet_wrap(Intervention + Event ~ Metric, scales = "free", ncol = 3) +
  geom_line() + theme_minimal()
ggsave(filename = "/Shared/Projects/concrete/scripts/sim-data/perf_plot.png", 
       device = "png", width = 1600, height = 1500, units = "px")

perf[, Time := as.factor(Time)] %>% dplyr::filter(Estimator != "gcomp") %>% 
  mutate(Estimator = factor(Estimator, 
                            levels = c("Aalen-Johansen", "tmle", "survtmle-1mo", 
                                       "survtmle-3mo", "survtmle-6mo"), 
                            labels = c("Aalen-Johansen", "contmle", "survtmle-1mo", 
                                       "survtmle-3mo", "survtmle-6mo"))) %>% 
  ggplot(aes(x = Time, y = `Pt Est`, colour = Estimator)) +
  facet_wrap(Intervention ~ Event, scales = "free", ncol = 2) +
  geom_errorbar(aes(ymin = l, ymax = u), width = 0.8, position = position_dodge2(width = 0.8)) +
  # geom_errorbar(aes(ymin = l.se, ymax = u.se), width = 0.8, 
  #               position = position_dodge2(width = 0.8), alpha = 0.5) +
  geom_point(size = 2, position = position_dodge2(width = 0.8)) + theme_minimal() + 
  geom_errorbar(aes(ymin = True, ymax = True))

# perf[, RelEff :=round(`se`[Estimator == "Aalen-Johansen"]^2 / `se`^2, 2),
#      by = c("Analysis", "Time", "Event", "Estimand")]

