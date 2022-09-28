library(survtmle); library(MOSS); library(SuperLearner)
surv_tmle.dir <- c("/Shared/Projects/Roadmap_CVOT/R/functions/")
x <- lapply(surv_tmle.dir, function(dir) lapply(list.files(dir, full.names = TRUE),
                                                function(x) try(source(x), silent = TRUE)))
library(survival); library(zoo)
contmle.dir <- c("/Shared/Projects/continuousTMLE/R", "~/research/SoftWare/continuousTMLE/R")
x <- lapply(contmle.dir, function(dir) lapply(list.files(dir, full.names = TRUE), source))

library(tidyverse); library(data.table); library(concrete)
source("./scripts/sim-data/sim_functions.R")
options(warn = -1)

# true risks --------------------------------------------------------------
# true risks for 3 competing risks
# 5e6 takes me ~60Gb memory
risks1 <- getTrueRisks(n = 5e6, assign_A = function(W, n) return(rep_len(1, n)))
gc()
risks0 <- getTrueRisks(n = 5e6, assign_A = function(W, n) return(rep_len(0, n)))
gc()

risks <- rbind(cbind("trt" = "A=1",
                     rbind(data.table("Event" = 1, "Time" = 1:nrow(risks1), True = risks1[[1]]),
                           data.table("Event" = 2, "Time" = 1:nrow(risks1), True = risks1[[2]]),
                           data.table("Event" = 3, "Time" = 1:nrow(risks1), True = risks1[[3]]))),
               cbind("trt" = "A=0",
                     rbind(data.table("Event" = 1, "Time" = 1:nrow(risks1), True = risks0[[1]]),
                           data.table("Event" = 2, "Time" = 1:nrow(risks1), True = risks0[[2]]),
                           data.table("Event" = 3, "Time" = 1:nrow(risks1), True = risks0[[3]]))))
risks %>% mutate(Event = factor(Event)) %>%
    ggplot(aes(x = Time, y = True, colour = Event, linetype = trt)) + geom_line(size = .8) +
    theme_minimal() + labs(title = "True Competing Risks")
rm(list = c("risks0", "risks1"))


# simulation --------------------------------------------------------------
set.seed(0)
seeds <- sample(0:12345678, size = 400)
library(foreach)
library(doParallel)
cl <- makeForkCluster(nnodes = 16, outfile = "")
registerDoParallel(cl)

Start <- Sys.time()
output <- foreach(i = seq_along(seeds),
                  .errorhandling = "remove") %dopar% {
                      out <- list()
                      out$estimates <- data.table()
                      out$concrete.notconverge <-
                          out$concrete.error <-
                          out$contmle.error <-
                          out$survtmle.6mo.notconverge <-
                          out$survtmle.6mo.error <-
                          out$survtmle.3mo.notconverge <-
                          out$survtmle.3mo.error <- 0

                      data <- simConCR(n = 8e2, random_seed = seeds[i])
                      TargetTime <- seq(730, 1460, length.out = 5)
                      target.event <- sort(setdiff(unique(data$EVENT), 0))

                      # concrete ------------------------------------------------------------------------------------
                      concreteArgs <- formatArguments(DataTable = data,
                                                      EventTime = "TIME",
                                                      EventType = "EVENT",
                                                      Treatment = "ARM",
                                                      ID = "id",
                                                      Intervention = 0:1,
                                                      TargetTime = TargetTime,
                                                      MaxUpdateIter = 50,
                                                      Verbose = FALSE)
                      concreteArgs <- formatArguments(concreteArgs)

                      concreteEst <- list();
                      iter <- 1
                      class(concreteEst) <- "try-error"
                      while (inherits(concreteEst, "try-error") & iter <= 3) {
                          Start <- Sys.time()
                          concreteEst <- try(doConcrete(concreteArgs))
                          if (inherits(concreteEst, "try-error")) {
                              out$concrete.error <- out$concrete.error + 1
                          }
                          iter <- iter + 1
                      }
                      concreteOut <- try(getOutput(concreteEst)$Risk)
                      Stop <- Sys.time()

                      if (!inherits(concreteOut, "try-error")) {
                          out$concrete.notconverge <- all(attr(concreteEst, "TmleConverged") == FALSE)
                          out$estimates <- rbind(out$estimates,
                                                 cbind(Package = "concrete", concreteOut))
                      }

                      attr(out, "concrete.time") <- Stop - Start
                      write_lines(paste0("Run ", i, " ; concrete ; ", Stop - Start),
                                  file = "scripts/sim-data/sim-diagnostics.txt", append = TRUE)
                      # contmle -------------------------------------------------------------------------------------
                      run <- list()
                      class(run) <- "try-error"

                      while (inherits(run, "try-error") & iter <= 3) {
                          Start <- Sys.time()
                          run <- try(contmle(
                              data, #-- dataset
                              target = target.event, #-- go after cause 1 and cause 2 specific risks
                              iterative = FALSE, #-- use one-step tmle to target F1 and F2 simultaneously
                              treat.effect = "ate", #-- target the ate directly
                              tau = TargetTime, #-- time-point of interest
                              estimation = list(
                                  "cens" = list(fit = "cox",
                                                model = Surv(TIME, EVENT == 0) ~ ARM+GEOGR1+SEX+AGE+EGFMDRBC+CREATBL+BMIBL+HBA1CBL+ETHNIC+MIFL+SMOKER+STROKSFL),
                                  "cause1" = list(fit = "cox",
                                                  model = Surv(TIME, EVENT == 1) ~ ARM+GEOGR1+SEX+AGE+EGFMDRBC+CREATBL+BMIBL+HBA1CBL+ETHNIC+MIFL+SMOKER+STROKSFL),
                                  "cause2" = list(fit = "cox",
                                                  model = Surv(TIME, EVENT == 2) ~ ARM+GEOGR1+SEX+AGE+EGFMDRBC+CREATBL+BMIBL+HBA1CBL+ETHNIC+MIFL+SMOKER+STROKSFL),
                                  "cause3" = list(fit = "cox",
                                                  model = Surv(TIME, EVENT == 3) ~ ARM+GEOGR1+SEX+AGE+EGFMDRBC+CREATBL+BMIBL+HBA1CBL+ETHNIC+MIFL+SMOKER+STROKSFL)),
                              treat.model = ARM ~ GEOGR1+SEX+AGE+EGFMDRBC+CREATBL+BMIBL+HBA1CBL+ETHNIC+MIFL+SMOKER+STROKSFL,
                              verbose = FALSE,
                              sl.models = NULL,
                              no.small.steps = 50))
                          if (inherits(run, "try-error")) {
                              out$contmle <- "contmle error"
                              out$contmle.error <- out$contmle.error + 1
                              iter <- iter + 1
                          } else {
                              out$contmle <- try(formatContmle(run))
                          }
                      }
                      Stop <- Sys.time()
                      attr(out, "contmle.time") <- Stop - Start
                      write_lines(paste0("Run ", i, " ; contmle ; ", Stop - Start),
                                  file = "scripts/sim-data/sim-diagnostics.txt", append = TRUE)
                      # survtmle 6mo ----------------------------------------------------------------
                      sl_lib_failure <-
                          sl_lib_censor <-
                          sl_lib_g <- c("SL.glmnet", "SL.xgboost")

                      intervals.per.year <- 2
                      target.time <- ceiling(TargetTime / 365 * intervals.per.year)

                      tmle_sl <- list()
                      iter <- 1
                      class(tmle_sl) <- "try-error"
                      while (inherits(tmle_sl, "try-error") & iter <= 3) {
                          Start <- Sys.time()
                          tmle_sl <- try(surv_tmle(
                              ftime = ceiling(data$TIME / 365 * intervals.per.year),
                              ftype = data$EVENT,
                              targets = target.time,
                              trt = data$ARM,
                              t0 = max(target.time),
                              adjustVars = as.data.frame(data[, !c("id", "TIME", "EVENT", "ARM")]),
                              SL.ftime = sl_lib_failure,
                              SL.ctime = sl_lib_censor,
                              SL.trt = sl_lib_g,
                              # glm.trt = glm_trt,
                              returnIC = TRUE,
                              returnModels = TRUE,
                              ftypeOfInterest = target.event,
                              trtOfInterest = c(1, 0),
                              maxIter = 50,
                              method = "hazard",
                              verbose = FALSE
                          ))
                          if (inherits(tmle_sl, "try-error")) {
                              out$survtmle.6mo.error <- out$survtmle.6mo.error + 1
                          }
                          iter <- iter + 1
                      }
                      if (!inherits(tmle_sl, "try-error")) {
                          if (!tmle_sl$converged) {
                              out$survtmle.6mo.notconverge <- out$survtmle.6mo.notconverge + 1
                          }

                          survtmle.est <- rbind(cbind("Estimator" = "gcomp6mo",
                                                      "Intervention" = rep(c("A=0", "A=1"), times = length(target.event)),
                                                      "Event" = rep(target.event, each = 2),
                                                      as.data.frame(tmle_sl$init_est)),
                                                cbind(Estimator = "survtmle6mo",
                                                      "Intervention" = rep(c("A=0", "A=1"), times = length(target.event)),
                                                      Event = rep(target.event, each = 2),
                                                      as.data.frame(tmle_sl$est))) %>% as.data.table() %>%
                              melt(., id.vars = c("Estimator", "Intervention", "Event"), value.name = "Risk", variable.name = "Time")
                          survtmle.est[, Time := as.numeric(gsub("\\D+", "", Time)) / intervals.per.year * 365]
                          survtmle.est <- full_join(survtmle.est,
                                                    cbind(Estimator = "survtmle6mo",
                                                          "Intervention" = rep(c("A=0", "A=1"), each = length(target.event) * length(target.time)),
                                                          Event = rep(target.event, each = length(target.time), times = 2),
                                                          Time = rep(TargetTime, times = 2 * length(target.event)),
                                                          data.table(se = sqrt(diag(tmle_sl$var)))))
                          out$estimates <- rbind(out$estimates,
                                                 cbind(Package = "survtmle",
                                                       survtmle.est))
                      }
                      Stop <- Sys.time()
                      attr(out, "survtmle6.time") <- Stop - Start
                      write_lines(paste0("Run ", i, " ; survtmle6 ; ", Stop - Start),
                                  file = "scripts/sim-data/sim-diagnostics.txt", append = TRUE)

                      # survtmle 3mo ----------------------------------------------------------------
                      intervals.per.year <- 4
                      target.time <- ceiling(TargetTime / 365 * intervals.per.year)
                      target.event <- order(setdiff(unique(data$EVENT), 0))

                      tmle_sl <- list()
                      iter <- 1
                      class(tmle_sl) <- "try-error"
                      while (inherits(tmle_sl, "try-error") & iter <= 3) {
                          Start <- Sys.time()
                          tmle_sl <- try(surv_tmle(
                              ftime = ceiling(data$TIME / 365 * intervals.per.year),
                              ftype = data$EVENT,
                              targets = target.time,
                              trt = data$ARM,
                              t0 = max(target.time),
                              adjustVars = as.data.frame(data[, !c("id", "TIME", "EVENT", "ARM")]),
                              SL.ftime = sl_lib_failure,
                              SL.ctime = sl_lib_censor,
                              SL.trt = sl_lib_g,
                              # glm.trt = glm_trt,
                              returnIC = TRUE,
                              returnModels = TRUE,
                              ftypeOfInterest = target.event,
                              trtOfInterest = c(1, 0),
                              maxIter = 50,
                              method = "hazard",
                              verbose = FALSE
                          ))
                          if (inherits(tmle_sl, "try-error")) {
                              out$survtmle.3mo.error <- out$survtmle.3mo.error + 1
                          }
                          iter <- iter + 1
                      }
                      if (!inherits(tmle_sl, "try-error")) {
                          if (!tmle_sl$converged) {
                              out$survtmle.3mo.notconverge <- out$survtmle.3mo.notconverge + 1
                          }
                          survtmle.est.3mo <- rbind(cbind("Estimator" = "gcomp3mo",
                                                          "Intervention" = rep(c("A=0", "A=1"), times = length(target.event)),
                                                          "Event" = rep(target.event, each = 2),
                                                          as.data.frame(tmle_sl$init_est)),
                                                    cbind(Estimator = "survtmle3mo",
                                                          "Intervention" = rep(c("A=0", "A=1"), times = length(target.event)),
                                                          Event = rep(target.event, each = 2),
                                                          as.data.frame(tmle_sl$est))) %>% as.data.table() %>%
                              melt(., id.vars = c("Estimator", "Intervention", "Event"), value.name = "Risk", variable.name = "Time")
                          survtmle.est.3mo[, Time := as.numeric(gsub("\\D+", "", Time)) / intervals.per.year * 365]
                          survtmle.est.3mo <- full_join(survtmle.est.3mo,
                                                        cbind(Estimator = "survtmle3mo",
                                                              "Intervention" = rep(c("A=0", "A=1"), each = length(target.event) * length(target.time)),
                                                              Event = rep(target.event, each = length(target.time), times = 2),
                                                              Time = rep(TargetTime, times = 2 * length(target.event)),
                                                              data.table(se = sqrt(diag(tmle_sl$var)))))
                          out$estimates <- rbind(out$estimates,
                                                 cbind(Package = "survtmle",
                                                       survtmle.est.3mo))
                      }

                      Stop <- Sys.time()
                      attr(out, "survtmle3.time") <- Stop - Start
                      write_lines(paste0("Run ", i, " ; survtmle3 ; ", Stop - Start),
                                  file = "scripts/sim-data/sim-diagnostics.txt", append = TRUE)

                      # return ------------------------------------------------------------------
                      out
                  }
registerDoSEQ()
stopCluster(cl)
Stop <- Sys.time()
Stop - Start


# Processing --------------------------------------------------------------

estimates <- do.call(rbind, lapply(seq_along(output), function(iter) {
    if (inherits(output[[iter]]$estimates, "data.frame"))
        return(cbind("iter" = iter, output[[iter]]$estimates))
    NULL
}))
estimates <- rbind(estimates,
                   estimates[, list("Intervention" = "RD",
                                    "Risk" = Risk[Intervention == "A=1"] - Risk[Intervention == "A=0"],
                                    "se" = sqrt(se[Intervention == "A=1"]^2 + se[Intervention == "A=0"]^2)),
                             by = c("iter", "Package", "Estimator", "Event", "Time")],
                   estimates[, list("Intervention" = "RR",
                                    "Risk" = Risk[Intervention == "A=1"] / Risk[Intervention == "A=0"],
                                    "se" = sqrt((se[Intervention == "A=1"] / Risk[Intervention == "A=0"])^2 +
                                                    (se[Intervention == "A=0"] * Risk[Intervention == "A=1"] /
                                                         Risk[Intervention == "A=0"]^2)^2)),
                             by = c("iter", "Package", "Estimator", "Event", "Time")]
)
estimates <- rbind(estimates,
                   do.call(rbind, lapply(seq_along(output), function(iter) {
                       if (inherits(output[[iter]]$contmle, "data.frame")) {
                           out <- setDT(output[[iter]]$contmle)
                           out[, Intervention := "RD"]
                           setnames(out, "RD", "Risk")
                           return(cbind("iter" = iter, out))
                       }
                       NULL
                   })))




perf <- estimates[, list(Risk = mean(Risk), se = mean(se),
                         l = quantile(Risk, 0.025), u = quantile(Risk, 0.975)),
          by = c("Package", "Intervention", "Estimator", "Event", "Time")]

perf[, Time := as.factor(Time)] %>%
    ggplot(aes(x = Time, y = Risk, colour = Estimator, shape = Package)) +
    facet_wrap(Intervention ~ Event, scales = "free") +
    geom_errorbar(aes(ymin = l, ymax = u), width = 0.8, position = position_dodge2(width = 0.8)) +
    geom_point(size = 2, position = position_dodge2(width = 0.8)) + theme_minimal()
