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
    data.table::setDTthreads(threads = 2)
    x <- try(source("/Shared/Projects/ConCR-TMLE/scripts/sim-data/sim_functions.R"))
    library(tidyverse)
    invisible(NULL)
}

# true risks --------------------------------------------------------------
# true risks for 3 competing risks
# 5e6 takes me ~60Gb memory
loadPackages()
data.table::setDTthreads(threads = parallel::detectCores(), restore_after_fork = TRUE)

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
B <- 400
seeds <- sample(0:12345678, size = B)
output <- list()
library(foreach)
library(doParallel)
j <- 1
n_cores <- min(parallel::detectCores(), 12)
while(j <= B) {
    seeds_j <- seeds[j:min(B, j + n_cores - 1)]
    j <- j + n_cores
    cl <- makeCluster(n_cores, type = "FORK")
    registerDoParallel(cl)
    Start <- Sys.time()
    sim.result <- try({
        foreach(i = seq_along(seeds_j),
                .errorhandling = "remove",
                .packages = "readr") %dopar% {
                    loadPackages()
                    DiagPath <- "scripts/sim-data/sim-diag/diagnostic.txt"
                    out <- list()
                    estimators <- c("concrete", "contmle", "survtmle-6mo", "survtmle-3mo")
                    out$estimates <- data.table()
                    out$meta <- list("time" = data.table("Estimator" = estimators, "Time" = NaN),
                                     "error" = list(), "notconverge" = list())

                    MaxIter <- 100
                    Data <- simConCR(n = 8e2, random_seed = seeds_j[i])
                    TargetTime <- seq(730, 1460, length.out = 5)
                    TargetEvent <- sort(setdiff(unique(Data$EVENT), 0))

                    # helper functions --------------------------------------------------------
                    simConcrete <- function(TargetTime, Data, MaxIter, DiagPath) {
                        require(concrete); require(data.table)
                        readr::write_lines(paste0("Run ", i, " ; concrete ; start"),
                                           file = DiagPath, append = TRUE)
                        concreteArgs <- formatArguments(DataTable = Data,
                                                        EventTime = "TIME",
                                                        EventType = "EVENT",
                                                        Treatment = "ARM",
                                                        ID = "id",
                                                        Intervention = 0:1,
                                                        TargetTime = TargetTime,
                                                        MaxUpdateIter = MaxIter,
                                                        Verbose = FALSE)
                        Start <- Sys.time()
                        concreteEst <- doConcrete(concreteArgs)
                        concreteOut <- getOutput(concreteEst)$Risk
                        Stop <- Sys.time()
                        estimates <- cbind(Package = "concrete", concreteOut)
                        attr(estimates, "converged") <- attr(concreteEst, "TmleConverged")$converged
                        attr(estimates, "time") <- difftime(Stop, Start, units = "mins")
                        readr::write_lines(paste0("Run ", i, " ; concrete ; ",
                                                  format(unclass(attr(estimates, "time")), digits = 3),
                                                  " ", attr(attr(estimates, "time"), "units")),
                                           file = DiagPath, append = TRUE)
                        return(estimates)
                    }

                    simContmle <- function(TargetTime, Data, MaxIter, DiagPath) {
                        required <- c("data.table", "zoo", "glmnet", "survival", "stringr",
                                      "nleqslv","prodlim", "Matrix", "coefplot", "hdnom")
                        z <- sapply(required, function(p) try(library(p, character.only = TRUE),
                                                              silent = TRUE))
                        readr::write_lines(paste0("Run ", i, " ; contmle ; start"),
                                           file = DiagPath, append = TRUE)
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
                                                    EGFMDRBC+CREATBL+BMIBL+HBA1CBL+ETHNIC+SMOKER+STROKSFL),
                                "cause3" = list(fit = "cox",
                                                model = Surv(TIME, EVENT == 3) ~ ARM+GEOGR1+SEX+AGE+MIFL+
                                                    EGFMDRBC+CREATBL+BMIBL+HBA1CBL+ETHNIC+SMOKER+STROKSFL)
                            ),
                            treat.model = ARM ~ GEOGR1+SEX+AGE+MIFL+EGFMDRBC+CREATBL+BMIBL+
                                HBA1CBL+ETHNIC+SMOKER+STROKSFL,
                            verbose = FALSE,
                            sl.models = NULL,
                            no.small.steps = MaxIter
                        )
                        Stop <- Sys.time()
                        attr(contmleRun, "time") <- difftime(Stop, Start, units = "mins")
                        readr::write_lines(paste0("Run ", i, " ; contmle ; ",
                                                  format(unclass(attr(contmleRun, "time")), digits = 3),
                                                  " ", attr(attr(contmleRun, "time"), "units")),
                                           file = DiagPath, append = TRUE)

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

                    simSurvtmle <- function(TargetTime, Bins, Data, MaxIter, DiagPath) {
                        require(SuperLearner)
                        survMonths <- paste0("survtmle", as.integer(12 / Bins), "mo")
                        readr::write_lines(paste0("Run ", i, " ; ", survMonths," ; start"),
                                           file = DiagPath, append = TRUE)
                        sl_lib_failure <- c("SL.glmnet", "SL.xgboost")
                        sl_lib_censor <- c("SL.glmnet", "SL.xgboost")
                        sl_lib_g <- c("SL.glmnet", "SL.xgboost")
                        TargetTime <- ceiling(TargetTime / 365 * Bins)
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
                                            ftypeOfInterest = TargetEvent,
                                            trtOfInterest = c(1, 0),
                                            method = "hazard",
                                            verbose = FALSE,
                                            tol = 1 / nrow(Data),
                                            maxIter = MaxIter,
                                            Gcomp = FALSE)
                        survtmleOut <- print(timepoints(tmleFit, times = TargetTime))
                        survtmleOut <- getSurvtmleTbl(survtmleOut, TargetEvent, Bins)

                        Stop <- Sys.time()
                        attr(survtmleOut, "time") <- difftime(Stop, Start, units = "mins")
                        readr::write_lines(paste0("Run ", i, " ; ", survMonths," ; ",
                                                  format(unclass(attr(survtmleOut, "time")), digits = 3), " ",
                                                  attr(attr(survtmleOut, "time"), "units")),
                                           file = DiagPath, append = TRUE)
                        return(survtmleOut)
                    }

                    getSurvtmleTbl <- function(survtmleOut, TargetEvent, Bins) {
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
                        })
                        do.call(dplyr::full_join, survtmle.list)
                    }

                    # concrete --------------------------------------------------------------------
                    if ("concrete" %in% estimators) {
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
                        rm(concreteOut)
                    }

                    # contmle ------------------------------------------------------------------------------
                    if ("contmle" %in% estimators) {
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
                    }

                    # survtmle 6mo ----------------------------------------------------------------
                    if ("survtmle-6mo" %in% estimators) {
                        Bins <- 2
                        survtmle6mo <- try(simSurvtmle(TargetTime, Bins, Data, MaxIter, DiagPath),
                                           silent = TRUE)

                        if (inherits(survtmle6mo, "try-error")) {
                            out$meta$error <- c(out$meta$error, paste0("survtmle", as.integer(12 / Bins), "mo"))
                        } else {
                            out$estimates <- rbind(out$estimates,
                                                   cbind(Package = "survtmle",
                                                         survtmle6mo))
                            out$meta$time[Estimator == paste0("survtmle-", as.integer(12 / Bins), "mo"),
                                          Time := attr(survtmle6mo, "time")]
                        }
                    }

                    # survtmle 3mo ----------------------------------------------------------------
                    if ("survtmle-3mo" %in% estimators) {
                        Bins <- 4
                        survtmle3mo <- try(simSurvtmle(TargetTime, Bins, Data, MaxIter, DiagPath),
                                           silent = TRUE)

                        if (inherits(survtmle3mo, "try-error")) {
                            out$meta$error <- c(out$meta$error, paste0("survtmle", as.integer(12 / Bins), "mo"))
                        } else {
                            out$estimates <- rbind(out$estimates,
                                                   cbind(Package = "survtmle",
                                                         survtmle3mo))
                            out$meta$time[Estimator == paste0("survtmle-", as.integer(12 / Bins), "mo"),
                                          Time := attr(survtmle3mo, "time")]
                        }
                    }

                    # return ------------------------------------------------------------------
                    return(out)
                }
    }, silent = TRUE)
    if (inherits(sim.result, "try-error")) {
        j <- j - n_cores
    } else {
        output <- c(output, sim.result)
    }
    registerDoSEQ()
    stopCluster(cl)
    rm(cl)
    Stop <- Sys.time()
    difftime(Stop, Start, units = "mins")
}






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
