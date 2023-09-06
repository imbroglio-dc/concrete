setwd(dir = "/Shared/Projects/novo_nordisk/")
library(tidyverse); library(mice); library(data.table); library(survival); library(concrete)
library(SuperLearner); library(survtmle)

surv_tmle.dir <- c("/Shared/Projects/novo_nordisk/R/functions/")
x <- lapply(surv_tmle.dir, function(dir) lapply(list.files(dir, full.names = TRUE),
                                                function(x) try(source(x), silent = TRUE)))

file_dir <- "data/SUSTAIN6 NN9535-3744/ADaM/"
set.seed(123445678)

# read in subject level covariate dataset
W <- haven::read_sas(paste0(file_dir, "adsl.sas7bdat")) %>%
    dplyr::select_at(.vars = vars(!starts_with("CRIT"))) %>%
    dplyr::select(-c("AGEG1N", "AGEG1", "AGEBLG1N", "AGEBLG1", "AGEBLG2N", "AGEBLG2", "SEX2",
                     "RACEGR1", "BMIBLU", "BMIG1BL", "BMIG1BLN", "TRTSEQPN", "TRTSEQAN",
                     "TRT30PN", "TRT30AN", "TR30PG1N", "TR30AG1N", "TR30PG2N", "TR30AG2N",
                     "TRT90AN", "ARMN", "ARMCD", "DURBLGR1", "DURBGR1N", "INTRLY", "ONTRTY",
                     "ONTWRY", "TRDURY", "AGEBLGR2", "BMIBLGR2", "A1CBLGR1", "A1CBGR1N",
                     "HBAG1BL", "MSSUBJID", "OUSUBJID", "RFICDT", "SUBJID", "AGEBL", "WDTRS",
                     "WDTFL", "DTRRS", "ACTARMN", "ACTARMCD", "STRRNLFL", "RNLIMPFL",
                     "RETINBSL", "RETBSLFL", "RETINSEV", "EGFRG1BL", "SMRIMPFL", "SRICKDFL",
                     "SMICKDFL", "TRTSDT")) %>%
    dplyr::select_if(~length(unique(.))>1)

outcomes <- haven::read_sas(paste0(file_dir, "adtte.sas7bdat")) %>%
    dplyr::select(-c("STUDYID", "ASEQ", "EACFL", "TOPIC_CD", "ARM", "ARMCD", "TRTP", "TRTPN",
                     "TRTA", "TRTAN", "TRTPG1N", "TRTAG1N", "TRTPG2N", "TRTAG2N", "TRTSEQP",
                     "TRTSEQA", "ATMF", "AVALC", "AVALU", "TRDURY", "TRTSDTF", "TRTEDTF",
                     "AGEU", "AGEBLU", "AGEBLG1N", "AGEBLG1", "AGEBLGR2", "SEX", "RACEGR1",
                     "HBA1CBLU", "A1CBLGR1", "A1CBGR1N", "BMIBLU", "BMIG1BL", "BMIG1BLN",
                     "BMIBLGR2", "RNLIMPFL", "SMRIMPFL", "SRICKDFL", "SMICKDFL", "EGFMDRDU",
                     "EGFEPBCU", "DIABDURU", "DURBLGR1", "DURBGR1N", "RETBSLFL", "ADT", "ADTM",
                     "age_category", "BWGTBLU", "COUNTRY2", "ETHNIC", "ETHNIC2", "GEOGR1",
                     "SUBJID", "TRTSDT", "TRTEDT", "STARTDT", "SITEID", "REGION", "RANDDT",
                     "RACEOTH", "RACEGR1N", "GEOGR1N")) %>%
    dplyr::select_if(~length(unique(.)) > 1)
setDT(outcomes)
# View(outcomes[, sum(CNSR != 1 & ANL02FL == "Y"), by = c("PARAM", "PARAMCD")])

doSus6Demo <- function(W, outcomes, targets = 365.25*seq(0.5, 2, 0.5),
                       analysis = c("All", "AdvSideEffects", "NonGastroASE"),
                       binwidth = c(1, 3, 6), seed = 123456789) {
    set.seed(seed)
    out <- lapply(analysis, function(definition) {
        obs <- outcomes %>%
            dplyr::filter(PARAMCD == "TMTDMIST", ANL02FL == "Y") %>%
            dplyr::select(USUBJID, CNSR, AVAL) %>%
            left_join(., W, by = "USUBJID") %>%
            setDT()
        
        definition <- "AdvSideEffects"
        reasons <- list("All" = setdiff(obs$DTRRSSP, ""),
                        "AdvSideEffects"= c("GASTROINTESTINAL TOLERABILITY",
                                            "ADVERSE EVENT OTHER THAN RELATED TO GASTROINTESTINAL TOLERABILITY"),
                        "NonGastroASE" = "ADVERSE EVENT OTHER THAN RELATED TO GASTROINTESTINAL TOLERABILITY")
        reasons <- reasons[[definition]]
        
        obs[, status := 0][, time := ONTRTD]
        obs[DTRRSSP %in% reasons, status := 2]
        obs[CNSR != 1 & AVAL <= ONTRTD, status := 1]
        obs[status == 1, time := AVAL]
        obs[, trt := as.numeric(TR30PG2 == "Sema")]
        obs[, event1 := as.numeric(status == 1)]
        obs[, event2 := 2*as.numeric(status == 2)]
        
        # Aalen-Johanssen -----------------------------------------------------------------------------
        
        aj <- survfit(Surv(time = time, event = as.factor(status)) ~ trt, data = obs, ctype = 1) %>%
            summary(., times = targets) %>% .$pstate %>% as.data.table() %>% .[, -1]
        setnames(aj, paste0("Risk", 1:ncol(aj)))
        aj <- survfit(Surv(time = time, event = as.factor(status)) ~ trt, data = obs, ctype = 1) %>%
            summary(., times = targets) %>% .$std.err %>% as.data.table() %>% .[, -1] %>%
            setnames(., new = paste0("se", 1:ncol(aj))) %>% cbind(aj, .) %>%
            cbind("Time" = targets, "A" = rep(0:1, each = length(targets))) %>%
            melt(., id.vars = c("Time", "A")) %>%
            .[, "Event" := stringr::str_extract(string = variable, pattern = "[0-9]+")] %>%
            .[, variable := stringr::str_extract(string = variable, pattern = "[^(0-9)]+")] %>%
            dcast(., Event + A + Time ~ variable) %>%
            .[, list("R1" = Risk[A == 1],
                     "R1_se" = se[A == 1],
                     "R0" = Risk[A == 0],
                     "R0_se" = se[A == 0],
                     "RD" = Risk[A == 1] -  Risk[A == 0],
                     "RD_se" = sqrt(se[A == 1]^2 + se[A == 0]^2),
                     "RR" = Risk[A == 1] / Risk[A == 0],
                     "RR_se" = sqrt((se[A == 1] / Risk[A == 0])^2 +
                                        se[A == 0]^2 * (Risk[A == 1] / Risk[A == 0]^2)^2)
            ),
            by = c("Event", "Time")] %>%
            melt(., id.vars = c("Event", "Time")) %>%
            .[, "Estimand" := stringr::str_extract(variable, "[A-Z0-9]+")] %>%
            .[, "type" := stringr::str_extract(variable, "[^(A-Z0-9_)]+")] %>%
            .[is.na(type), type := "Pt Est"] %>%
            .[, -"variable"] %>%
            dcast(., Event + Time + Estimand ~ type, value.var = "value") %>%
            .[, c("CI Low", "CI Hi", "Analysis", "Estimator") :=
                  list(`Pt Est` - 1.96*se, `Pt Est` + 1.96*se, "CompRisk", "NPMLE")]
        
        # KM ------------------------------------------------------------------------------------------
        
        km <- survfit(Surv(time = time, event = event1) ~ trt, data = obs) %>% summary(., times = targets) %>%
            with(., expr = cbind(Event = 1, A = rep(0:1, each = length(targets)), Risk = 1 - surv, se = std.err))
        km <- survfit(Surv(time = time, event = event2==2) ~ trt, data = obs) %>% summary(., times = targets) %>%
            with(., expr = cbind(Event = 2, A = rep(0:1, each = length(targets)), Risk = 1 - surv, se = std.err)) %>%
            rbind(km, .) %>% as.data.table() %>% cbind(Time = rep(targets)) %>%
            .[, list("R1" = Risk[A == 1],
                     "R1_se" = se[A == 1],
                     "R0" = Risk[A == 0],
                     "R0_se" = se[A == 0],
                     "RD" = Risk[A == 1] -  Risk[A == 0],
                     "RD_se" = sqrt(se[A == 1]^2 + se[A == 0]^2),
                     "RR" = Risk[A == 1] / Risk[A == 0],
                     "RR_se" = sqrt((se[A == 1] / Risk[A == 0])^2 +
                                        se[A == 0]^2 * (Risk[A == 1] / Risk[A == 0]^2)^2)),
              by = c("Event", "Time")]  %>%
            melt(., id.vars = c("Event", "Time")) %>%
            .[, "Estimand" := stringr::str_extract(variable, "[A-Z0-9]+")] %>%
            .[, "type" := stringr::str_extract(variable, "[^(A-Z0-9_)]+")] %>%
            .[is.na(type), type := "Pt Est"] %>%
            .[, -"variable"] %>%
            dcast(., Event + Time + Estimand ~ type, value.var = "value") %>%
            .[, c("CI Low", "CI Hi", "Analysis", "Estimator") :=
                  list(`Pt Est` - 1.96*se, `Pt Est` + 1.96*se, "RightCens", "NPMLE")]
        
        
        # data processing -----------------------------------------------------------------------------
        
        
        smallobs <- dplyr::select(obs, trt, time, status, event1, event2,
                                  SEX, RACE, SMOKER, HBA1CBL, BMIBL, DIABDUR, EVCVDFL, INSTRTBS,
                                  AGE, BMIBL, DIABDUR, EGFMDRD, HBA1CBL, SUMONTFL, SUFL, ETHNICFL,
                                  EGFRTXT, PMYCINFL, PSSTIAFL, PRREVAFL, PSTN50FL, PSCOHDFL,
                                  PASCISFL, PCKIDSFL, PMICPRFL, PHLVHYFL, PLVSDSFL, PABRINFL,
                                  CV90DYFL, STRATUM, STRATUM2, STRATUM3, STRATUM4, STRATUM5,
                                  CREAG1BL, ACRG1BL, MALBG1BL, ANTDBFL, EGFMDRD, EGFREPBC, CKDEPTXT,
                                  RETISEVN, AGE, COUNTRY, SITEID) %>%
            mutate(ACRG1BL = case_when(ACRG1BL == "" ~ "NA",
                                       T ~ ACRG1BL)) %>%
            # removing covariates with high missingness
            dplyr::select(-c("MALBG1BL", "PABRINFL", "PASCISFL", "PCKIDSFL", "PHLVHYFL", "PLVSDSFL",
                             "PMICPRFL", "PMYCINFL", "PRREVAFL", "PSCOHDFL", "PSSTIAFL", "PSTN50FL")) %>%
            # removing covariates with high collinearity
            dplyr::select(-c(STRATUM, STRATUM2, STRATUM3, STRATUM4, STRATUM5)) %>%
            dplyr::mutate_if(.predicate = is.character, ~ as.factor(.))
        setDT(smallobs)
        set.seed(0)
        imp <- mice::mice(smallobs[, .SD, .SDcols = !c("trt", "time", "status", "event1", "event2")],
                          m = 1, maxit = 10, printFlag = FALSE)
        smallobs[, BMIIMP := as.numeric(is.na(BMIBL))][, BMIBL := complete(imp)[["BMIBL"]]]
        
        
        # survtmle ------------------------------------------------------------------------------------
        
        
        screeners <- "All"
        
        sl_lib_g <- expand.grid(c("SL.glmnet", "SL.xgboost"), screeners)
        sl_lib_g <- lapply(1:nrow(sl_lib_g),
                           function(i) as.character(unlist(sl_lib_g[i,])))
        
        sl_lib_censor <- expand.grid(c("SL.glmnet", "SL.xgboost"), screeners)
        sl_lib_censor <- lapply(1:nrow(sl_lib_censor),
                                function(i) as.character(unlist(sl_lib_censor[i,])))
        
        sl_lib_failure <- expand.grid(c("SL.glmnet", "SL.xgboost"), screeners)
        sl_lib_failure <- lapply(1:nrow(sl_lib_failure),
                                 function(i) as.character(unlist(sl_lib_failure[i,])))
        
        survtmle_cr <- lapply(binwidth, function(bin) {
            tmle_sl <- try(try-error)
            errors <- 0
            tmp <- smallobs
            smallobs <- smallobs[1:300, ]
            while(inherits(tmle_sl, "try-error")) {
                window <- 365.25 / 12 * bin
                debugonce(stats::optim)
                tmle_sl <- try(surv_tmle(
                    ftime = ceiling(smallobs$time/window),
                    ftype = smallobs$status,
                    targets = ceiling(targets/window),
                    trt = smallobs$trt,
                    t0 = max(ceiling(targets/window)),
                    adjustVars = smallobs[, !c("time", "status", "trt", "event1", "event2")],
                    SL.ftime = sl_lib_failure,
                    SL.ctime = sl_lib_censor,
                    SL.trt = sl_lib_g,
                    returnIC = TRUE,
                    returnModels = TRUE,
                    ftypeOfInterest = c(1, 2),
                    trtOfInterest = c(1, 0),
                    maxIter = 500,
                    method = "hazard",
                    verbose = TRUE))
                if (inherits(tmle_sl, "try-error"))
                    errors <- errors + 1
            }
            attr(tmle_sl, "errors") <- errors
            return(tmle_sl)
        })
        
        survtmle_j1 <- lapply(binwidth, function(bin) {
            tmle_sl <- try(try-error)
            errors <- 0
            while(inherits(tmle_sl, "try-error")) {
                window <- 365.25 / 12 * bin
                tmle_sl <- try(surv_tmle(
                    ftime = ceiling(smallobs$time/window),
                    ftype = smallobs$event1,
                    targets = ceiling(targets/window),
                    trt = smallobs$trt,
                    t0 = max(ceiling(targets/window)),
                    adjustVars = smallobs[, !c("time", "status", "trt", "event1", "event2")],
                    SL.ftime = sl_lib_failure,
                    SL.ctime = sl_lib_censor,
                    SL.trt = sl_lib_g,
                    returnIC = TRUE,
                    returnModels = TRUE,
                    ftypeOfInterest = c(1),
                    trtOfInterest = c(1, 0),
                    maxIter = 500,
                    method = "hazard",
                    verbose = TRUE))
                if (inherits(tmle_sl, "try-error"))
                    errors <- errors + 1
            }
            attr(tmle_sl, "errors") <- errors
            return(tmle_sl)
        })
        
        survtmle_j2 <- lapply(binwidth, function(bin) {
            tmle_sl <- try(try-error)
            errors <- 0
            while(inherits(tmle_sl, "try-error")) {
                window <- 365.25 / 12 * bin
                tmle_sl <- try(surv_tmle(
                    ftime = ceiling(smallobs$time/window),
                    ftype = smallobs$event2,
                    targets = ceiling(targets/window),
                    trt = smallobs$trt,
                    t0 = max(ceiling(targets/window)),
                    adjustVars = smallobs[, !c("time", "status", "trt", "event1", "event2")],
                    SL.ftime = sl_lib_failure,
                    SL.ctime = sl_lib_censor,
                    SL.trt = sl_lib_g,
                    returnIC = TRUE,
                    returnModels = TRUE,
                    ftypeOfInterest = c(1),
                    trtOfInterest = c(1, 0),
                    maxIter = 500,
                    method = "hazard",
                    verbose = TRUE))
                if (inherits(tmle_sl, "try-error"))
                    errors <- errors + 1
            }
            attr(tmle_sl, "errors") <- errors
            return(tmle_sl)
        })
        
        survtmle.est <- rbind(cbind("Estimator" = "d.gcomp",
                                    "A" = rep(0:1, times = length(target.event)),
                                    "Event" = rep(target.event, each = 2),
                                    as.data.frame(tmle_sl$init_est)),
                              cbind(Estimator = "survtmle",
                                    A = rep(0:1, times = length(target.event)),
                                    Event = rep(target.event, each = 2),
                                    as.data.frame(tmle_sl$est))) %>% as.data.table() %>%
            melt(., id.vars = c("Estimator", "A", "Event"), value.name = "Risk", variable.name = "Time")
        survtmle.est[["Time"]] <- rep(targets, each = length(target.event)*4)
        survtmle.est <- full_join(survtmle.est,
                                  cbind(Estimator = "survtmle",
                                        A = rep(0:1, each = length(target.event) * length(targets)),
                                        Event = rep(target.event, each = length(targets), times = 2),
                                        Time = rep(targets, times = 2 * length(target.event)),
                                        data.table(se = sqrt(diag(tmle_sl$var))))
        )
        
        survtmle.est <- dcast(survtmle.est, ... ~ A, value.var = c("Risk", "se"))
        
        # concrete sample -----------------------------------------------------------------------------
        
        cr.args <- formatArguments(DataTable = smallobs[, !c("event1", "event2")],
                                   EventTime = "time",
                                   EventType = "status",
                                   Treatment = "trt",
                                   TargetTime = targets,
                                   Intervention = 0:1)
        cr.args$Model[["0"]][["coxnet"]] <- cr.args$Model[["1"]][["coxnet"]] <-
            cr.args$Model[["2"]][["coxnet"]] <- "coxnet"
        cr.est <- doConcrete(ConcreteArgs = cr.args)
        cr.out <- getOutput(ConcreteEst = cr.est)
        
        event1.args <- formatArguments(DataTable = smallobs[, !c("status", "event2")],
                                       EventTime = "time",
                                       EventType = "event1",
                                       Treatment = "trt",
                                       TargetTime = targets,
                                       Intervention = 0:1)
        event1.args$Model[["0"]][["coxnet"]] <- event1.args$Model[["1"]][["coxnet"]] <- "coxnet"
        event1.est <- doConcrete(event1.args)
        event1.out <- getOutput(ConcreteEst = event1.est)
        
        event2.args <- formatArguments(DataTable = smallobs[, !c("status", "event1")],
                                       EventTime = "time",
                                       EventType = "event2",
                                       Treatment = "trt",
                                       TargetTime = targets,
                                       Intervention = 0:1)
        event2.args$Model[["0"]][["coxnet"]] <- event2.args$Model[["2"]][["coxnet"]] <- "coxnet"
        event2.est <- doConcrete(event2.args)
        event2.out <- getOutput(ConcreteEst = event2.est)
        
        comp <- rbind(cbind("Analysis" = "CompRisk", cr.out),
                      cbind("Analysis" = "RightCens", event1.out),
                      cbind("Analysis" = "RightCens", event2.out))
        comp_plot <- dplyr::filter(comp, Event %in% 1:2, Estimator == "tmle")
        
        comp_plot <- comp_plot %>%
            mutate(Time = factor(Time),
                   Event = case_when(Event == 1 ~ "MACE", T ~ "Trt Disc")) %>%
            ggplot(aes(x = Time, y = `Pt Est`, colour = Analysis)) +
            theme_minimal() + facet_wrap(Intervention~Event, ncol = 2, scales = "free") +
            geom_point(position = position_dodge(.2)) +
            geom_errorbar(aes(ymin = `CI Low`, ymax = `CI Hi`),
                          position = position_dodge(.2)) + ylab("MACE Risk Ratio (TMLE)") +
            labs(title = "SUSTAIN6 Prelim Analysis: Competing Risks vs. Right Censoring",
                 subtitle = "Effect of Sema on MACE & Trt Discontinuation")
        
        est <- comp[order(Analysis, Estimator)] %>%
            mutate(Estimand = case_when(Intervention == "[A=1] - [A=0]" ~ "RD",
                                        Intervention == "[A=1] / [A=0]" ~ "RR",
                                        Intervention == "A=1" ~ "R1",
                                        Intervention == "A=0" ~ "R0",
                                        T ~ "NA")) %>%
            dplyr::select(Analysis, Estimand, Estimator, Time, Event, `Pt Est`,
                          se, `CI Low`, `CI Hi`) %>%
            rbind(., aj, km)
        comp_plot2 <- est %>% filter(Estimator != "gcomp") %>%
            mutate(Time = factor(Time),
                   Event = case_when(Event == 1 ~ "MACE", T ~ "Trt Disc")) %>%
            ggplot(aes(x = Time, y = `Pt Est`)) +
            theme_minimal() + facet_wrap(Estimand ~ Event, scales = "free", ncol = 2) +
            geom_point(aes(shape = Estimator, colour = Analysis), position = position_dodge(.5)) +
            geom_errorbar(aes(ymin = `CI Low`, ymax = `CI Hi`, linetype = Estimator, colour = Analysis),
                          position = position_dodge(.5)) +
            ylab("Risk Estimands") +
            labs(title = "SUSTAIN6 Prelim Analysis: Competing Risks vs. Right Censoring",
                 subtitle = "Effect of Sema on MACE & Trt Discontinuation")
        
        return(list("est" = est, "plot1" = comp_plot, "plot2" = comp_plot2,
                    "etc" = list(cr.est, event1.est, event2.est)))
    })
    names(out) <- analysis
    return(out)
}
Sus6Demo <- doSus6Demo(W = W, outcomes =  outcomes, analysis = c("AdvSideEffects"))


# MI ------------------------------------------------------------------------------------------
# tmp <- Sus6Demo[[2]]$est
# tmp %>% filter(Event == 1, Estimator != "gcomp") %>%
#     mutate(Time = as.factor(Time),
#            Analysis = factor(Analysis, levels = c("RightCens", "CompRisk"))) %>%
#     ggplot(aes(x = Time, y = `Pt Est`, colour = Analysis, shape = Estimator, linetype = Estimator)) +
#     geom_point(position = position_dodge2(.5)) +
#     geom_errorbar(aes(ymin = `CI Low`, ymax = `CI Hi`), position = position_dodge(.5)) +
#     theme_bw() + labs(title = "SUSTAIN-6: Effect of Sema on MI with Death as Censoring vs. Competing",
#                       y = "Risk Difference", x = "Days")


# MACE: ITT -----------------------------------------------------------------------------------
obs <- outcomes %>%
    dplyr::filter(PARAMCD == "TMTDMIST", ANL01FL == "Y") %>%
    dplyr::select(USUBJID, CNSR, AVAL) %>%
    left_join(W, ., by = "USUBJID") %>%
    as.data.table()
obs[, time := AVAL][, status := 1-CNSR]
obs[, trt := as.numeric(TR30PG2 == "Sema")]

targets <- 365.25*seq(0.5, 2, 0.5)
# KM ----

km <- survfit(Surv(time = time, event = status) ~ trt, data = obs) %>% summary(., times = targets) %>%
    with(., expr = cbind(Event = 1, A = rep(0:1, each = length(targets)), Risk = 1 - surv, se = std.err)) %>%
    as.data.table() %>% cbind(Time = rep(targets)) %>%
    .[, list("R1" = Risk[A == 1],
             "R1_se" = se[A == 1],
             "R0" = Risk[A == 0],
             "R0_se" = se[A == 0],
             "RD" = Risk[A == 1] -  Risk[A == 0],
             "RD_se" = sqrt(se[A == 1]^2 + se[A == 0]^2),
             "RR" = Risk[A == 1] / Risk[A == 0],
             "RR_se" = sqrt((se[A == 1] / Risk[A == 0])^2 +
                                se[A == 0]^2 * (Risk[A == 1] / Risk[A == 0]^2)^2)),
      by = c("Event", "Time")]  %>%
    melt(., id.vars = c("Event", "Time")) %>%
    .[, "Estimand" := stringr::str_extract(variable, "[A-Z0-9]+")] %>%
    .[, "type" := stringr::str_extract(variable, "[^(A-Z0-9_)]+")] %>%
    .[is.na(type), type := "Pt Est"] %>%
    .[, -"variable"] %>%
    dcast(., Event + Time + Estimand ~ type, value.var = "value") %>%
    .[, c("CI Low", "CI Hi", "Analysis", "Estimator") :=
          list(`Pt Est` - 1.96*se, `Pt Est` + 1.96*se, "IntenttoTrt", "NPMLE")]

# concrete ----
smallobs <- dplyr::select(obs, trt, time, status,
                          SEX, RACE, SMOKER, HBA1CBL, BMIBL, DIABDUR, EVCVDFL, INSTRTBS,
                          AGE, BMIBL, DIABDUR, EGFMDRD, HBA1CBL, SUMONTFL, SUFL, ETHNICFL,
                          EGFRTXT, PMYCINFL, PSSTIAFL, PRREVAFL, PSTN50FL, PSCOHDFL,
                          PASCISFL, PCKIDSFL, PMICPRFL, PHLVHYFL, PLVSDSFL, PABRINFL,
                          CV90DYFL, STRATUM, STRATUM2, STRATUM3, STRATUM4, STRATUM5,
                          CREAG1BL, ACRG1BL, MALBG1BL, ANTDBFL, EGFMDRD, EGFREPBC, CKDEPTXT,
                          RETISEVN, AGE, COUNTRY, SITEID) %>%
    mutate(ACRG1BL = case_when(ACRG1BL == "" ~ "NA",
                               T ~ ACRG1BL)) %>%
    # removing covariates with high missingness
    dplyr::select(-c("MALBG1BL", "PABRINFL", "PASCISFL", "PCKIDSFL", "PHLVHYFL", "PLVSDSFL",
                     "PMICPRFL", "PMYCINFL", "PRREVAFL", "PSCOHDFL", "PSSTIAFL", "PSTN50FL")) %>%
    # removing covariates with high collinearity
    dplyr::select(-c(STRATUM, STRATUM2, STRATUM3, STRATUM4, STRATUM5)) %>%
    dplyr::mutate_if(.predicate = is.character, ~ as.factor(.))
setDT(smallobs)
set.seed(0)
imp <- mice::mice(smallobs[, .SD, .SDcols = !c("trt", "time", "status")],
                  m = 1, maxit = 10, printFlag = FALSE)
smallobs[, BMIIMP := as.numeric(is.na(BMIBL))][, BMIBL := complete(imp)[["BMIBL"]]]

maceargs <- formatArguments(DataTable = smallobs,
                            EventTime = "time",
                            EventType = "status",
                            Treatment = "trt",
                            TargetTime = targets,
                            Intervention = 0:1)
maceargs$Model$`0`[["coxnet"]] <- maceargs$Model$`1`[["coxnet"]] <- "coxnet"
maceest <- doConcrete(ConcreteArgs = maceargs)

tmp <- as_tibble(getOutput(maceest)) %>%
    mutate(Estimand = case_when(Intervention == "[A=1] - [A=0]" ~ "RD",
                                Intervention == "[A=1] / [A=0]" ~ "RR",
                                Intervention == "A=1" ~ "R1",
                                Intervention == "A=0" ~ "R0",
                                T ~ "NA")) %>%
    subset(., select = c("Time", "Event", "Estimand",
                         "Estimator", "Pt Est", "se", "CI Low", "CI Hi")) %>%
    cbind("Analysis" = "IntenttoTrt", .) %>%
    rbind(km) %>%
    rbind(subset(Sus6Demo[[2]]$est,
                 select = c("Time", "Event", "Analysis", "Estimand",
                            "Estimator", "Pt Est", "se", "CI Low", "CI Hi")),
          .) %>%
    as.data.table()

tmp %>% filter(Event == 1, Estimator != "gcomp") %>%
    dplyr::filter(Estimand == "RR") %>%
    mutate(Time = factor(Time, labels = targets / 365.25 * 12,
                         levels = targets),
           Analysis = factor(Analysis, levels = c("IntenttoTrt", "RightCens", "CompRisk"),
                             labels = c("Intent to Treat", "Right Censoring", "Competing Risk")),
           Estimator = case_when(Estimator == "tmle" ~ "TMLE",
                                 Estimator == "NPMLE" ~ "KM/AJ",
                                 TRUE ~ Estimator)) %>%
    ggplot(aes(x = Time, y = `Pt Est`, colour = Analysis, shape = Estimator, linetype = Estimator)) +
    geom_point(position = position_dodge2(.6)) +
    # facet_wrap(~ Estimand, scales = "free") +
    geom_errorbar(aes(ymin = `CI Low`, ymax = `CI Hi`), position = position_dodge(.6)) +
    theme_bw() + labs(# title = "SUSTAIN-6: Effect of Sema on MACE in the presence of Treatment Discontinuation",
        y = "Risk Ratio", x = "Months")
# ggsave(filename = "/Shared/Projects/novo_nordisk/R/competing_risks/SUSTAIN6/RR-combined.png",
#        device = "png", path = "", width = 4, height = 3, units = "in")

# tmp %>% filter(Event == 2, Estimator != "gcomp") %>%
#     mutate(Time = as.factor(Time),
#            Analysis = factor(Analysis, levels = c("IntenttoTrt", "RightCens", "CompRisk"))) %>%
#     ggplot(aes(x = Time, y = `Pt Est`, colour = Analysis, shape = Estimator, linetype = Estimator)) +
#     geom_point(position = position_dodge2(.5)) + facet_wrap(~ Estimand, scales = "free") +
#     geom_errorbar(aes(ymin = `CI Low`, ymax = `CI Hi`), position = position_dodge(.5)) +
#     theme_bw() + labs(title = "SUSTAIN-6: Effect of Sema on Treatment Discontinuation in the Presence of MACE",
#                       y = "Risk Estimands", x = "Days")

tmp %>% filter(Estimator != "gcomp", Estimand %in% c("RD"),
               Analysis == "CompRisk")  %>%
    mutate(Time = factor(Time, labels = targets / 365.25 * 12,
                         levels = targets),
           Analysis = factor(Analysis, levels = c("IntenttoTrt", "RightCens", "CompRisk"),
                             labels = c("Intent to Treat", "Right Censoring", "Competing Risk")),
           Estimator = case_when(Estimator == "tmle" ~ "TMLE",
                                 Estimator == "NPMLE" ~ "KM/AJ",
                                 TRUE ~ Estimator)) %>%
    ggplot(aes(x = Time, y = `Pt Est`, colour = Estimator)) +
    geom_point(position = position_dodge2(.5)) + facet_wrap(Estimand ~ Event, scales = "free") +
    geom_errorbar(aes(ymin = `CI Low`, ymax = `CI Hi`), position = position_dodge(.5)) +
    theme_minimal() +
    labs(# title = "SUSTAIN-6: Effect of Sema on MACE in the presence of Treatment Discontinuation",
        y = "Risk Difference", x = "Months")

tmp[, RelEff :=  round(`se`[Estimator == "NPMLE"]^2 / `se`^2, 2),
    by = c("Analysis", "Time", "Event", "Estimand")]
tmp %>%
    .[Estimator != "gcomp" & Event == 1, ] %>% .[order(RelEff)] %>%
    dplyr::filter(Estimand == "R0") %>%
    dcast(Analysis + Time ~ Estimator, value.var = c("Pt Est", "RelEff")) %>%
    dplyr::select(-RelEff_NPMLE) %>%
    dplyr::mutate_at(.vars = c("Pt Est_NPMLE", "Pt Est_tmle"), ~round(1 - ., 2))
