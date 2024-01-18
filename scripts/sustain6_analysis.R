setwd(dir = "/Shared/Projects/novo_nordisk/")
library(tidyverse); library(mice); library(data.table); library(survival); library(concrete)
library(SuperLearner); library(survtmle)

surv_tmle.dir <- c("/Shared/Projects/novo_nordisk/R/functions/")
x <- lapply(surv_tmle.dir, function(dir) lapply(list.files(dir, full.names = TRUE),
                                                function(x) try(source(x), silent = TRUE)))

glmnets <- create.Learner(base_learner = "SL.glmnet",
                          tune = list("alpha" = c(0, 0.5, 1)),
                          detailed_names = T)

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
  dplyr::select_if(~length(unique(.))>1) %>% 
  dplyr::mutate(SMOKER = factor(SMOKER, 
                                levels = c("Never smoked", "NK", "Previous smoker", "Current smoker"), 
                                labels = c("Never smoked", "NK", "Previous smoker", "Current smoker"))) %>% 
  dplyr::mutate_if(is.character, as.factor)

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
                       binwidth = c(1, 3, 6), seed = 123456789, 
                       cens_schemes = c("consistent", "early_cens"), 
                       boot = FALSE) {
  set.seed(seed)
  out <- lapply(analysis, function(definition) {
    obs <- outcomes %>%
      dplyr::filter(PARAMCD == "TMTDMIST", ANL02FL == "Y") %>%
      dplyr::select(USUBJID, CNSR, AVAL) %>%
      left_join(., W, by = "USUBJID")  %>% 
      mutate_if(is.character, ~factor(.)) %>%
      setDT()
    
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
    
    if (boot) {
      obs <- obs[sample(1:nrow(obs), nrow(obs), replace = TRUE), ]
    }
    
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
    getSurvtmleTbl <- function(survtmleOut, TargetEvent, TargetTime, Binwidth) {
      require(dplyr)
      x <- c("Risk", "se")
      survtmle.list <- lapply(seq_along(x), function(i) {
        survtmle.est <- cbind("Estimator" = paste0("survtmle-", Binwidth, "mo"),
                              "Intervention" = rep(c("A=0", "A=1"),
                                                   times = length(TargetEvent)),
                              "Event" = rep(TargetEvent, each = 2),
                              survtmleOut[[i]])
        survtmle.est <- melt(setDT(survtmle.est),
                             id.vars = c("Estimator", "Intervention", "Event"),
                             value.name = x[i], variable.name = "Time")
        survtmle.est[, Time := rep(TargetTime, each = length(TargetEvent) * 2)]
        return(survtmle.est)
      })
      out <- do.call(dplyr::full_join, survtmle.list)
      setDT(out)
      setnames(out, "Risk", "Pt Est")
    }
    
    screeners <- "All"
    
    sl_lib_g <- expand.grid(c("SL.glm", glmnets$names, "SL.bayesglm"), screeners)
    sl_lib_g <- lapply(1:nrow(sl_lib_g),
                       function(i) as.character(unlist(sl_lib_g[i,])))
    
    sl_lib_censor <- expand.grid(c("SL.glm", glmnets$names, "SL.bayesglm", "SL.xgboost"), screeners)
    sl_lib_censor <- lapply(1:nrow(sl_lib_censor),
                            function(i) as.character(unlist(sl_lib_censor[i,])))
    
    sl_lib_failure <- expand.grid(c("SL.glm", glmnets$names, "SL.bayesglm", "SL.xgboost"), screeners)
    sl_lib_failure <- lapply(1:nrow(sl_lib_failure),
                             function(i) as.character(unlist(sl_lib_failure[i,])))
    
    survtmle_out <- bind_rows(lapply(cens_schemes, function(cens_scheme) {
      survtmle_out <- lapply(binwidth, function(bin) {
        
        survtmleOut <- list("CompRisk" = smallobs$status, 
                            "J = 1" = smallobs$event1, 
                            "J = 2" = smallobs$event2)
        
        survtmleOut <- lapply(survtmleOut, function(y) {
          tmle_sl <- try(try-error, silent = TRUE)
          errors <- 0
          while(inherits(tmle_sl, "try-error") & errors < 10) {
            window <- 365.25 / 12 * bin
            
            z <- ceiling(smallobs$time / window)
            z_targ <- ceiling(targets / window)
            if (cens_scheme == "early_cens") {
              z[y != 0] <- z[y != 0] + 1
              z_targ <- z_targ + 1
            }
            
            tmle_sl <- try(survtmle(ftime = z,
                                    ftype = y,
                                    trt = smallobs$trt,
                                    adjustVars = smallobs[, !c("time", "status", "trt", "event1", "event2")],
                                    t0 = max(z_targ),
                                    SL.ftime = sl_lib_failure,
                                    SL.ctime = sl_lib_censor,
                                    SL.trt = sl_lib_g,
                                    returnIC = TRUE,
                                    returnModels = TRUE,
                                    method = "hazard",
                                    verbose = FALSE,
                                    maxIter = 100,
                                    Gcomp = FALSE))
            if (inherits(tmle_sl, "try-error"))
              errors <- errors + 1
          }
          if (!inherits(tmle_sl, "try-error")) {
            tmle_sl <- print(timepoints(tmle_sl, times = z_targ))
            tmle_sl <- getSurvtmleTbl(survtmleOut = tmle_sl, TargetEvent = c(1, 2), 
                                      TargetTime = targets, Binwidth = bin)
          }
          attr(tmle_sl, "errors") <- errors
          return(tmle_sl)
        })
      }) %>% 
        lapply(., function(x) {
          rbind(cbind("Analysis" = "CompRisk", x$CompRisk),
                cbind("Analysis" = "RightCens", x$`J = 1`[Event == 1, ]),
                cbind("Analysis" = "RightCens", x$`J = 2`[Event == 2, ]))
        }) %>% 
        bind_rows() %>% 
        .[, se := sqrt(se)]
      
      survtmle_out <- rbind(survtmle_out,
                            survtmle_out[, list("Intervention" = "RD",
                                                "Pt Est" = `Pt Est`[Intervention == "A=1"] - `Pt Est`[Intervention == "A=0"],
                                                "se" = sqrt(se[Intervention == "A=1"]^2 + se[Intervention == "A=0"]^2)),
                                         by = c("Analysis", "Estimator", "Event", "Time")], 
                            survtmle_out[, list("Intervention" = "RR",
                                                "Pt Est" = `Pt Est`[Intervention == "A=1"] / `Pt Est`[Intervention == "A=0"],
                                                "se" = sqrt((se[Intervention == "A=1"] / `Pt Est`[Intervention == "A=0"])^2 +
                                                              (se[Intervention == "A=0"] * `Pt Est`[Intervention == "A=1"] /
                                                                 `Pt Est`[Intervention == "A=0"]^2)^2)),
                                         by = c("Analysis", "Estimator", "Event", "Time")], 
                            fill = TRUE) %>% 
        mutate(Estimand = case_when(Intervention == "A=1" ~ "R1",
                                    Intervention == "A=0" ~ "R0",
                                    T ~ Intervention), 
               "CI Low" = `Pt Est` - 1.96*se, 
               "CI Hi"  = `Pt Est` + 1.96*se) %>% 
        dplyr::select(-Intervention)
      if (cens_scheme == "early_cens") 
        survtmle_out[, Estimator := paste0(Estimator, "-lagged")]
      return(survtmle_out)
    }))
    
    
    # concrete -----------------------------------------------------------------------------
    
    cr.args <- formatArguments(DataTable = smallobs[, !c("event1", "event2")],
                               EventTime = "time",
                               EventType = "status",
                               Treatment = "trt",
                               TargetTime = targets,
                               Intervention = 0:1)
    cr.args$Model[["0"]][["coxnet"]] <- cr.args$Model[["1"]][["coxnet"]] <-
      cr.args$Model[["2"]][["coxnet"]] <- "coxnet"
    formatArguments(cr.args)
    cr.est <- doConcrete(ConcreteArgs = cr.args)
    cr.out <- getOutput(ConcreteEst = cr.est, 
                        GComp = TRUE,
                        Estimand = c("RR", "RD", "Risk"),
                        Simultaneous = FALSE)
    
    event1.args <- formatArguments(DataTable = smallobs[, !c("status", "event2")],
                                   EventTime = "time",
                                   EventType = "event1",
                                   Treatment = "trt",
                                   TargetTime = targets,
                                   Intervention = 0:1)
    event1.args$Model[["0"]][["coxnet"]] <- event1.args$Model[["1"]][["coxnet"]] <- "coxnet"
    formatArguments(event1.args)
    event1.est <- doConcrete(event1.args)
    event1.out <- getOutput(ConcreteEst = event1.est,
                            GComp = TRUE,
                            Estimand = c("RR", "RD", "Risk"),
                            Simultaneous = FALSE)
    
    event2.args <- formatArguments(DataTable = smallobs[, !c("status", "event1")],
                                   EventTime = "time",
                                   EventType = "event2",
                                   Treatment = "trt",
                                   TargetTime = targets,
                                   Intervention = 0:1)
    event2.args$Model[["0"]][["coxnet"]] <- event2.args$Model[["2"]][["coxnet"]] <- "coxnet"
    formatArguments(event2.args)
    event2.est <- doConcrete(event2.args)
    event2.out <- getOutput(ConcreteEst = event2.est,
                            GComp = TRUE,
                            Estimand = c("RR", "RD", "Risk"),
                            Simultaneous = FALSE)
    
    # get output --------------------------------------------------------------
    
    est <- rbind(cbind("Analysis" = "CompRisk", cr.out),
                 cbind("Analysis" = "RightCens", event1.out),
                 cbind("Analysis" = "RightCens", event2.out)) %>%
      dplyr::filter(Event %in% 1:2)
    
    comp_plot <- est %>%
      mutate(Time = factor(Time),
             Event = case_when(Event == 1 ~ "MACE", T ~ "Trt Disc")) %>%
      ggplot(aes(x = Time, y = `Pt Est`, shape = Estimator, colour = Analysis, group = Analysis)) +
      theme_minimal() + facet_wrap(Intervention~Event, ncol = 2, scales = "free") +
      geom_point(position = position_dodge(.2)) + 
      geom_errorbar(aes(ymin = `CI Low`, ymax = `CI Hi`),
                    position = position_dodge(.2), width = 0.2) + 
      ylab("MACE Risk Ratio (TMLE)") +
      labs(title = "SUSTAIN6 Prelim Analysis: Competing Risks vs. Right Censoring",
           subtitle = "Effect of Sema on MACE & Trt Discontinuation")
    
    est <- est[order(Analysis, Estimator)] %>%
      mutate(Estimand = case_when(Intervention == "[A=1] - [A=0]" ~ "RD",
                                  Intervention == "[A=1] / [A=0]" ~ "RR",
                                  Intervention == "A=1" ~ "R1",
                                  Intervention == "A=0" ~ "R0",
                                  T ~ "NA")) %>%
      dplyr::select(Analysis, Estimand, Estimator, Time, Event, `Pt Est`,
                    se, `CI Low`, `CI Hi`) %>%
      rbind(., aj, km, survtmle_out)
    comp_plot2 <- est %>% filter(Estimator != "gcomp") %>%
      mutate(Time = factor(Time),
             Event = case_when(Event == 1 ~ "MACE", T ~ "Trt Disc")) %>%
      ggplot(aes(x = Time, y = `Pt Est`)) +
      theme_minimal() + facet_wrap(Estimand ~ Event, scales = "free", ncol = 2) +
      geom_point(aes(colour = Estimator, shape = Analysis), position = position_dodge2(.5)) +
      geom_errorbar(aes(ymin = `CI Low`, ymax = `CI Hi`, colour = Estimator, linetype = Analysis),
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

if (file.exists("/Shared/Projects/concrete/scripts/S6.RDS")) {
  Sus6Demo <- readRDS("/Shared/Projects/concrete/scripts/S6.RDS")
} else {
  Sus6Demo <- doSus6Demo(W = W, outcomes =  outcomes, analysis = c("AdvSideEffects"),
                         seed = 123456)
  saveRDS(Sus6Demo, file = "/Shared/Projects/concrete/scripts/S6.RDS")
}



# MACE: ITT -----------------------------------------------------------------------------------
doITT <- function(W, outcomes, targets = 365.25*seq(0.5, 2, 0.5),
                  binwidth = c(1, 3, 6), seed = 123456789, 
                  cens_schemes = c("consistent", "early_cens"), boot = FALSE) {
  set.seed(seed)
  obs <- outcomes %>%
    dplyr::filter(PARAMCD == "TMTDMIST", ANL01FL == "Y") %>%
    dplyr::select(USUBJID, CNSR, AVAL) %>%
    left_join(W, ., by = "USUBJID") %>%
    as.data.table()
  obs[, time := AVAL][, status := 1-CNSR]
  obs[, trt := as.numeric(TR30PG2 == "Sema")]
  
  if (boot) {
    obs <- obs[sample(1:nrow(obs), nrow(obs), replace = TRUE), ]
  }
  
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
  
  # set.seed(0)
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
  formatArguments(maceargs)
  maceest <- doConcrete(ConcreteArgs = maceargs)
  
  
  # survtmle ----------------------------------------------------------------
  getSurvtmleTbl <- function(survtmleOut, TargetEvent, TargetTime, Binwidth) {
    require(dplyr)
    x <- c("Risk", "se")
    survtmle.list <- lapply(seq_along(x), function(i) {
      survtmle.est <- cbind("Estimator" = paste0("survtmle-", Binwidth, "mo"),
                            "Intervention" = rep(c("A=0", "A=1"),
                                                 times = length(TargetEvent)),
                            "Event" = rep(TargetEvent, each = 2),
                            survtmleOut[[i]])
      survtmle.est <- melt(setDT(survtmle.est),
                           id.vars = c("Estimator", "Intervention", "Event"),
                           value.name = x[i], variable.name = "Time")
      survtmle.est[, Time := rep(TargetTime, each = length(TargetEvent) * 2)]
      return(survtmle.est)
    })
    out <- do.call(dplyr::full_join, survtmle.list)
    setDT(out)
    setnames(out, "Risk", "Pt Est")
  }
  
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
  survtmle_out <- bind_rows(lapply(cens_schemes, function(cens_scheme) {
    survtmle_out <- lapply(binwidth, function(bin) {
      survtmleOut <- list("IntenttoTrt" = smallobs$status)
      survtmleOut <- lapply(survtmleOut, function(y) {
        tmle_sl <- try(try-error, silent = TRUE)
        errors <- 0
        while(inherits(tmle_sl, "try-error") & errors < 10) {
          window <- 365.25 / 12 * bin
          
          z <- ceiling(smallobs$time / window)
          z_targ <- ceiling(targets / window)
          if (cens_scheme == "early_cens") {
            z[y != 0] <- z[y != 0] + 1
            z_targ <- z_targ + 1
          }
          
          tmle_sl <- try(survtmle(ftime = z,
                                  ftype = y,
                                  trt = smallobs$trt,
                                  adjustVars = smallobs[, !c("time", "status", "trt", "event1", "event2")],
                                  t0 = max(z_targ),
                                  SL.ftime = sl_lib_failure,
                                  SL.ctime = sl_lib_censor,
                                  SL.trt = sl_lib_g,
                                  returnIC = TRUE,
                                  returnModels = TRUE,
                                  method = "hazard",
                                  verbose = FALSE,
                                  maxIter = 100,
                                  Gcomp = FALSE))
          if (inherits(tmle_sl, "try-error")) errors <- errors + 1
        }
        if (!inherits(tmle_sl, "try-error")) {
          tmle_sl <- print(timepoints(tmle_sl, times = z_targ))
          tmle_sl <- getSurvtmleTbl(survtmleOut = tmle_sl, TargetEvent = c(1, 2), 
                                    TargetTime = targets, Binwidth = bin)
        }
        attr(tmle_sl, "errors") <- errors
        return(tmle_sl)
      })
    })
    survtmle_out <- bind_rows(lapply(survtmle_out, function(x) {
      cbind("Analysis" = "IntenttoTrt", "Errors" = attr(x, "errors"), x$IntenttoTrt)
    }))
    survtmle_out[, se := sqrt(se)]
    survtmle_out <- rbind(survtmle_out,
                          survtmle_out[, list("Intervention" = "RD",
                                              "Pt Est" = `Pt Est`[Intervention == "A=1"] - `Pt Est`[Intervention == "A=0"],
                                              "se" = sqrt(se[Intervention == "A=1"]^2 + se[Intervention == "A=0"]^2)),
                                       by = c("Analysis", "Estimator", "Event", "Time")], 
                          survtmle_out[, list("Intervention" = "RR",
                                              "Pt Est" = `Pt Est`[Intervention == "A=1"] / `Pt Est`[Intervention == "A=0"],
                                              "se" = sqrt((se[Intervention == "A=1"] / `Pt Est`[Intervention == "A=0"])^2 +
                                                            (se[Intervention == "A=0"] * `Pt Est`[Intervention == "A=1"] /
                                                               `Pt Est`[Intervention == "A=0"]^2)^2)),
                                       by = c("Analysis", "Estimator", "Event", "Time")], 
                          fill = TRUE)
    survtmle_out <- survtmle_out %>% 
      mutate(Estimand = case_when(Intervention == "A=1" ~ "R1",
                                  Intervention == "A=0" ~ "R0",
                                  T ~ Intervention), 
             "CI Low" = `Pt Est` - 1.96*se, 
             "CI Hi"  = `Pt Est` + 1.96*se) %>% 
      dplyr::select(-Intervention)
    if (cens_scheme == "early_cens") 
      survtmle_out[, Estimator := paste0(Estimator, "-lagged")]
    return(survtmle_out)
  }))
  output <- list("maceest" = maceest, "survtmle_out" = survtmle_out, "km" = km)
  return(output)
}


# combine -----------------------------------------------------------------
if (file.exists("/Shared/Projects/concrete/scripts/S6-combined.RDS")) {
  tmp <- readRDS("/Shared/Projects/concrete/scripts/S6-combined.RDS")
} else {
  ITT_out <- doITT(W = W, outcomes =  outcomes)
  tmp <- as_tibble(getOutput(ITT_out$maceest, Estimand = c("RR", "RD", "Risk"), 
                             Simultaneous = FALSE)) %>%
    mutate(Estimand = case_when(Intervention == "[A=1] - [A=0]" ~ "RD",
                                Intervention == "[A=1] / [A=0]" ~ "RR",
                                Intervention == "A=1" ~ "R1",
                                Intervention == "A=0" ~ "R0",
                                T ~ "NA")) %>%
    subset(., select = c("Time", "Event", "Estimand",
                         "Estimator", "Pt Est", "se", "CI Low", "CI Hi")) %>%
    cbind("Analysis" = "IntenttoTrt", .) %>%
    rbind(., ITT_out$km, ITT_out$survtmle_out) %>%
    rbind(subset(Sus6Demo[["AdvSideEffects"]]$est,
                 select = c("Time", "Event", "Analysis", "Estimand",
                            "Estimator", "Pt Est", "se", "CI Low", "CI Hi")),
          .) %>%
    as.data.table() %>% 
    distinct()
  saveRDS(tmp, file = "/Shared/Projects/concrete/scripts/S6-combined.RDS")
}




# plot --------------------------------------------------------------------

estimands <- c("RR" = "Risk Ratio", "RD" = "Risk Difference", 
               "R0" = "Placebo Risk", "R1" = "Treated Risk")
plots <- lapply(seq_along(estimands), function(i) {
  tmp %>% filter(Event == 1, Estimator != "gcomp", Estimand == names(estimands)[i]) %>%
    # dplyr::filter(Estimand == "RR") %>%
    mutate(Time = factor(Time, labels = sort(unique(tmp$Time)) / 365.25 * 12,
                         levels = sort(unique(tmp$Time))),
           Analysis = factor(Analysis, levels = c("IntenttoTrt", "RightCens", "CompRisk"),
                             labels = c("Intent to Treat", "Right Censoring", "Competing Risk")),
           Estimator = case_when(Estimator == "tmle" ~ "Cont-TMLE",
                                 Estimator == "NPMLE" ~ "K-M / A-J",
                                 Estimator == "survtmle-1mo" ~ "Disc-TMLE 1mo", 
                                 Estimator == "survtmle-3mo" ~ "Disc-TMLE 3mo", 
                                 Estimator == "survtmle-6mo" ~ "Disc-TMLE 6mo", 
                                 Estimator == "survtmle-1mo-lagged" ~ "Disc-TMLE 1mo -", 
                                 Estimator == "survtmle-3mo-lagged" ~ "Disc-TMLE 3mo -", 
                                 Estimator == "survtmle-6mo-lagged" ~ "Disc-TMLE 6mo -"), 
           Estimator = factor(Estimator, levels = c("K-M / A-J", "Cont-TMLE", "Disc-TMLE 1mo", 
                                                    "Disc-TMLE 3mo", "Disc-TMLE 6mo", 
                                                    "Disc-TMLE 1mo -", "Disc-TMLE 3mo -", "Disc-TMLE 6mo -"))) %>%
    ggplot(aes(x = Time, y = `Pt Est`, color = Estimator)) +
    geom_point(position = position_dodge2(width = .5)) +
    facet_wrap(~ Analysis, scales = "free", ncol = 1) +
    geom_errorbar(aes(ymin = `CI Low`, ymax = `CI Hi`), position = position_dodge(width = .5)) +
    theme_bw() + labs(# title = "SUSTAIN-6: Effect of Sema on MACE in the presence of Treatment Discontinuation",
      y = estimands[i], x = "Months")
})

plots_combined <- 
  tmp %>% filter(Event == 1, Estimator != "gcomp") %>%
  mutate(Time = factor(Time, labels = sort(unique(tmp$Time)) / 365.25 * 12,
                       levels = sort(unique(tmp$Time))),
         Analysis = factor(Analysis, levels = c("IntenttoTrt", "RightCens", "CompRisk"),
                           labels = c("Intent to Treat", "Right Censoring", "Competing Risk")),
         Estimator = case_when(Estimator == "tmle" ~ "Cont-TMLE",
                               Estimator == "NPMLE" ~ "K-M / A-J",
                               Estimator == "survtmle-1mo" ~ "Disc-TMLE 1mo", 
                               Estimator == "survtmle-3mo" ~ "Disc-TMLE 3mo", 
                               Estimator == "survtmle-6mo" ~ "Disc-TMLE 6mo"), 
         Estimator = factor(Estimator, levels = c("K-M / A-J", "Cont-TMLE", "Disc-TMLE 1mo", 
                                                  "Disc-TMLE 3mo", "Disc-TMLE 6mo"))) %>% 
  dplyr::filter(Estimator != "NA") %>% 
  ggplot(aes(x = Time, y = `Pt Est`, color = Estimator, shape = Analysis, linetype = Analysis)) +
  geom_point(position = position_dodge2(width = .8)) +
  facet_wrap(~ Estimand, scales = "free", ncol = 1) +
  geom_errorbar(aes(ymin = `CI Low`, ymax = `CI Hi`), position = position_dodge(width = .8)) +
  theme_bw() + labs(y = "Estimates", x = "Months")

releff_tbl <- tmp[, RelEff :=  round(se[Estimator == "NPMLE"]^2 / se^2, 2),
                  by = c("Analysis", "Time", "Event", "Estimand")] %>% 
  mutate(Analysis = factor(Analysis, levels = c("IntenttoTrt", "RightCens", "CompRisk"),
                           labels = c("Intent to Treat", "Right Censoring", "Competing Risk"))) %>% 
  .[Estimator != "gcomp" & Event == 1, ] %>% .[order(RelEff)] %>%
  # dplyr::filter(Estimand %in% c("R0", "R1")) %>% 
  dplyr::select(Time, Analysis, Estimand, Estimator, `Pt Est`, RelEff) %>% 
  arrange(Estimator, Analysis, Estimand, Time) %>% 
  dcast(Estimator + Analysis + Time ~ Estimand, value.var = "RelEff")
releff_tbl %>% dcast(Analysis + Time ~ Estimator, value.var = "RD")


ptest_tbl <- tmp[, `Pt Est` :=  round(`Pt Est`, 2),
                 by = c("Analysis", "Time", "Event", "Estimand")] %>% 
  mutate(Analysis = factor(Analysis, levels = c("IntenttoTrt", "RightCens", "CompRisk"),
                           labels = c("Intent to Treat", "Right Censoring", "Competing Risk"))) %>% 
  .[Estimator != "gcomp" & Event == 1, ] %>% .[order(RelEff)] %>%
  # dplyr::filter(Estimand %in% c("R0", "R1")) %>% 
  dplyr::select(Time, Analysis, Estimand, Estimator, `Pt Est`, RelEff) %>% 
  arrange(Estimator, Analysis, Estimand, Time) %>% 
  dcast(Estimator + Analysis + Time ~ Estimand, value.var = "Pt Est")
ptest_tbl %>% dcast(Analysis + Time ~ Estimator, value.var = "RR")

tmp[Estimator == "NPMLE" & Estimand == "RR" & Event == 1, ] %>% 
  mutate(`CI Low` = round(`CI Low`, 2), 
         `CI Hi` = round(`CI Hi`, 2))

tmp[, RelEff :=  round(se[Estimator == "NPMLE"]^2 / se^2, 2),
    by = c("Analysis", "Time", "Event", "Estimand")] %>% 
  mutate(Output = paste0(round(`Pt Est`, 2), " (", round(`CI Low`, 2), " - ", 
                         round(`CI Hi`, 2), ") [", round(RelEff, 2), "]")) %>% 
  dcast(Analysis + Time + Estimand + Event ~ Estimator, value.var = "Output")

tmp %>% mutate(Time = factor(Time)) %>% 
  dplyr::filter(Event == 1) %>% 
  ggplot(aes(x = Time, y = RR, ymin = `CI Low`, ymax = `CI Hi`, colour = Estimator)) + 
  geom_point(position = position_dodge2(0.3)) + 
  geom_errorbar(position = position_dodge2(0.3), width = 0.3) + 
  theme_minimal() + 
  facet_wrap(Analysis ~ Event, scales = "free", ncol = 1)

tmp[, RelEff :=  round(se[Estimator == "NPMLE"]^2 / se^2, 2),
    by = c("Analysis", "Time", "Event", "Estimand")] %>% 
  mutate(Output = paste0(round(`Pt Est`, 2), " (", round(`CI Low`, 2), " - ", 
                         round(`CI Hi`, 2), ") [", round(RelEff, 2), "]")) %>% 
  dplyr::filter(Event == 1, Estimator %in% c("NPMLE", "tmle", "survtmle-1mo", "survtmle-3mo", "survtmle-6mo"), 
                Estimand == "RR", Analysis == "CompRisk") %>% 
  dcast(Analysis + Estimand + Time ~ Estimator, value.var = "Output") %>% 
  dplyr::select(NPMLE, tmle, `survtmle-1mo`, `survtmle-3mo`, `survtmle-6mo`) %>% 
  knitr::kable(format = "latex")


# log-transformed delta method se -----------------------------------------

tmp %>% 
  dplyr::filter(Estimator != "gcomp") %>% 
  pivot_wider(id_cols = c("Time", "Event", "Analysis", "Estimator"), 
              names_from = "Estimand", values_from = c("Pt Est", "se")) %>% 
  group_by(Time, Event, Analysis, Estimator) %>% 
  transmute(
    rr = `Pt Est_RR`, 
    logse = sqrt( (se_R1 / `Pt Est_R1`)^2 + (se_R0 / `Pt Est_R0`)^2 ), 
    ci_ll = exp(log(`Pt Est_RR`) - 1.96*logse), 
    ci_ld = `Pt Est_RR` - 1.96*se_RR,
    ci_ul = exp(log(`Pt Est_RR`) + 1.96*logse), 
    ci_ud = `Pt Est_RR` + 1.96*se_RR,
    ci_l = ci_ul - ci_ll,
    ci_d = ci_ud - ci_ld
  ) %>% 
  group_by(Time, Event, Analysis) %>% 
  arrange(Event, Analysis, Time) %>% 
  mutate(relci_l = ci_l / head(ci_l, 1),
         relci_d = ci_d / head(ci_d, 1), 
         rel_ci_ld = (ci_l / ci_d)) %>% 
  dplyr::filter(Event == 1, 
                !grepl("lagged", Estimator)) %>% 
  pivot_longer(cols = c("relci_l", "relci_d", "rel_ci_ld")) %>% 
  mutate(name = case_when(name == "relci_l" ~ "Standardized Log-transformed CI widths",
                          name == "relci_d" ~ "Standardized Direct CI widths",
                          name == "rel_ci_ld" ~ "Relative Log-transformed/Direct CI widths")) %>% 
  dplyr::filter(name == "Relative Log-transformed/Direct CI widths") %>% 
  ggplot(aes(x= Time, y = value, colour = Estimator)) + 
  theme_minimal() + facet_wrap(Analysis~name, scales = "free", nrow = 1) + 
  geom_line()


# bootstrap ---------------------------------------------------------------
# bootstrap
library(doParallel); library(foreach)
cl <- makeForkCluster(nnodes = 10)
registerDoParallel(cl)
set.seed(12345678)
seeds <- sample(0:1e9, size = 100, replace = FALSE)
foreach(i = seq_along(seeds), 
        .packages = c("SuperLearner", "tidyverse", "data.table")) %dopar% {
          S6_boot <- doSus6Demo(W = W, outcomes =  outcomes,
                                analysis = c("AdvSideEffects"),
                                cens_schemes = "consistent",
                                seed = seeds[i], boot = TRUE)
          ITT_boot <- doITT(W = W, outcomes =  outcomes, seed = seeds[i], 
                            cens_schemes = "consistent", boot = TRUE)
          boot_out <- as_tibble(getOutput(ITT_boot$maceest, Estimand = c("RR", "RD", "Risk"), 
                                          Simultaneous = FALSE)) %>%
            mutate(Estimand = case_when(Intervention == "[A=1] - [A=0]" ~ "RD",
                                        Intervention == "[A=1] / [A=0]" ~ "RR",
                                        Intervention == "A=1" ~ "R1",
                                        Intervention == "A=0" ~ "R0",
                                        T ~ "NA")) %>%
            subset(., select = c("Time", "Event", "Estimand",
                                 "Estimator", "Pt Est", "se", "CI Low", "CI Hi")) %>%
            cbind("Analysis" = "IntenttoTrt", .) %>%
            rbind(., ITT_boot$km, ITT_boot$survtmle_out) %>%
            rbind(subset(S6_boot[["AdvSideEffects"]]$est,
                         select = c("Time", "Event", "Analysis", "Estimand",
                                    "Estimator", "Pt Est", "se", "CI Low", "CI Hi")),
                  .) %>%
            as.data.table() %>% 
            distinct()
          saveRDS(object = boot_out, file = paste0("/Shared/Projects/concrete/scripts/S6b", i, ".RDS"))
          NULL
        }

boot <- lapply(seq_along(seeds), function(i) {
  readRDS(file = paste0("/Shared/Projects/concrete/scripts/S6b", i, ".RDS"))
}) %>% bind_rows() %>% 
  mutate(coverage = `Pt Est` <= `CI Hi` & `Pt Est` >= `CI Low`) %>% 
  group_by(Time, Event, Analysis, Estimand, Estimator) %>% 
  summarise(PtEst = mean(`Pt Est`), 
            se_boot = sqrt(mean((`Pt Est` - mean(`Pt Est`))^2))) %>% 
  left_join(., 
            dplyr::select(tmp, Time, Event, Analysis, Estimand, Estimator, se)) %>% 
  arrange(Event, Analysis, Estimand, Time, Estimator) %>% 
  group_by(Time, Event, Analysis, Estimand) %>% 
  mutate(rel_ci = se / head(se, 1), 
         rel_ciboot = se_boot / head(se_boot, 1), 
         relse_ic_boot = se / se_boot)

boot %>% 
  dplyr::filter(Event == 1, Estimator != "gcomp") %>% 
  pivot_longer(cols = c("relse_ic_boot")) %>% 
  ggplot(aes(x = Time, y = value, colour = Estimator)) + 
  theme_minimal() + facet_wrap(Analysis ~ Estimand, scales = "free") + 
  geom_line() + 
  labs(y = "Relative SE (IC / Boot)")
