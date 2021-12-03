

concr_tmle <- function(data, target_times, target_events, models, cv_args, ...) 
{
  # parameter checking --------------------------------------------------------------------------
  
  #
  # target time(s)
  # target event(s)
  # binary treatment
  # stopping criteria
  # ...
  
  
  # initial estimation / superlearning ----------------------------------------------------------
  
  ## cross validation setup ----
  # stratifying cv so that folds are balanced for treatment assignment & outcomes
  strata_ids <- as.numeric(factor(paste0(data$ARM, ":", data$EVENT)))
  cv_folds <- origami::make_folds(data, origami::folds_vfold, strata_ids = strata_ids)
  
  ## propensity score ----
  trt_task <- make_sl3_Task(
    data = data[, -c("TIME", "EVENT", "id")],
    covariates = colnames(data[, -c("TIME", "EVENT", "id", "ARM")]), 
    outcome = "ARM"
  )
  
  trt_sl <- Lrnr_sl$new(learners = models[["A"]], 
                        folds = cv_folds)
  trt_fit <- trt_sl$train(trt_task)
  propensity <- trt_fit$predict()
  
  if (min(propensity) < 0.05 | max(propensity) > 0.95) {
    warning("practical positivity violation likely, truncating propensity +-0.05\n")
    propensity[propensity < 0.05] <- 0.05
    propensity[propensity > 0.95] <- 0.95
  }
  
  ## hazards: events & censoring ----
  sl_models <- list()
  for (j in grep("\\d+", names(models), value = T)) {
    models_j <- models[[j]]
    library_risk <- data.table(matrix(NaN, nrow = nrow(data), ncol = length(models_j)))
    colnames(library_risk) <- names(models_j)
    
    for (fold_v in cv_folds) {
      train_indices <- fold_v$training_set
      val_indices <- fold_v$validation_set
      train_data <- data[train_indices, -c("id")]
      val_data <- data[val_indices, ][order(-TIME)]
      
      for (k in 1:length(models_j)) {
        ## train model ----
        coxph_args <- list("formula" = models_j[[k]], "data" = train_data)
        model_fit <- do.call(coxph, coxph_args)
        
        ## validation loss (-log partial likelihood) ----
        val_data[, fit.lp := predict(model_fit, type="lp", newdata=val_data)]
        val_data[, at.risk := cumsum(exp(fit.lp))]
        val_data[at.risk == 0, at.risk := 1]
        val_data[, names(models_j)[k] := (EVENT == j) * (fit.lp - log(at.risk))]
      }
      library_risk[val_indices, names(models_j) := subset(val_data, select = names(models_j))]
    }
    ## metalearner (discrete selector) ----
    sl_models[[j]] <- list("cv_risks" = -colSums(library_risk), 
                           "model" = models_j[[which.min(-colSums(library_risk))]])
  }
  
  ## fit sl selection ----
  
  sl_fits <- lapply(sl_models, function(slmod) {
    ## fit sl model ----
    model_fit <- do.call(coxph, list("formula" = slmod$model, "data" = data))
    
    haz <- rbind(data.table(time = 0, hazard = 0),
                 suppressWarnings(setDT(basehaz(model_fit, centered = TRUE))))
    colnames(haz) <- c("TIME", "haz")
    
    return(list("fit" = model_fit, "bhaz" = haz))
  })
  
  ## output initial estimates ----
  haz_times <- unique(c(target_times, data$TIME))
  haz_times <- haz_times[haz_times <= max(target_times)]
  
  hazards <- data.table("TIME" = c(0, haz_times))[order(TIME)]
  hazards <- rbind(cbind(hazards, A = 0), cbind(hazards, A = 1))
  
  ## baseline hazards ----
  for (j in grep("\\d+", names(models), value = T)) {
    hazards <- merge(hazards, sl_fits[[j]][["bhaz"]], by = "TIME", all.x = T)
    hazards[, haz := zoo::na.locf(haz), by = "A"]
    if (j == "0") {
      hazards[, haz := c(0, haz[-.N]), by = "A"]
      setnames(hazards, "haz", "cumhaz.C.t-")
    } else {
      # hazards[, paste0("cumbhaz.j", j) := haz]
      hazards[, haz := c(0, diff(haz)), by = "A"]
      setnames(hazards, "haz", paste0("bhaz.j", j))
    }
  }
  # important to order by id - clever covar A depends on this
  Ht <- rbind(data, copy(data)[, ARM := 1 - ARM])[order(id, ARM)] 
  setcolorder(Ht, c("id", "ARM", "TIME", "EVENT"))
  
  ## clever covariate A component ----
  Ht[ARM == 1, Ht_a := 1 / propensity]
  Ht[ARM == 0, Ht_a := 1 / (1 - propensity)]
  
  ## clever covariate hazard components ----
  
  for (j in grep("\\d+", names(models), value = T)) {
    if (j == "0") {
      Ht[, cox.C := predict(sl_fits[["0"]]$fit, newdata = Ht, type = "risk")]
    } else {
      Ht[, paste0("cox.j", j) := predict(sl_fits[[j]]$fit, 
                                         newdata = Ht, type = "risk")]
    }
  }
  setnames(Ht, "ARM", "A")
  setnames(Ht, "TIME", "T_tilde")
  Ht <- merge(data[, c("id", "ARM")], 
              Ht, 
              by = "id")
  Ht_colnames <- c("id", "ARM", "A", "T_tilde", "EVENT", "Ht_a",
                   grep("cox", colnames(Ht), value = T))
  Ht <- Ht[, ..Ht_colnames]
  
  Ht <- merge(Ht, hazards, by = "A", allow.cartesian = T)[order(id, TIME, A)]
  setcolorder(Ht, c("id", "TIME", "A"))
  
  ### clever covariate C component ----
  Ht[, "S.C.a.t-" := exp(-`cumhaz.C.t-` * cox.C)]
  if (min(Ht$`S.C.a.t-`) < 0.05) {
    warning("targeting a time where probability of being censored >95%, aim lower\n")
  }
  Ht[, "Ht.g" := Ht_a / `S.C.a.t-`]
  Ht <- Ht[, -c("Ht_a", "cox.C", "cumhaz.C.t-", "S.C.a.t-")]
  
  ### event-free survival ----
  Ht[, S.t := 1]
  for (j in target_events) {
    Ht[, S.t := S.t * exp(-cumsum(get(paste0("bhaz.j", j)) * 
                                    get(paste0("cox.j", j)))), 
       by = c("id", "A")]
  }
  Ht[, "S.t-" := c(1, S.t[-.N]), by = c("id", "A")]
  
  ### cause-specific risks ----
  for (j in target_events) {
    Ht[, paste0("F.j", j, ".t") := cumsum(S.t * get(paste0("bhaz.j", j)) * get(paste0("cox.j", j))), 
       by = c("id", "A")]
  }
  
  ### target event (j), target time (k), cause (l) clever covariate ----
  for (k in target_times) {
    for (j in target_events) {
      #### event j risk at target time k ----
      Ht[, fjtau := get(paste0("F.j", j, ".t"))[TIME == k], by = c("id", "A")]
      for (l in target_events) {
        Ht[, hjlt := 1] # why is this the default?
        #### event component ----
        Ht[S.t > 0, hjlt := (l ==j) - (fjtau - get(paste0("F.j", j, ".t"))) / S.t, 
           by = c("id", "A")]
        #### full clever covariate term ----
        Ht[S.t > 0, hjlt := hjlt * Ht.g * (TIME <= k) * (TIME <= T_tilde), 
           by = c("id", "A")]
        setnames(Ht, "hjlt", paste0("h.j", j, "l", l, ".t", k))
      }
    }
  }
  
  # Update step ---------------------------------------------------------------------------------
  
  # one-step tmle loop
  {
    ## EIC ----
    for (k in target_times) {
      for (j in target_events) {
        for (l in target_events) {
          Ht[, EIC := sum(get(paste0("h.j", j, "l", l, ".t", k)) * 
               ((TIME == T_tilde) * (EVENT == j) -
                  get(paste0("cox.j", j)) * get(paste0("bhaz.j", j)))), 
             by = c("id", "A")]
        }
        Ht[, EIC := sum(EIC)]
        setnames(Ht, "EIC", paste0("EIC.j", j, ".t", k))
      }
    }
    ## fluctuation
    
    ## update
  }
  
  
  
  
  # output --------------------------------------------------------------------------------------
  
  # g-comp (sl estimate)
  # unadjusted cox model
  # tmle & ic
  
  
}


