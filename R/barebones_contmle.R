
concr_tmle <- function(data, target_times, target_events, models, cv_args = NULL, 
                       no.small.steps = 30, onestep_eps = 0.5, verbose = FALSE, ...) 
{
  # parameter checking --------------------------------------------------------------------------
  
  events <- unique(data$EVENT)
  events <- sort(events[events > 0])
  
  PnEIC_norm_fun <- function(x, y) {
    return(sqrt(sum(unlist(x) * unlist(y))))
  }
  
  PnEIC_wt_fun <- function(PnEIC) {
    return(PnEIC)
  }
  
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
      Ht[, (paste0("cox.j", j)) := predict(sl_fits[[j]]$fit, 
                                         newdata = Ht, type = "risk")]
    }
  }
  setnames(Ht, "ARM", "A") # A is the target trts, ARM is the assigned
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
    warning(paste0("targeting a time where probability of being censored >95%,",
                   "cens survival truncated at min 0.05, but you should really", 
                   "aim lower\n"))
  }
  Ht[, "Ht.g" := Ht_a / `S.C.a.t-`]
  Ht <- Ht[, -c("Ht_a", "cox.C", "cumhaz.C.t-", "S.C.a.t-")]
  
  ### event-free survival ----
  Ht[, S.t := 1]
  
  ## expand out bhaz and cox to js, then do by = c("id", "A", "J")
  
  for (j in target_events) {
    Ht[, S.t := S.t * exp(-cumsum(get(paste0("bhaz.j", j)) * 
                                    get(paste0("cox.j", j)))), 
       by = c("id", "A")]
  }
  if (min(Ht[, S.t]) <= 0) {
    stop("you're targeting a time when there are no survivors, aim lower\n")
  }
  Ht[, "S.t-" := c(1, S.t[-.N]), by = c("id", "A")]
  
  ### cause-specific risks ----
  for (j in target_events) {
    Ht[, (paste0("F.j", j, ".t")) := cumsum(`S.t-` * 
                                            get(paste0("bhaz.j", j)) * 
                                            get(paste0("cox.j", j))
    ), 
    by = c("id", "A")]
  }
  
  
  ### cause-specific risks at target times ----
  ### target event (j), target time (k), cause (l) clever covariate ----
  for (k in target_times) {
    Ht[, Ht.g := Ht.g * (TIME <= k) * (TIME <= T_tilde), (ARM == A)]
    for (j in target_events) {
      #### event j risk at target time k ----
      Ht[, fjtau := get(paste0("F.j", j, ".t"))[TIME == k], by = c("id", "A")]
      for (l in events) {
        Ht[, hjlt := 1] # why was 1 the default, a 0 should squash oversights no?
        #### event component ----
        Ht[, hjlt := (l == j) - (fjtau - get(paste0("F.j", j, ".t"))) / S.t]
        setnames(Ht, "hjlt", paste0("h.j", j, ".l", l, ".t", k))
      } #### mean treated/control event j risk at target time k ----
      Ht[, paste0("Psi.j", j, ".t", k) := mean(fjtau), by = "A"]
      setnames(Ht, "fjtau", paste0("F.j", j, ".t", k))
    }
  }
  
  ## initial estimator (g-computation) ----
  Psi_init <- Ht[, mget(c("A", grep("Psi.+", colnames(Ht), value = T)))]
  Psi_init <- Psi_init[, lapply(.SD, unique), by = "A"]
  Psi_init <- melt(Psi_init, id.vars = "A")
  Psi_init[, c("dummy", "event", "time") := tstrsplit(variable, "\\..")]
  Psi_init <- Psi_init[, -c("variable", "dummy")]
  Psi_init <- dcast(Psi_init, A + time ~ event, value.var = "value")
  Psi_init[, S := do.call(sum, mget(paste0(1:3))), by = c("A", "time")]
  Psi_init <- Psi_init[, lapply(.SD, as.numeric)][order(time, A)]
  
  
  # Update step ---------------------------------------------------------------------------------
  
  ## EIC ----
  Ht_eic <- data.table(id = data$id)
  for (k in target_times) {
    for (j in target_events) {
      Ht[, EIC := 0]
      for (l in events) {
        Ht[, EIC := EIC + get(paste0("h.j", j, ".l", l, ".t", k)) * Ht.g * 
             ((TIME == T_tilde ) * (EVENT == l) - 
                get(paste0("cox.j", l)) * get(paste0("bhaz.j", l)))]
      }
      Ht_eic <- merge(Ht_eic, 
                      dcast(Ht[, .(EIC = sum(EIC), 
                                   Fjt = unique(get(paste0("F.j", j, ".t", k))), 
                                   Psijt = unique(get(paste0("Psi.j", j, ".t", k)))), 
                               by = c("id", "A")][, EIC := EIC + Fjt - Psijt], 
                            id ~ A, value.var = "EIC"), 
                      by = "id")
      setnames(Ht_eic, "0", paste0("EIC.a0.j", j, ".t", k))
      setnames(Ht_eic, "1", paste0("EIC.a1.j", j, ".t", k))
    }
  }
  
  init_ic <- list("ic" = Ht_eic[, -c("id")], 
                  "pnEIC" = colMeans(Ht_eic[, -c("id")]),
                  "seEIC" = sqrt(diag(var(Ht_eic[, -c("id")]))))
  
  PnEIC <- init_ic$pnEIC
  PnEIC <- as.data.table(cbind("A" = 1:0,
                               rbind(PnEIC[grep("a1", names(PnEIC))], 
                                     PnEIC[grep("a0", names(PnEIC))])))
  eica <- do.call(paste0, expand.grid("PnEIC.j", target_events, ".t", target_times))
  setnames(PnEIC, 2:ncol(PnEIC), eica)
  
  PnEIC_wtd <- PnEIC_wt_fun(PnEIC)
  
  Ht <- merge(Ht, PnEIC_wtd, by = "A")
  PnEIC_norm <- PnEIC_norm_fun(PnEIC, PnEIC_wtd)
  
  ## one-step tmle loop (one-step) ----
  
  ## 5.2 set update epsilon down-scaling factor --------------------------------
  gc()
  
  for (step in 1:no.small.steps) {
    cat("starting step", step, "with update epsilon =", onestep_eps, "\n")
    PnEIC_prev <- PnEIC
    PnEIC_wtd_prev <- PnEIC_wtd
    PnEIC_norm_prev <- PnEIC_norm
    
    ## make backups ----
    ## save current cause-specific hazards and clever covs in case the update
    ## makes things worse
    for (j in target_events) { # loop over target events
      for (k in target_times) { # loop over target times
        for (l in events) { # nested loop over events
          Ht[, (paste0("h.j", j, ".l", l, ".t", k, ".tmp")) :=
               get(paste0("h.j", j, ".l", l, ".t", k))]
        }
      }
      Ht[, (paste0("cox.j", j, ".tmp")) := get(paste0("cox.j", j))]
    }
    
    ## 5.1 calculate update step direction -------------------------------------
    
    for (l in events) { # loop over event hazards
      Ht[, (paste0("delta.j", l, ".dx")) := 0] # fluctuation of causes-specific hazard
      for (j in target_events) { # loop over target events
        for (k in target_times) { # loop over target times
          Ht[, (paste0("delta.j", l, ".dx")) :=
               get(paste0("delta.j", l, ".dx")) + (TIME <= k) * Ht.g *
               get(paste0("h.j", j, ".l", l, ".t", k)) * 
               get(paste0("PnEIC.j", j, ".t", k)) / PnEIC_norm]
        }
      }
    }
    
    ## 5.3 update cause specific hazards ----------------------------------------
    Ht[, S.t := 1]
    for (l in events) { # loop over competing events
      Ht[, (paste0("cox.j", l)) := get(paste0("cox.j", l)) *
           exp(onestep_eps * get(paste0("delta.j", l, ".dx")))]
      #### truncation cox fit +- 500 why ? -----
      Ht[, (paste0("cox.j", l)) := pmax(-500, pmin(500, get(paste0("cox.j", l))))]
      # Does changing the order of event updates matter? -------------------------------
      Ht[, S.t := S.t * exp(-cumsum(get(paste0("bhaz.j", l)) * 
                                      get(paste0("cox.j", l)))), 
         by = c("id", "A")]
      # when does cox go to Infinity? -------
      if (any(is.infinite(Ht[, get(paste0("cox.j", l))])))
        stop(paste0("cox.fit for j=", j, " went to infinity. wtf happened?\n"))
      if (min(Ht[, S.t]) <= 0)
        stop(paste0("Survival at some time is leq 0. wtf happened?\n"))
    }
    Ht[, "S.t-" := c(1, S.t[-.N]), by = c("id", "A")]
    
    ## update clever covariates ----
    for (j in target_events) {
      Ht[, (paste0("F.j", j, ".t")) := cumsum(`S.t-` * 
                                              get(paste0("bhaz.j", j)) * 
                                              get(paste0("cox.j", j))), 
         by = c("id", "A")]
      for (k in target_times) {
        # Ht[, (paste0("S.t", k)) := S.t[TIME == k], by = c("id", "A")]
        Ht[, fjtau := get(paste0("F.j", j, ".t"))[TIME == k], by = c("id", "A")]
        for (l in events) {
          Ht[, hjlt := 1] #### ASK HELENE: why is 1 the default? ---------
          #### event component ----
          Ht[, hjlt := (l == j) - (fjtau - get(paste0("F.j", j, ".t"))) / S.t]
          
          Ht[round(S.t, 8) == 0, hjlt := 0 - (l == j)] ##### ASK HELENE WHY -------
          
          Ht[, (paste0("h.j", j, ".l", l, ".t", k)) := hjlt]
        } #### mean treated/control event j risk at target time k ----
        Ht[, (paste0("Psi.j", j, ".t", k)) := mean(fjtau), by = "A"]
        Ht[, (paste0("F.j", j, ".t", k)) := fjtau]
      }
    }
    Ht <- Ht[, -c("hjlt", "fjtau")]
    ## recalculate eic -----
    Ht_eic <- data.table(id = data$id)
    for (k in target_times) {
      for (j in target_events) {
        Ht[, EIC := 0]
        for (l in events) {
          Ht[, EIC := EIC + get(paste0("h.j", j, ".l", l, ".t", k)) * Ht.g * 
               ((TIME == T_tilde ) * (EVENT == l) - 
                  get(paste0("cox.j", l)) * get(paste0("bhaz.j", l)))]
        }
        Ht_eic <- merge(Ht_eic, 
                        dcast(Ht[, .(EIC = sum(EIC), 
                                     Fjt = unique(get(paste0("F.j", j, ".t", k))), 
                                     Psijt = unique(get(paste0("Psi.j", j, ".t", k)))), 
                                 by = c("id", "A")][, EIC := EIC + Fjt - Psijt], 
                              id ~ A, value.var = "EIC"), 
                        by = "id")
        setnames(Ht_eic, "0", paste0("EIC.a0.j", j, ".t", k))
        setnames(Ht_eic, "1", paste0("EIC.a1.j", j, ".t", k))
      }
    }
    
    ic <- list("ic" = Ht_eic, 
               "pnEIC" = colMeans(Ht_eic[, -c("id")]),
               "seEIC" = sqrt(diag(var(Ht_eic[, -c("id")]))))
    
    PnEIC <- ic$pnEIC
    PnEIC <- as.data.table(cbind("A" = 1:0,
                                 rbind(PnEIC[grep("a1", names(PnEIC))], 
                                       PnEIC[grep("a0", names(PnEIC))])))
    eica <- do.call(paste0, expand.grid("PnEIC.j", target_events, ".t", target_times))
    setnames(PnEIC, 2:ncol(PnEIC), eica)
    
    PnEIC_wtd <- PnEIC_wt_fun(PnEIC)
    PnEIC_norm <- PnEIC_norm_fun(PnEIC, PnEIC_wtd)
    
    if (PnEIC_norm_prev <= PnEIC_norm) {
      warning(paste0("update overshot! one-step update epsilon of ", 
                     onestep_eps, " will be halved\n"))
      onestep_eps <- onestep_eps / 2
      ## Revert to previous cshaz & clevcovs if update made things worse
      for (j in target_events) { # loop over target events
        for (k in target_times) { # loop over target times
          for (l in events) { # nested loop over events
            Ht[, (paste0("h.j", j, ".l", l, ".t", k)) :=
                 get(paste0("h.j", j, ".l", l, ".t", k, ".tmp"))]
          }
        }
        Ht[, (paste0("cox.j", j)) := get(paste0("cox.j", j, ".tmp"))]
      }
      
      PnEIC <- PnEIC_prev
      PnEIC_wtd <- PnEIC_wtd_prev
      PnEIC_norm <- PnEIC_norm_prev
    } else {
      # update PnEIC columns in Ht with updated values
      ht_eic_cols <- grep("PnEIC", colnames(Ht), value = T)
      Ht <- Ht[, !..ht_eic_cols]
      Ht[PnEIC, on = .(A = A), (ht_eic_cols) := mget(sprintf("i.%s", ht_eic_cols))]
      gc() # otherwise memory blows up
      
      summ_eic <- melt(ic$ic, id.vars = "id")
      summ_eic[, c("EIC", "A", "J", "T") := tstrsplit(variable, "\\.")]
      summ_eic <- summ_eic[, -c("variable", "EIC")]
      summ_eic <- dcast(summ_eic, id + A + `T` ~ J, value.var = "value")[
        , S := do.call(sum, mget(paste0("j", 1:3))), 
        by = c("id", "A", "T")]
      summ_eic <- dcast(summ_eic, id ~ ..., 
                        value.var = c("j1", "j2", "j3", "S"), 
                        sep = ".")
      
      onestep_stop <- abs(colMeans(summ_eic[, -"id"])) <= 
        sqrt(colMeans(summ_eic[, -"id"]^2)) / ( sqrt(nrow(data)) * log(nrow(data)))
      
      if (all(onestep_stop) | step == no.small.steps) {
        if (step == no.small.steps) {
          message("Warning: Algorithm did not converge")
        }
        
        if (verbose) {
          print(paste0("converged", " at ", step, "th step"))
          print(paste0("PnEIC = ", colMeans(summ_eic[, -"id"])))
        }
        
        ## format output ----
        Psi_tmle <- Ht[, mget(c("A", grep("Psi.+", colnames(Ht), value = T)))]
        Psi_tmle <- Psi_tmle[, lapply(.SD, unique), by = "A"]
        Psi_tmle <- melt(Psi_tmle, id.vars = "A")
        Psi_tmle[, c("dummy", "event", "time") := tstrsplit(variable, "\\..")]
        Psi_tmle <- Psi_tmle[, -c("variable", "dummy")]
        Psi_tmle <- dcast(Psi_tmle, A + time ~ event, value.var = "value")
        Psi_tmle[, S := do.call(sum, mget(paste0(1:3))), by = c("A", "time")]
        Psi_tmle <- Psi_tmle[, lapply(.SD, as.numeric)][order(time, A)]
        
        tmle_ic <- summ_eic[, -c("id")]
        tmle_se <- sqrt(diag(var(tmle_ic)) / nrow(data))
        
        break
      }
    }
  }

# output --------------------------------------------------------------------------------------

# g-comp (sl estimate)
# unadjusted cox model
# tmle & ic
  
  return(list(estimates = list("tmle" = Psi_tmle, 
                               "g-comp" = Psi_init), 
              ic = tmle_ic, 
              se = tmle_se))
  
}




