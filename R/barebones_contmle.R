
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
  
  
  ### cause-specific risks at target times ----
  ### target event (j), target time (k), cause (l) clever covariate ----
  for (k in target_times) {
    for (j in target_events) {
      #### event j risk at target time k ----
      Ht[, fjtau := get(paste0("F.j", j, ".t"))[TIME == k], by = c("id", "A")]
      for (l in target_events) {
        Ht[, hjlt := 1] # why is this the default?
        #### event component ----
        Ht[S.t > 0, hjlt := (l == j) - (fjtau - get(paste0("F.j", j, ".t"))) / S.t, 
           by = c("id", "A")]
        #### full clever covariate term ----
        Ht[S.t > 0, hjlt := hjlt * Ht.g * (TIME <= k) * (TIME <= T_tilde) * (ARM == A), 
           by = c("id", "A")]
        setnames(Ht, "hjlt", paste0("h.j", j, ".l", l, ".t", k))
      }#### mean treated and control event j risk at target time k ----
      Ht[, paste0("Psi.j", j, ".t", k) := mean(fjtau), by = "A"]
      setnames(Ht, "fjtau", paste0("F.j", j, ".t", k))
    }
  }
  
  ## initial estimator (g-computation) ----
  Psi_init <- c("A", grep("Psi", colnames(Ht), value = T))
  Psi_init <- Ht[, ..Psi_init][, lapply(.SD, unique), by = "A"]
  Psi_init <- melt(Psi_init, "A", variable.name = "dummy")[, c("Psi", "J", "t") := tstrsplit(dummy, "\\.\\w")][, -c("dummy")]
  Psi_init <- dcast(Psi_init, `t` ~ ..., value.var = "value")
  setnames(Psi_init, 2:ncol(Psi_init), 
           do.call(paste0, expand.grid("Psi.a", 0:1, ".j", target_events)))
  setcolorder(Psi_init, c("t", grep("a1", colnames(Psi_init), value = T)))
  Psi_init <- Psi_init[, `t` := as.numeric(`t`)][order(`t`)]
  
  
  # Update step ---------------------------------------------------------------------------------
  
  # one-step tmle loop
  {
    ## EIC ----
    Ht_eic <- copy(Ht[, .(id = unique(id))])
    for (k in target_times) {
      for (j in target_events) {
        Ht[, EIC := 0]
        for (l in target_events) {
          Ht[, EIC := EIC + get(paste0("h.j", j, ".l", l, ".t", k)) *  
               ((TIME == T_tilde ) * (EVENT == j) - 
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
    
    init_ic <- list("ic" = copy(Ht_eic), 
                    "se" = tail(apply(Ht_eic, 2, function(eic) {
                      sqrt(mean(eic^2)/nrow(Ht_eic))
                    }), -1))
    
    PnEIC <- mean(colMeans(Ht_eic[, -c("id")]))
    
    ## update steps (one-step) ----
    for (step in 1:no.small.steps) { 
      # Need to switch L and J? ---------
      ## making duplicate columns for cause-specific hazards and clever covs
      for (j in target_events) { # loop over target events
        mat[, (paste0("cox.j", j, ".tmp")) := get(paste0("cox.j", j))]
        for (k in target_times) { # loop over target times
          if (target.S) {
            mat[, (paste0("h", ".j", j, ".t", k, ".tmp")) :=
                  get(paste0("h", ".j", j, ".t", k))]
          }
          for (l in target_events) { # nested loop over events
            mat[, (paste0("h.j", j, ".l", l, ".t", k, ".tmp")) :=
                  get(paste0("h.j", j, ".l", l, ".t", k))]
          }
        }
      }
      
      if (target.S) {
        for (each in outcome.index) {
          fit.delta <- estimation[[each]][["event"]]
          mat[, (paste0("delta", fit.delta, ".dx")) := 0]
          for (kk in 1:length(tau)) {
            mat[, (paste0("delta", fit.delta, ".dx")) :=
                  get(paste0("delta", fit.delta, ".dx")) +
                  (get(time.var) <= tau[kk]) * Ht * ((get(
                    paste0("Ht", ".lambda", fit.delta, ".", kk)
                  )) * Pn.eic2[kk]) / Pn.eic.norm]
          }
        }
      } else { 
        ## 5.1 calculate update step direction? -------------------------------------
        
        for (each in outcome.index) { # loop over event hazards
          fit.delta <- estimation[[each]][["event"]]
          mat[, (paste0("delta", fit.delta, ".dx")) := 0] # fluctuation of cs-haz
          for (each2 in outcome.index[target]) { # loop over target events
            fit.delta2 <- estimation[[each2]][["event"]]
            for (kk in 1:length(tau)) { # loop over target times
              mat[, (paste0("delta", fit.delta, ".dx")) :=
                    get(paste0("delta", fit.delta, ".dx")) + 
                    (get(time.var) <= tau[kk]) * Ht * 
                    ((get( paste0( "Ht", fit.delta2, ".lambda", fit.delta, ".", kk ))) *
                       Pn.eic2[each2 == outcome.index[target]][[1]][kk]) / Pn.eic.norm]
            }
          }
        }
        
      }
      
      
      ## 5.2 set update epsilon down-scaling factor ----------------------------------------
      # why not just use deps.size?
      deps <- deps.size
      
      ## 5.3 update cause specific hazards ----------------------------------------
      mat[, surv.t := 1]
      for (each in outcome.index) { # loop over target events
        fit.delta <- estimation[[each]][["event"]]
        mat[, (paste0("fit.cox", fit.delta)) :=
              get(paste0("fit.cox", fit.delta)) *
              exp(deps * get(paste0("delta", fit.delta, ".dx")))]
        mat[get(paste0("fit.cox", fit.delta)) > 500, 
            (paste0("fit.cox", fit.delta)) := 500]
        # Does changing the order of event updates matter? -------------------------------
        mat[, surv.t := surv.t * exp(-cumsum(get(paste0("dhaz", fit.delta)) *
                                               get(paste0("fit.cox", fit.delta)))),
            by = c("id", A.name)]
        mat[get(paste0("fit.cox", fit.delta)) == Inf, surv.t := 0]
      }
      mat[, surv.t1 := c(0, surv.t[-.N]), by = c("id", A.name)]
      for (kk in 1:length(tau)) {
        mat[, (paste0("surv.tau", kk)) :=
              surv.t[get(time.var) == max(get(time.var)[get(time.var) <= tau[kk]])],
            by = c("id", A.name)]
      }
      for (each in outcome.index[target]) {
        fit.delta <- estimation[[each]][["event"]]
        mat[, (paste0("F", fit.delta, ".t")) := cumsum(surv.t * get(paste0("dhaz", fit.delta)) *
                                                         get(paste0("fit.cox", fit.delta))),
            by = c("id", A.name)]
        for (kk in 1:length(tau)) {
          mat[, (paste0("F", fit.delta, ".tau", kk)) :=
                get(paste0("F", fit.delta, ".t"))[get(time.var) == max(get(time.var)[
                  get(time.var) <= tau[kk]])], by = c("id", A.name)]
          for (each2 in outcome.index) {
            fit.delta2 <- estimation[[each2]][["event"]]
            if (fit.delta == fit.delta2) {
              mat[surv.t > 0, (paste0("Ht", fit.delta, ".lambda", fit.delta2, ".", kk)) :=
                    -(1 - (get(paste0("F", fit.delta, ".tau", kk)) - 
                             get(paste0("F", fit.delta, ".t"))) / surv.t)]
              mat[round(surv.t, 8) == 0, 
                  (paste0("Ht", fit.delta, ".lambda", fit.delta2, ".", kk)) := -1]
            } else {
              mat[surv.t > 0, (paste0("Ht", fit.delta, ".lambda", fit.delta2, ".", kk)) :=
                    (get(paste0("F", fit.delta, ".tau", kk)) - 
                       get(paste0("F", fit.delta, ".t"))) / surv.t]
              mat[round(surv.t, 8) == 0, 
                  (paste0("Ht", fit.delta, ".lambda", fit.delta2, ".", kk)) := 0]
            }
          }
        }
      }
      if (target.S) {
        for (each2 in outcome.index) {
          fit.delta2 <- estimation[[each2]][["event"]]
          for (kk in 1:length(tau)) {
            mat[, (paste0("Ht", ".lambda", fit.delta2, ".", kk)) := 0]
            for (each in outcome.index) {
              fit.delta <- estimation[[each]][["event"]]
              mat[, (paste0("Ht", ".lambda", fit.delta2, ".", kk)) := 
                    get(paste0("Ht", ".lambda", fit.delta2, ".", kk)) -
                    get(paste0("Ht", fit.delta, ".lambda", fit.delta2, ".", kk))]
            }
          }
        }
      }
      
      Pn.eic <- Pn.eic.fun(mat)
      Pn.eic2 <- Pn.eic2.fun(Pn.eic)
      Pn.eic.norm <- Pn.eic.norm.fun(Pn.eic2, Pn.eic)
      
      if (Pn.eic.norm.prev <= Pn.eic.norm) {
        if (cr) {
          for (each in outcome.index) {
            fit.delta <- estimation[[each]][["event"]]
            mat[, (paste0("fit.cox", fit.delta)) := get(paste0("fit.cox", fit.delta, ".tmp"))]
            for (kk in 1:length(tau)) {
              if (target.S)
                mat[, (paste0("Ht", ".lambda", fit.delta, ".", kk)) :=
                      get(paste0("Ht", ".lambda", fit.delta, ".", kk, ".tmp"))]
              for (each2 in outcome.index[target]) {
                fit.delta2 <- estimation[[each2]][["event"]]
                mat[, (paste0("Ht", fit.delta2, ".lambda", fit.delta, ".", kk)) :=
                      get(paste0("Ht", fit.delta2, ".lambda", fit.delta, ".", kk, ".tmp"))]
              }
            }
          }
        }
        
        Pn.eic <- Pn.eic.fun(mat)
        Pn.eic2 <- Pn.eic2.fun(Pn.eic)
        
        Pn.eic.norm <- Pn.eic.norm.fun(Pn.eic2, Pn.eic)
        deps.size <- 0.5 * deps.size # 0.1 * deps.size
        
      } else {
        Pn.eic.norm.prev <- Pn.eic.norm
        
        if (target.S) {
          Pn.eic3 <-
            sapply(1:length(Pn.eic), function(kk)
              Pn.eic[[kk]] / (ifelse(
                any(unlist(S.se.init) == 0), S.se.init[kk] + 0.001, S.se.init[kk]
              ) * sqrt(n)))
        } else {
          Pn.eic3 <-
            lapply(1:length(Pn.eic), function(kk)
              Pn.eic[[kk]] / (ifelse(any(unlist(init.ic) == 0), 
                                     init.ic[[kk]] + 0.001, 
                                     init.ic[[kk]]
              ) * sqrt(n)))
        }
      }
      
      if (target.S) {
        Pn.eic3 <-
          sapply(1:length(Pn.eic), function(kk)
            Pn.eic[[kk]] / (ifelse(
              any(unlist(S.se.init) == 0), S.se.init[kk] + 0.001, S.se.init[kk]
            ) * sqrt(n)))
      } else {
        Pn.eic3 <-
          lapply(1:length(Pn.eic), function(kk)
            Pn.eic[[kk]] / (ifelse(
              any(unlist(init.ic) == 0), init.ic[[kk]] + 0.001, init.ic[[kk]]
            ) * sqrt(n)))
      }
      
      if (cr &
          length(target) == length(outcome.index) & !target.S) {
        #--- here in fact want to check that we solve survival eic well enough!
        #--- i.e., we add it to our vector of Pn.eic3;
        Pn.eic3[[length(Pn.eic3) + 1]] <-
          sapply(1:length(Pn.eic3[[1]]), function(jj) {
            sum(sapply(Pn.eic, function(xx)
              xx[[jj]]))
          }) / (S.se.init * sqrt(n))
        check.sup.norm <- max(abs(unlist(Pn.eic3))) <= (criterion)
      } else {
        check.sup.norm <- max(abs(unlist(Pn.eic3))) <= (criterion)
      }
      
      if (no.update.if.solved) {
        # if (target.S) {
        #    Pn.eic3 <- sapply(1:length(Pn.eic), function(kk) Pn.eic[[kk]] / 
        #                        (ifelse(any(unlist(S.se.init) == 0), S.se.init[kk] + 0.001, 
        #                                S.se.init[kk])*sqrt(n)))
        # } else {
        #    Pn.eic3 <- lapply(1:length(Pn.eic), function(kk) Pn.eic[[kk]] / 
        #                        (ifelse(any(unlist(init.ic) == 0), init.ic[[kk]] + 0.001, 
        #                                init.ic[[kk]])*sqrt(n)))
        # }
        
        if (length(target) > 1) {
          Pn.eic2 <- lapply(1:length(Pn.eic2), function(kk) {
            tmp <- Pn.eic2[[kk]]
            tmp[abs(Pn.eic3[[kk]]) <= criterion] <- 0
            return(tmp)
          })
        } else {
          Pn.eic2[abs(Pn.eic3) <= criterion] <- 0
        }
        #Pn.eic.norm <- Pn.eic.norm.fun(Pn.eic2, Pn.eic)
        
        if (all(unlist(Pn.eic2) == 0))
          check.sup.norm <-
            TRUE
        else
          Pn.eic.norm <- Pn.eic.norm.fun(Pn.eic2, Pn.eic)
        print(Pn.eic2)
        
      }
      
      if (check.sup.norm |
          step == no.small.steps) {
        # (left.criterion <= criterion) {
        if (step == no.small.steps) {
          message("Warning: Algorithm did not converge")
        }
        
        if (verbose)
          print(paste0("converged", " at ", step, "th step"))
        if (verbose)
          print(paste0("eic = ", Pn.eic.fun(mat)))
        
        #-- 12c -- compute sd:
        if (cr) {
          if (treat.effect[1] == "stochastic") {
            final.fit <- lapply(target, function(each) {
              sapply(1:length(tau), function(kk) {
                mean(rowSums(sapply(a, function(aa)
                  (
                    mat[get(A.name) == aa, pi.star[1] *
                          get(paste0("F", estimation[[outcome.index[each]]][["event"]],
                                     ".tau", kk))[1],
                        by = "id"][, 2][[1]]
                  ))))
              })
            })
          } else {
            final.fit <- lapply(outcome.index[target], function(each) {
              sapply(1:length(tau), function(kk) {
                mean(rowSums(sapply(a, function(aa)
                  (2 * (aa == a[1]) - 1) * (mat[get(A.name) == aa,
                                                get(paste0("F", estimation[[each]][["event"]],
                                                           ".tau", kk))[1],
                                                by = "id"][, 2][[1]]))))
              })
            })
          }
          names(final.fit) <-
            paste0("F", sapply(outcome.index[target], function(each)
              estimation[[each]][["event"]]))
          final.ic <-
            eval.ic(mat, final.fit, target.index = outcome.index[target])
        }
        
        final.list <-
          lapply(1:length(final.fit), function(each.index) {
            out <- rbind(tmle.est = final.fit[[each.index]],
                         tmle.se = final.ic[[each.index]])
            colnames(out) <- paste0("tau=", tau)
            return(out)
          })
        
        names(final.list) <-
          paste0("F", sapply(outcome.index[target], function(each)
            estimation[[each]]["event"]))
        
        if (length(target) == length(outcome.index) &
            length(outcome.index) > 1) {
          S.se <-
            sapply(tau, function(tt)
              sqrt(mean(rowSums(
                sapply(1:length(outcome.index), function(target11) {
                  unlist(
                    eval.ic( mat,
                             unlist(lapply(final.list, function(xx)
                               xx["tmle.est", paste0("tau=", tt)]
                             )),
                             target.index = outcome.index[target11],
                             tau.values = tt,
                             survival = TRUE
                    )
                  )
                })
              )^2) / n))
          S.fit <- sapply(tau, function(tt)
            (treat.effect != "ate") - sum(sapply(final.list, function(fl)
              fl["tmle.est", paste0("tau=", tt)])))
          final.S <- rbind(S.fit, S.se)
          colnames(final.S) <- paste0("tau=", tau)
          rownames(final.S) <- c("tmle.est", "tmle.se")
          final.list$S <- final.S
          if (target.S) {
            tmle.list$tmle <- final.list[names(final.list) == "S"]
          } else {
            tmle.list$tmle <- final.list
          }
        } else {
          tmle.list$tmle <- final.list
        }
        
        if (!second.round) {
          if (check.sup)
            tmle.list$check.sup.norm <- list(
              check.sup.norm = check.sup.norm,
              lhs = max(abs(unlist(Pn.eic3))),
              rhs = criterion / sqrt(length(target) * length(tau))
            )
          tmle.list$convergenced.at.step <- step
          
          if (!push.criterion)
            break
          else
            criterion <- criterion / sqrt(n)
          second.round <- TRUE
        } else {
          if (target.S) {
            tmle.list$tmle.second.round <- final.list[names(final.list) == "S"]
          } else {
            tmle.list$tmle.second.round <- final.list
          }
          if (check.sup)
            tmle.list$check.sup.norm.second.round <- list(
              check.sup.norm =
                max(abs(unlist(Pn.eic3))) <= (criterion),
              lhs = max(abs(unlist(Pn.eic3))),
              rhs = criterion
            )
          tmle.list$convergenced.at.step.second.round <- step
          break
        }
      }
    }
    
  }
  
  
  
  
  # output --------------------------------------------------------------------------------------
  
  # g-comp (sl estimate)
  # unadjusted cox model
  # tmle & ic
  
  
}


