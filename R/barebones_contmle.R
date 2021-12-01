

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
  
  
  # loop
  for (i in 1:length(cv_folds))
  {
    train_indices <- cv_folds[[i]]$training_set
    val_indices <- cv_folds[[i]]$validation_set
    train_data <- data[train_indices, -c("id")]
    val_data <- data[val_indices, -c("id")]
    
    for (models_j in models) {
      library_predict <- matrix(nrow = nrow(data), ncol = length(models_j))
      colnames(library_predict) <- names(models_j)
      for (model in models_j) {
        coxph_args <- list("formula" = model, "data" = train_data)
        model_fit <- do.call(coxph, coxph_args)
        
        model_predict <- predict(model_fit, newdata = val_data, type = "lp")
      }
    }
      
      ## train ----
      
      
      ## predict ----
      
    })
  }
  
  ## metalearner ----
  
  ## output initial estimates ----
  
  
  
  
  # Update step ---------------------------------------------------------------------------------
  
  # one-step tmle loop
  {
    ## clever covariate (A and C) ----
    
    ## clever covariate (Events) ----
    
    ## fluctuation
    
    ## update
  }
  
  
  
  
  # output --------------------------------------------------------------------------------------
  
  # g-comp (sl estimate)
  # unadjusted cox model
  # tmle & ic
  
  
}


