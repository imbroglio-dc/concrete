##' bla
##'
##' blub
##' @title atitle
##' @param loss.fun should just be set to cox.loss.fun.
##' @param dt dataset.
##' @param V number of folds in cross-validation.
##' @param seed random seed :). 
##' @param delta.var name of event type variable. 
##' @param delta.value type of event of interest here-
##' @param treatment name of treatment variable. 
##' @param change.points specified if there is a changepoint in the effect of treatment across time.
##' @param cox.models a list of Cox models to be compared with cross-validation.
##' @return a
##' @seealso b
##' @examples c
##' @export 
##' @author Helene C. W. Rytgaard <hely@@biostat.ku.dk>
cox.sl <- function(loss.fun, dt, V=5, seed=19192,
                   delta.var=NULL, delta.value=NULL, treatment="A",
                   cox.models=list(mod1=list(Surv(time, delta==1)~A+L1+L2+L3),
                                   mod2=list(Surv(time, delta==1)~A+L1.squared+L2+L3),
                                   mod3=list(Surv(time, delta==1)~L2.squared+A+L1.squared+L2+L3),
                                   mod4=list(Surv(time, delta==1)~A+L1.squared),
                                   mod5=list(Surv(time, delta==1)~A*L1+L2+L3),
                                   mod6=list(Surv(time, delta==1)~A*L1.squared+L2+L3))) {
    
    cox.cve <- lapply(cox.models, function(cox.model) {
        return(cv.fun(loss.fun=cox.loss.fun, dt=dt, cox.model=cox.model[[1]],
                      delta.var=delta.var, delta.value=delta.value,
                      change.point=change.points, treatment=treatment))
    })
    
    picked.model <- unlist(cox.cve)[unlist(cox.cve)==min(unlist(cox.cve))]
    picked.cox.model <- list(form=cox.models[[gsub("\\.cve", "", names(picked.model))]][[1]],
                             cve=picked.model[[1]])
    
    return(list(picked.cox.model=picked.cox.model, 
                cox.cve.all=cox.cve))
    
}
