##' bla
##'
##' blub
##' @title atitle
##' @param loss.fun should be set to cox.loss.fun if comparing cox models; 
##' @param dt dataset. 
##' @param V number of folds in cross-validation. 
##' @param seed random seed :). 
##' @param X design matrix if one of c("coxnet", "cv.glmnet", "glmnet") is used.
##' @param Y outcome object (Surv)  if one of c("coxnet", "cv.glmnet", "glmnet") is used.
##' @param offset offset used to construct poisson hal. 
##' @param method.risk option to pick different cross-validation schemes for cox models; basically it picks the
##' risk set for the partial likelihood. Should be chosen as 'test'.
##' @param time.var name of time variable. 
##' @param penalty.factor variable to specify if certain (groups of) variables should not be penalized. 
##' @param delta.var name of event type variable. 
##' @param delta.value type of event of interest here-
##' @param change.point specified if there is a changepoint in the effect of treatment across time.
##' @param treatment name of treatment variable. 
##' @param cox.model model to compute the cve for. 
##' @param lambda.cv grid over which to choose penalization if one of c("coxnet", "cv.glmnet", "glmnet") is used.
##' @return a
##' @seealso b
##' @examples c
##' @export 
##' @author Helene C. W. Rytgaard <hely@@biostat.ku.dk>
cv.fun <- function(loss.fun, dt, V=5, seed=19192, X=NULL, Y=NULL, offset=NULL,
                   method.risk=c("test","train","VvH"), time.var=NULL,
                   penalty.factor=rep(1, ncol(X)), delta.var="delta", delta.value=1,
                   change.point=NULL, treatment=NULL, 
                   cox.model=NULL, lambda.cvs=c(sapply(1:5, function(jjj) (9:1)/(10^jjj)))) {
    
    set.seed(seed)
    
    n <- length(dt[, unique(id)])
    cv.split <- matrix(sample(1:n, size=n), ncol=V)
    
    cve.list <- lapply(1:V, function(vv) {
        
        test.set <- cv.split[,vv]
        train.set <- dt[, id][!dt[, id] %in% test.set]
        
        dt.train <- dt[id %in% train.set]
        
        if (length(X)==0) {
            train.fit <- coxph(cox.model, data=dt.train)
            return(loss.fun(train.fit=train.fit, dt=dt, risk.set=test.set, test.set=test.set, 
                            delta.var=delta.var, delta.value=delta.value, change.point=change.point))
        } else {
            
            if (length(offset)>0) {
                
                dt.train <- dt[id %in% train.set]
                
                dt.train[, D:=sum(time.obs==grid.time & get(delta.var)==delta.value), by="x"]
                dt.train[, risk.time:=grid.time-c(0, grid.time[-.N]), by="id"]
                dt.train[, RT:=sum(risk.time), by="x"]
                
                tmp.dt <- unique(dt.train[, c("x", "D", "RT")])
                Y.train <- tmp.dt[RT>0, D]
                offset.train <- tmp.dt[RT>0, log(RT)]
                X.train <- unique.matrix(X[dt$id %in% train.set,])[tmp.dt$RT>0,]
                
                if (length(lambda.cvs)>0) {
                    
                    train.fit <- glmnet(x=as.matrix(X.train), y=Y.train,
                                        offset=offset.train,
                                        family="poisson",
                                        maxit=1000,
                                        penalty.factor=penalty.factor,
                                        lambda=lambda.cvs)
                    
                } else {
                    
                    train.fit <- cv.glmnet(x=as.matrix(X.train), y=Y.train,
                                           offset=offset.train,
                                           family="poisson",
                                           penalty.factor=penalty.factor,
                                           maxit=1000)
                    
                    lambda.cvs <- train.fit$lambda.1se
                    
                }
                
                return(sapply(lambda.cvs, function(lambda.cv) {
                    loss.fun(train.fit=train.fit, dt=dt, test.set=test.set,
                             time.var=time.var, 
                             X=X, lambda.cv=lambda.cv, delta.var=delta.var, delta.value=delta.value)
                }))
                
            } else {
                
                Y.train <- Y[dt$id %in% train.set]
                X.train <- X[dt$id %in% train.set,]
                
                if (length(lambda.cvs)>0) {
                    
                    train.fit <- glmnet(x=as.matrix(X.train), y=Y.train,
                                        family="cox",
                                        maxit=1000,
                                        penalty.factor=penalty.factor,
                                        lambda=lambda.cvs)
                    
                } else {
                    
                    train.fit <- cv.glmnet(x=as.matrix(X.train), y=Y.train,
                                           family="cox",
                                           penalty.factor=penalty.factor,
                                           maxit=1000)
                    
                    lambda.cvs <- train.fit$lambda.1se
                    
                }
                
                return(sapply(lambda.cvs, function(lambda.cv)
                    loss.fun(train.fit=train.fit, dt=dt, risk.set=test.set, test.set=test.set,
                                 X=X, lambda.cv=lambda.cv, delta.var=delta.var, delta.value=delta.value)))
            }
            
            
        }       
    })
    
    if (length(X)==0 | length(lambda.cvs)==0) {
        return(list(cve=sum(unlist(cve.list))))
    } else {
        cve <- unlist(lapply(1:length(lambda.cvs), function(mm) {
            sum(unlist(lapply(cve.list, function(out) out[[mm]])))
        }))
        lambda.cv <- min(lambda.cvs[cve==min(cve)])
        return(list(min=list(lambda.cv=lambda.cv,
                             cve=unique(cve[cve==min(cve)])),
                    all=cbind(lambda=lambda.cvs, cve=cve)))
    } 
    
}
