### prepare-pbc.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Apr 20 2022 (14:58) 
## Version: 
## Last-Updated: Apr 20 2022 (14:59) 
##           By: Thomas Alexander Gerds
##     Update #: 1
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
## preparePBC <- function(){
data <- as.data.table(survival::pbc)
set.seed(12345)
data[is.na(trt), trt := sample(1:2, sum(is.na(trt)), replace = TRUE)][, trt := trt - 1]
data[, status := as.numeric(status >= 1)]
intervention <- list("A == 1" = list("intervention" = function(a, L) {rep_len(1, length(a))},
                                     "g.star" = function(a, L) {as.numeric(a == 1)}),
                     "A == 0" = list("intervention" = function(a, L) {rep_len(0, length(a))},
                                     "g.star" = function(a, L) {as.numeric(a == 0)}))
target.time <- 500 * (2:4)
target.event <- sort(unique(data[status > 0, status]))
logreg <- make_learner(Lrnr_glm)
lasso <- make_learner(Lrnr_glmnet) # alpha default is 1
ridge <- Lrnr_glmnet$new(alpha = 0)
e_net <- make_learner(Lrnr_glmnet, alpha = 0.5)
a_lrnrs <- logreg # make_learner(Stack, logreg, lasso, ridge, e_net)
model <- list("Trt" = a_lrnrs,
              "0" = list(mod1 = Surv(time, status == 0) ~ trt + age + sex),
              "1" = list(mod1 = Surv(time, status == 1) ~ trt + age + sex))
estimation <- list("cause1" = list(fit = "cox",
                                   model = Surv(time, status == 1) ~ trt + age + sex),
                   "cens" = list(fit = "cox",
                                 model = Surv(time, status == 0) ~ trt + age + sex))
## }

######################################################################
### prepare-pbc.R ends here
