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
## }

######################################################################
### prepare-pbc.R ends here
