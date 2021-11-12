library(here); library(tidyverse); library(mice)
# read in subject level covariate dataset
W <- haven::read_sas("data/ADaM/adsl.sas7bdat")

# 9341th subject in subject level dataset, not in the time-to-event dataset
#   FASFL = N, so supposed to be taken out?
W <- subset(W, subset = USUBJID != unique(W$USUBJID)[ # remove the subject who is not in the time-to-event dataset
    which(!(unique(W$USUBJID) %in% unique(haven::read_sas("data/ADaM/adtte.sas7bdat")$USUBJID)))])

# REMOVE UNITS, ID VARIABLES, AND POST-BASELINE VARIABLES
W <- dplyr::select(W, -c(STUDYID, SUBJID, DIABDURU, WSTCRBLU,
                         BMIBLU, WGTBLU, WGTTBLU, PULSEBLU, SYSBPBLU,
                         DIABPBLU, HBA1CBLU, HDL1BLU, LDL1BLU,
                         CHOL1BLU, TRIG1BLU, HBC2BLU, HDL2BLU,
                         LDL2BLU, CHOL2BLU, TRIG2BLU, CREATBLU,
                         EGFREPIU, EGFREPBU, HGTTBLU, AGEU, HGTBLU)) %>%
    dplyr::select(-c(ISHRGRL, ISHRGRN, DTHDT, EOSDT, EOTDT, LSTCONDT, LSTSVDT,
                     RANDDT, BRTHDT, EOSSTT, TOTREAT, STDURY, TRDURD, TRDURY,
                     TRDU15OD, TRDUROY, TRDU15OY))

# REMOVE DEGENERATE COLUMNS
W <- select_if(W, ~length(unique(.))>1)

# REMOVE DUPLICATE COLUMNS & LINEARLY DEPENDENT COLUMNS
W <- dplyr::select(W, -c(PREVTFL, WGTBL, HBC2BL, ARMCD,
                         HDL2BL, LDL2BL, CHOL2BL, TRIG2BL,
                         GERDCMFL, # near dup GERDBLFL, off by 2
                         STROKSFL, # remove the sensitivity analysis stroke flag?,
                         CVHIFL, CVMEDFL, CVRISKN, age_category))

# Collinearity? - 4 renal failure columns (RENFSEV, RENFCKD, RNFSCSEV, RNFSCCKD)
# dplyr::select(W, RENFSEV, RENFCKD, RNFSCSEV, RNFSCCKD) %>% group_by_all() %>% count() %>% view()

# excessive missingness (ISMHIBFL [99.3%], ISMHGPFL [99.9%], RETINSEV [77.9%], NYHACLAS[82.3%])
# skimr::skim(dplyr::select(W, ISMHIBFL, ISMHGPFL, RETINSEV, NYHACLAS))

# remove vars with >99% missingness
W <- dplyr::select(W, -c(ISMHIBFL, ISMHGPFL))

# how to deal with missing retinopathy severity & NYHA class?
W <- W %>%
    mutate(RETINSEV = case_when(RETINSEV == "" ~ "NA",
                                T ~ RETINSEV),
           RETINSEV = factor(RETINSEV, ordered = T,
                             levels = c("NA", "no retinopathy", "non-proliferative", "proliferative")),
           NYHACLAS = case_when(NYHACLAS == "" ~ "NA",
                                T ~ NYHACLAS),
           NYHACLAS = factor(NYHACLAS, ordered = T,
                             levels = c("NA", "NYHA CLASS I", "NYHA CLASS II", "NYHA CLASS III")))

# make logical variables logicals
W <- W %>%
    mutate_if(~mean(unique(.) %in% c("Y", "N")) == 1, ~case_when(. == "Y" ~ T,
                                                                 T ~ F)) %>%
    mutate(SEX = (SEX == "M"),
           ARM = (ARM == "Liraglutide"),
           ETHNIC = ETHNIC == "HISPANIC OR LATINO") %>%
    rename(MALE = SEX,
           HISPANIC = ETHNIC)

# make factor variables factors, make SMOKER and renal disease cols ordered
W <- W %>% mutate(SMOKER = factor(SMOKER, ordered = T,
                                  levels = c("NEVER SMOKED", "PREVIOUS SMOKER",
                                             "CURRENT SMOKER"))) %>%
    mutate_at(vars(starts_with("RNF"), starts_with("RENF")),
              ~factor(., ordered = T,
                      levels = c("Normal (EGFR>=90)", "Mild (EGFR<90)",
                                 "Moderate (EGFR<60)", "Severe (EGFR<30)"))) %>%
    mutate_if(is.character, as_factor)


# imputation --------------------------------------------------------------

# imputed <- mice(W, m = 20, maxit = 1, visitSequence = "monotone")

# skimr::skim(W)
# view(W[1:10, ])
