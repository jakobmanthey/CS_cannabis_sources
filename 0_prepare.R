# ==================================================================================================================================================================
# ==================================================================================================================================================================
# ==================================================================================================================================================================

# PROJECT TITLE:  CANNASTREET
# CODE AUTHOR:    JAKOB
# DATE STARTED:   240227

# ==================================================================================================================================================================
# ==================================================================================================================================================================
# ==================================================================================================================================================================

# 0) ESSENTIALS
# ______________________________________________________________________________________________________________________

# clean workspace
rm(list=ls())

packages <- c("data.table", "ggplot2", "ggthemes",
              "Hmisc", "tidyr", "nnet",
              "tableone", "knitr", "kableExtra")

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}
# Load packages
invisible(lapply(packages, library, character.only = TRUE))

# current date:
DATE <- format(Sys.Date(), "%Y%m%d")

# themes and options
theme_set( theme_gdocs() )
options(scipen = 999)

# ==================================================================================================================================================================
# ==================================================================================================================================================================
# ==================================================================================================================================================================

# 1) LOAD DATA (version 20 Mar 2025)
# ______________________________________________________________________________________________________________________

##  Germany WAVE 2
# -------------------------------------------------------

path <- file.path("data", "20250321","full_GER_wave2_cleaned_20250321.rds")
input <- data.table(readRDS(file = path))

nrow(input) # 11578

# GSZB2 = Grundstichprobe Zwischenbericht 2
# ZS = Zielstichprobe
# AS1 = Auswertungsstichprobe 1
# AS2 = Auswertungsstichprobe 2 (DUIC)
# ASkom = Auswertungsstichprobe kombiniert

# ==================================================================================================================================================================
# ==================================================================================================================================================================
# ==================================================================================================================================================================

## 2) PREPARE DATA
# ______________________________________________________________________________________________________________________

## 2.1) get data
#-------------------------------------------------------

data <- copy(input[,.(ID, weights,
                      sex, 
                      age = age_in_years,
                      agegroup = agegroup_basicsample,
                      edu = edu_group,
                      gisd_k,
                      DEGURBA = DGURBA,
                      can_freq,
                      ndays_flowers,ndays_resin,ndays_edibles,ndays_concentrate_liquids,
                      QUANT.1,QUANT.2,QUANT.3,QUANT.4,
                      spend_monthly,
                      age_firstuse,
                      MEDICALUSE.01,MEDICALUSE.02,
                      CAST_score_full,
                      alcohol_freq,AUDITC2,AUDITC3,
                      HEALTH.01,HEALTH.02,
                      DISTRESS_total,
                      SOURCE.1,SOURCE.2,SOURCE.3,SOURCE.4,SOURCE.5_1,SOURCE.5_2,
                      SOURCE.6,SOURCE.7,SOURCE.8,SOURCE.9,SOURCE.11.open)])

# convert haven to standard factors
data[, (names(data)) := lapply(.SD, function(x) if(haven::is.labelled(x)) haven::as_factor(x) else x), .SDcols = names(data)]


## 2.2) define OUTCOMES
#-------------------------------------------------------

# source_pharma
data$SOURCE.1
data$source_pharma <- data$SOURCE.1 == "ausgewählt"
data[, table(SOURCE.1,source_pharma, useNA = "always")]
#data$SOURCE.1 <- NULL

# source_social
data$SOURCE.2
data$source_social <- data$SOURCE.2 == "ausgewählt"
data[, table(SOURCE.2, source_social, useNA = "always")]
#data$SOURCE.2 <- NULL

# source_illegal
data$SOURCE.3
data$SOURCE.4
data$SOURCE.7
data$SOURCE.8
data$source_illegal <- haven::as_factor(data$SOURCE.3) == "ausgewählt" | data$SOURCE.4 == "ausgewählt" | data$SOURCE.7 == "ausgewählt" | data$SOURCE.8 == "ausgewählt"
data[, table(SOURCE.3, source_illegal, useNA = "always")]
data[, table(SOURCE.4, source_illegal, useNA = "always")]
data[, table(SOURCE.7, source_illegal, useNA = "always")]
data[, table(SOURCE.8, source_illegal, useNA = "always")]
#data$SOURCE.3 <- data$SOURCE.4 <- data$SOURCE.7 <- data$SOURCE.8 <- NULL

# source_homegrow (including AV)
data$SOURCE.5_1
data$SOURCE.5_2
data$SOURCE.6
data$source_homegrow <- data$SOURCE.5_1 == "ausgewählt" | data$SOURCE.5_2 == "ausgewählt" | data$SOURCE.6 == "ausgewählt"
data[, table(SOURCE.5_1, source_homegrow, useNA = "always")]
data[, table(SOURCE.5_2, source_homegrow, useNA = "always")]
data[, table(SOURCE.6, source_homegrow, useNA = "always")]
#data$SOURCE.5_1 <- data$SOURCE.5_2 <- data$SOURCE.6 <- NULL

# source_other
data$SOURCE.9
data$SOURCE.11.open
data$source_other <- data$SOURCE.9 == "ausgewählt" | data$SOURCE.11.open %like% "[a-zA-Z]"
data[, table(SOURCE.11.open, source_other, useNA = "always")]
data[, table(SOURCE.9, source_other, useNA = "always")]
data[ID %in% c(223997,211379)] # only other but no information
data[ID %in% c(223997,211379), source_other := T]

#data$SOURCE.9 <- NULL
#data$SOURCE.11.open <- NULL


##  all as numerics
data$source_illegal_num <- as.numeric(data$source_illegal)+1
data$source_social_num <- as.numeric(data$source_social)+1
data$source_homegrow_num <- as.numeric(data$source_homegrow)+1
data$source_pharma_num <- as.numeric(data$source_pharma)+1
data$source_other_num <- as.numeric(data$source_other)+1

data$s_pharmacy <- as.numeric(data$SOURCE.1)
data$s_friends <- as.numeric(data$SOURCE.2)
data$s_unknowns <- as.numeric(data$SOURCE.3)
data$s_knowndealer <- as.numeric(data$SOURCE.4)
data$s_homegrown_self <- as.numeric(data$SOURCE.5_1)
data$s_homegrown_others <- as.numeric(data$SOURCE.5_2)
data$s_cannabisclub <- as.numeric(data$SOURCE.6)
data$s_online <- as.numeric(data$SOURCE.7)
data$s_socialmedia <- as.numeric(data$SOURCE.8)
data$s_no_possess <- as.numeric(data$SOURCE.9)
data$s_other <- as.numeric(data$SOURCE.11.open %like% "[a-zA-Z]")+1

##  ONLY homegrow and no other source
data$temp <- rowSums(data[,.SD, .SDcols = names(data)[names(data) %like% "^s_"]])
data[, only_homegrow := ifelse(s_homegrown_self == 2 & temp == 12, T, F)]
data[, only_medical := ifelse(s_pharmacy == 2 & temp == 12, T, F)]
data[, only_homegrow_or_medical := ifelse(s_homegrown_self == 2 & s_pharmacy == 2 & temp == 13, T, F)]
data[, table(only_homegrow,only_medical)]
data[, table(only_homegrow,only_homegrow_or_medical)]


## 2.3) define COVARIATES - SOCIODEMOGRAPHICS
#-------------------------------------------------------

##  SEX
data$sex
levels(data$sex) <- c("women","men","other")

##  AGEGROUPS
data$agegroup
data[, table(age,agegroup, useNA = "always")]

##  EDUCATION
data$edu

##  DEGURBA
data$DEGURBA
data[, prop.table(table(can_freq, DEGURBA),2)]

##  GISD
input[, table(gisd_5,gisd_k)]
input[, mean(gisd_score), by = gisd_k] # 1 = low deprivation, 3 = high deprivation
data$gisd_k <- factor(data$gisd_k, levels = c(1,2,3), labels = c("low","mid","high"))
data[, prop.table(table(can_freq, gisd_k),2)]



## 2.4) define COVARIATES - CANNABIS
#-------------------------------------------------------

##  FREQ12M
data$can_freq
data[can_freq != "gar nicht", freq12m := factor(can_freq)]
levels(data$freq12m) <- c("less than monthly","1+ monthly", "1+ weekly", "(near) daily")
data[, can12m := can_freq != "gar nicht"]
data[, table(can12m)] # 1490

##  FREQ30D
table(data$ndays_flowers)
table(data$ndays_resin)
table(data$ndays_edibles)
table(data$ndays_concentrate_liquids)
select <- names(data)[names(data) %like% "^ndays"]

data[, ndays_any := apply(.SD, 1, function(x) any(!is.na(x))), .SDcols = select]
data[, (select) := lapply(.SD, function(x) as.numeric(gsub("[a-zA-Z]", "", x))), .SDcols = select]

data[ndays_any == T, freq30d := apply(.SD, 1, max, na.rm = T), .SDcols = select]
data[, freq30d := ifelse(freq30d == "-Inf",0,freq30d)]
data[, table(freq30d)]
data[, table(ndays_flowers,freq30d)]
data[, can30d := ifelse(is.na(freq30d), F, freq30d > 0)]
data$ndays_flowers <- data$ndays_resin <- data$ndays_edibles <- data$ndays_concentrate_liquids <- data$ndays_any <- NULL

##  FREQ30D GROUPING (freq30d_group)
cuts <- data[!is.na(freq30d), quantile(freq30d,c(0,0.2,0.4,0.6,0.8,1))]
data[, freq30d_group := factor(cut(freq30d,cuts, include.lowest = T))]
data[, table(freq30d,freq30d_group, useNA = "always")]
levels(data$freq30d_group) <- gsub(","," to <=",levels(data$freq30d_group))
levels(data$freq30d_group) <- gsub("\\(",">",levels(data$freq30d_group))
levels(data$freq30d_group) <- gsub("\\]|\\[","",levels(data$freq30d_group))

##  QUANTITY (quant_tot)
data$QUANT.1 # total 30 days
data$QUANT.2 # total weekly
data$QUANT.3 # total daily
data$QUANT.4 # do not know
data[, quant_tot := ifelse(!is.na(QUANT.1), QUANT.1,
                           ifelse(!is.na(QUANT.2), QUANT.2 * 4.25,
                                  ifelse(!is.na(QUANT.3), QUANT.3 * 30, NA)))]
data[can30d == T, table(quant_tot)]
data[can30d == T, quant_tot := ifelse(quant_tot == 0, NA, quant_tot)]

##  QUANTITY PER USE DAY CORRECTION --> not more than 150g per 30 days (5g per day)
data[, table(quant_tot > 150)] # 14
data[, quant_tot := ifelse(quant_tot > 150, 150, quant_tot)]

##  QUANTITY GROUPING (quant_group)
data[!is.na(quant_tot), summary(quant_tot)]
cuts <- data[quant_tot>0, quantile(quant_tot,c(0,0.25,0.5,0.75,1))]
data[, quant_group := factor(cut(quant_tot,cuts,include.lowest = T,dig.lab = 4))]
data[, table(quant_tot,quant_group, useNA = "always")]
levels(data$quant_group) <- gsub(","," to <=",levels(data$quant_group))
levels(data$quant_group) <- gsub("\\(",">",levels(data$quant_group))
levels(data$quant_group) <- gsub("\\]|\\[","",levels(data$quant_group))

##  QUANTITY PER USE DAY (quant_pd)
data[can30d == T, quant_pd := quant_tot/freq30d]
data[, table(quant_pd,can30d)]
#data[quant_pd == 0, quant_pd := NA]
data$QUANT.1 <- data$QUANT.2 <- data$QUANT.3 <- data$QUANT.4 <- NULL

##  QUANTITY PER USE DAY GROUPING (quant_pd_group)
cuts <- data[can30d == T & !is.na(quant_pd), quantile(quant_pd,c(0,0.25,0.5,0.75,1))]
data[, quant_pd_group := factor(cut(quant_pd,cuts, include.lowest = T))]
data[, table(quant_pd,quant_pd_group, useNA = "always")]
levels(data$quant_pd_group) <- gsub(","," to <=",levels(data$quant_pd_group))
levels(data$quant_pd_group) <- gsub("\\(",">",levels(data$quant_pd_group))
levels(data$quant_pd_group) <- gsub("\\]|\\[","",levels(data$quant_pd_group))

##  SPEND
data$spend_monthly
summary(data$spend_monthly)

##  AGE ONSET
data$age_firstuse
data[, earlyonset := age_firstuse <=16]
data[, table(age_firstuse,earlyonset)]

##  MEDICAL
data$MEDICALUSE.01
data$purpose <- factor(data$MEDICALUSE.01)
levels(data$purpose)  <- c("only medical","medical and non-medical", "only non-medical")
data$purpose <- factor(data$purpose, levels = c("only non-medical","medical and non-medical","only medical"))
data$MEDICALUSE.01 <- NULL
table(data$purpose)

##  PRESCRIPTION
data$MEDICALUSE.02
#data[, prescribed := factor(MEDICALUSE.02)]
data[!is.na(MEDICALUSE.02), prescribed := ifelse(MEDICALUSE.02 %like% "Ja, mir wurde medizinisches", T, F)]
data[, table(prescribed)]

##  CAST
data[, table(CAST_score_full)]
data[, castrisk := CAST_score_full >= 7]

##  HEALTH
data[, table(HEALTH.01, HEALTH.02)]
data[, prop.table(table(HEALTH.01, HEALTH.02),1)]
data[, health_good := HEALTH.01 <= 2]
data[, health_chronic := HEALTH.02 == 1]
data[, prop.table(table(health_good, health_chronic),1)]

##  DISTRESS
data[, table(DISTRESS_total)]
data[!is.na(DISTRESS_total), distress := DISTRESS_total >= 13]
data$DISTRESS_total <- NULL

## 2.5) define SAMPLE FOR REGRESSION
#-------------------------------------------------------

##  sample for regressions:
data[,regsample := can30d == T & complete.cases(quant_pd_group)]
data[, table(regsample)] # 723

##  AUDIT-C???


# ==================================================================================================================================================================
# ==================================================================================================================================================================
# ==================================================================================================================================================================

## 3) OUTPUT
# ______________________________________________________________________________________________________________________

saveRDS(data, paste0("data/data_sources_prepared_",DATE,".rds"))

