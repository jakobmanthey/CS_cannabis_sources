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
              "tableone", "knitr", "kableExtra",
              "survey", "gtsummary", "survey", "dplyr")

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}
# Load packages
invisible(lapply(packages, library, character.only = TRUE))

# current date:
DATE <- format(Sys.Date(), "%y%m%d")

# themes and options
theme_set( theme_gdocs() )
options(scipen = 999)


# ==================================================================================================================================================================
# ==================================================================================================================================================================
# ==================================================================================================================================================================

# 1) PREPARE DATA
# ______________________________________________________________________________________________________________________

##  get data (after running LCA)
# -------------------------------------------------------

#source("1_LCA.R")

data <- readRDS(paste0("data/data_sources_prepared_after LCA_",DATE,".rds"))
lcamod <- readRDS(paste0("data/LCA model_",DATE,".rds"))

##  helper functions
# -------------------------------------------------------

# Custom test functions for survey data
custom_ttest <- function(data, variable, by, ...) {
  svyttest(as.formula(paste(variable, "~", by)), design = data) %>% 
    broom::tidy()
}

custom_chisq <- function(data, variable, by, ...) {
  svychisq(as.formula(paste("~", variable, "+", by)), design = data) %>% 
    broom::tidy()
}



# ==================================================================================================================================================================
# ==================================================================================================================================================================
# ==================================================================================================================================================================

## 2) RESULTS
# ______________________________________________________________________________________________________________________

## 2.0) Methods
#-------------------------------------------------------

#input[, min(START)]
#input[, max(END)]
data[,summary(weights)]

## 2.1) Descriptives general
#-------------------------------------------------------

data[, table(can12m)] # 1478
data[, table(can30d)] # 1073

data[, wtd.mean(can12m, weights)]
data[, wtd.mean(can30d, weights)]

svy <- svydesign(ids = ~ID, weights = ~weights, data = data)
svyciprop(~can12m, svy)
svyciprop(~can30d, svy)


## 2.2) A priori types of cannabis sourcing
#-------------------------------------------------------

# more than 1 source type:
data$source_atleast2 <- NA
data[can12m == T]$source_atleast2 <- rowSums(data[can12m == T, .(source_illegal, source_social, 
                                                                 source_homegrow, source_pharma, source_other)]) >1
# (378+92+4)/1478 = 32.1%
prop.table(xtabs(weights ~ source_atleast2, data = data)) # 35.0%
svy <- svydesign(ids = ~ID, weights = ~weights, data = data)
svyciprop(~source_atleast2, svy) # ci

# home cultivation more often in subsample?
data[can12m == T, weights::wtd.chi.sq(var1 = source_homegrow,
                                      var2 = regsample, weight = weights)]



# strictly legal
data[can12m == T, .(sum(source_homegrow == T | source_pharma == T),
                    mean(source_homegrow == T | source_pharma == T),
                    wtd.mean(source_homegrow == T | source_pharma == T, weights))] # 50.7%
data[can12m == T, .(sum(only_homegrow), 
                    mean(only_homegrow == T),
                    wtd.mean(only_homegrow == T, weights))]
data[can12m == T, .(sum(only_cannabisclub),
                    mean(only_cannabisclub == T),
                    wtd.mean(only_cannabisclub == T, weights))]
data[can12m == T, .(sum(only_medical), 
                    mean(only_medical == T),
                    wtd.mean(only_medical == T, weights))]

##  Cannabis sourcing by medical/recreational use purpose
data[, .(.N, mean = wtd.mean(source_pharma, weights)), by = purpose][order(purpose)]
data[, weights::wtd.chi.sq(var1 = source_pharma,
                           var2 = purpose, weight = weights)]



## 2.3) Multinomial regression
#-------------------------------------------------------

vars <- names(data)[names(data) %like% "sex|agegr|edu|DEGURBA|gisd|distress|health|freq12m|purpose|prescribed|castrisk"]
f <- as.formula(paste0("class ~ ", paste0(vars, collapse = " + ")))

sub <- copy(data[can12m == T & sex != "other"])
sub$sex <-factor(sub$sex)

mnmod <- multinom(formula = f,data = sub)
summary(mnmod)

set.seed(924)

# Risk Ratios
rs <- as.data.frame(exp(coef(mnmod))) 
rs <- data.table(class = rownames(rs), rs)
rs <- melt(rs, id.vars = c("class"), value.name = "riskratio")[order(class,variable)]

# confidence intervals:
ci <- as.data.table(round(exp(confint(mnmod)),3))
names(ci) <- c("variable","bound","class","value")
ci <- dcast(ci, class + variable ~ bound)
ci$ci <- paste0(ci$'2.5 %'," to ", ci$'97.5 %')

# p values
sm <- summary(mnmod)
z <- summary(mnmod)$coefficients/summary(mnmod)$standard.errors
p <- (1 - pnorm(abs(z), 0, 1)) * 2
p <- data.table(class = rownames(p), p)
p <- melt(p, id.vars = c("class"), value.name = "p")[order(class,variable)]

regout <- merge(rs, ci, by = c("class","variable"))
regout <- merge(regout, p, by = c("class","variable"))

##  var labels for Table and Figure
regout[, varlab := dplyr::recode_factor(variable,
                                        "(Intercept)" = "Intercept",
                                        "sexmen" = "Sex: men",
                                        "agegroup25-34" = "Age: 25-34",
                                        "agegroup35-44" = "Age: 35-44",
                                        "agegroup45-54" = "Age: 45-54",
                                        "agegroup55-64" = "Age: 55-64",
                                        "edumid" = "Education: mid",
                                        "eduhigh" = "Education: high",
                                        "gisd_kmid" = "Regional deprivation: mid",
                                        "gisd_khigh" = "Regional deprivation: high",
                                        "DEGURBArural" = "DEGURBA: rural",
                                        "DEGURBAtowns_suburbs" = "DEGURBA: towns and suburbs",
                                        "freq12m1+ monthly"  = "Use frequency: at least monthly",
                                        "freq12m1+ weekly"  = "Use frequency: at least weekly",
                                        "freq12m(near) daily"  = "Use frequency: (near) daily",
                                        "purposemedical and non-medical" = "Use purpose: medical and non-medical",
                                        "purposeonly medical" = "Use purpose: only medical",
                                        "prescribedTRUE" = "Has prescription for medical cannabis",
                                        "castriskTRUE" = "CAST: high CUD risk",
                                        "distressTRUE" = "K6: high distress",
                                        "health_goodTRUE" = "Self-reported good health",
                                        "health_chronicTRUE" = "Self-reported chronic disease")]
regout[, table(varlab, useNA = "always")]

kable(regout[order(class,varlab),.(class,varlab,riskratio,ci,p)]) %>%
  kable_styling() %>%
  save_kable(file = paste0("tabs/sup_tab5_",DATE,".html"))


## 2.4) Cannabis sourcing and use quantities
#-------------------------------------------------------

# total quantities
data[regsample == T, sum(quant_tot)] # 11724
data[regsample == T, sum(quant_tot*weights)] # 14033

# distribution of classes in subsample
data[regsample == T, .(weighted_sum = sum(weights)), by = class][
  , .(class,weighted_prop = weighted_sum / sum(weighted_sum))][order(class)]

# mean/median quantity
data[regsample == T, .(mean = wtd.mean(quant_tot, weights), median = median(quant_tot)), by = class][order(mean)]

ggplot(data[regsample == T], aes(x = quant_tot)) + 
  geom_histogram()

mean(data[regsample == T]$quant_tot)
var(data[regsample == T]$quant_tot)

# negative binomial
testquant <- glm.nb(quant_tot ~ class, data[regsample == T])
summary(testquant)

# total quantities by class
temp  <- data[regsample == T,.(ID, class, quant_tot, weights)]
temp[, wt_tot := sum(quant_tot*weights)]
temp[, wt_class := sum(quant_tot*weights), by = class]
temp <- unique(temp[, .(class,prop = wt_class/wt_tot)])[order(prop)]

#social+pharmacy+homegrow:
temp[class %like% "SOCIAL|PHARMACY|HOME CULTIVATION", sum(prop)]

# mean/median €
data[regsample == T, .(mean = wtd.mean(spend_monthly, weights), median = median(spend_monthly)), by = class][order(mean)]

ggplot(data[regsample == T], aes(x = spend_monthly)) + 
  geom_histogram()









# ==================================================================================================================================================================
# ==================================================================================================================================================================
# ==================================================================================================================================================================

## 3) TABLES
# ______________________________________________________________________________________________________________________

## 3.1) TABLE 1
#-------------------------------------------------------

# Define variables and create survey design
vars <- names(data)[names(data) %like% "sex|age|edu|gisd|DEGURBA|freq12m|freq30d|quant_pd|earlyonset|purpose|prescribed|castrisk|distress|health"]
subtab <- data[can12m == T, .SD, .SDcols = c(vars, "weights", "regsample")]
subtab_svy <- svydesign(ids = ~1, weights = ~weights, data = subtab)

# Total sample table
table1a <- tbl_svysummary(
  subtab_svy,
  include = vars,
  statistic = list(
    all_continuous() ~ "{mean} ({sd})",
    all_categorical() ~ "{n_unweighted} ({p}%)"
  ),
  digits = list(all_continuous() ~ 2, all_categorical() ~ c(0, 1))
)

# Subgroup table with custom tests
table1b <- tbl_svysummary(
  subtab_svy,
  by = regsample,
  include = vars,
  statistic = list(
    all_continuous() ~ "{mean} ({sd})",
    all_categorical() ~ "{n_unweighted} ({p}%)"
  ),
  digits = list(all_continuous() ~ 2, all_categorical() ~ c(0, 1))
) %>%
  add_p(
    test = list(
      all_continuous() ~ "custom_ttest",
      all_categorical() ~ "custom_chisq"
    ),
    pvalue_fun = ~ style_pvalue(., digits = 3)
  )

# Merge and export
combined_table <- tbl_merge(
  tbls = list(table1a, table1b),
  tab_spanner = c("**Overall**", "**By Subgroup**")
)

combined_table %>%
  as_kable_extra() %>%
  kableExtra::kable_styling(
    full_width = FALSE,
    font_size = 16,
    html_font = "Arial"
  ) %>%
  kableExtra::save_kable(file = paste0("tabs/tab1_",DATE,".html"))

rm(subtab, subtab_svy, table1a, table1b, combined_table)

############




## 3.2) TABLE 2
#-------------------------------------------------------

# total sample
data[can12m == T, .(sum(source_homegrow),wtd.mean(source_homegrow, weights))] # 39.7%
data[can12m == T, .(sum(source_illegal), wtd.mean(source_illegal, weights))] # 38.4%
data[can12m == T, .(sum(source_social), wtd.mean(source_social, weights))] # 34.5%
data[can12m == T, .(sum(source_pharma), wtd.mean(source_pharma, weights))] # 17.8%
data[can12m == T, .(sum(source_other), wtd.mean(source_other, weights))] # 11.9%

# subsample
data[regsample == T, .(sum(source_homegrow),wtd.mean(source_homegrow, weights))] # 39.7%
data[regsample == T, .(sum(source_illegal), wtd.mean(source_illegal, weights))] # 38.4%
data[regsample == T, .(sum(source_social), wtd.mean(source_social, weights))] # 34.5%
data[regsample == T, .(sum(source_pharma), wtd.mean(source_pharma, weights))] # 17.8%
data[regsample == T, .(sum(source_other), wtd.mean(source_other, weights))] # 11.9%

# STRICTLY LEGAL
data[can12m == T, .(sum(only_homegrow == T | only_cannabisclub == T | only_medical == T),
                    wtd.mean(only_homegrow == T | only_cannabisclub == T | only_medical == T, weights))] # 16.7%
data[regsample == T, .(sum(only_homegrow == T | only_cannabisclub == T | only_medical == T),
                       wtd.mean(only_homegrow == T | only_cannabisclub == T | only_medical == T, weights))] # 18.6%




# ==================================================================================================================================================================
# ==================================================================================================================================================================
# ==================================================================================================================================================================

## 3) SUPPLEMENTARY TABLES
# ______________________________________________________________________________________________________________________

## 3.1) SUPPLEMENTARY TABLE 1
#-------------------------------------------------------

# Define variables and create survey design
vars <- names(data)[names(data) %like% "sex|age$|agegroup|edu|gisd|DEGURBA|can12m|alc12m|tob12m"]
subtab <- data[, .SD, .SDcols = c(vars, "weights", "can12m")]
names(subtab)[11] <- "subsample"
subtab_svy <- svydesign(ids = ~1, weights = ~weights, data = subtab)

# Total sample table
table1a <- tbl_svysummary(
  subtab_svy,
  include = vars,
  statistic = list(
    all_continuous() ~ "{mean} ({sd})",
    all_categorical() ~ "{n_unweighted} ({p}%)"
  ),
  digits = list(all_continuous() ~ 2, all_categorical() ~ c(0, 1))
)

# Subgroup table (all 12m can users) with custom tests
table1b <- tbl_svysummary(
  subtab_svy,
  by = subsample,
  include = vars,
  statistic = list(
    all_continuous() ~ "{mean} ({sd})",
    all_categorical() ~ "{n_unweighted} ({p}%)"
  ),
  digits = list(all_continuous() ~ 2, all_categorical() ~ c(0, 1))
) %>%
  add_p(
    test = list(
      all_continuous() ~ "custom_ttest",
      all_categorical() ~ "custom_chisq"
    ),
    pvalue_fun = ~ style_pvalue(., digits = 3)
  )

# Merge and export
combined_table <- tbl_merge(
  tbls = list(table1a, table1b),
  tab_spanner = c("**Overall**", "**By Can-Use**")
)

combined_table %>%
  as_kable_extra() %>%
  kableExtra::kable_styling(
    full_width = FALSE,
    font_size = 16,
    html_font = "Arial"
  ) %>%
  kableExtra::save_kable(file = paste0("tabs/sup_tab1_",DATE,".html"))

rm(subtab, subtab_svy, table1a, table1b, combined_table)


# ==================================================================================================================================================================
# ==================================================================================================================================================================
# ==================================================================================================================================================================

## 4) FIGURES
# ______________________________________________________________________________________________________________________

# colors
color_2 <- c("#FF6666", "#66ffff")
color_4 <- c("#9b3721","#489b21","#21859b","#74219b")
color_6 <- c("#f39c12","#2ecc71","#16a085","#3498db","#8e44ad","#e74c3c")
color_5 <- c("#4e0a39", "#c2154a", "#eb8102", "#f8cb00", "#929d0e")
color_6 <- c("#4e0a39", "#c2154a", "#eb8102", "#f8cb00", "#929d0e", "#007c5f")
color_7 <- c("#4e0a39", "#c2154a", "#eb8102", "#f8cb00", "#929d0e", "#007c5f", "#004b7d")


## 4.1) LCA plot
#-------------------------------------------------------

# get class_props
class_props <- data[!is.na(class), .(weighted_sum = sum(weights)), by = class][
  , .(class,weighted_prop = weighted_sum / sum(weighted_sum))][order(class)]

# prepare data
pdat <- as.data.frame(lcamod$probs)
pdat <- melt(data.table(class  = rownames(pdat), pdat), id.vars = "class")
pdat <- pdat[stringr::str_sub(variable, -2,-2) == 2]
pdat[, class := stringr::str_sub(class, 7,7)]
pdat[, source := stringr::str_sub(variable, 3,-7)]
pdat[, prob := round(value,10)]
pdat <- pdat[,.(class,source,prob)]
pdat$source <- factor(pdat$source, c("knowndealer","unknowns","online","socialmedia",
                                     "friends", 
                                     "cannabisclub",
                                     "homegrown_self", "homegrown_others", 
                                     "pharmacy",
                                     "no_possess", "other"))

pdat$class <- factor(pdat$class, levels = 1:6, 
                     labels = c(
                       "HOME CULTIVATION",
                       "NO-POSESS",
                       "ILLEGAL",
                       "MIX",
                       "SOCIAL",
                       "PHARMACY"))

#class_p <- lcamod_all$P
pdat$class_lab <- NA_character_
for(c in 1:length(class_props$class)){
  
  lab <- levels(pdat$class)[c]
  pdat[class == levels(pdat$class)[c], 
       ':=' (class_lab = paste0(lab,
                                " (",
                                round(class_props[class == lab]$weighted_prop,3)*100,"%)"),
             class_prop = class_props[class == lab]$weighted_prop)]    
  rm(lab)
  
}

pdat$class_lab
levelorder <- unique(pdat[order(class_prop, decreasing = T)]$class_lab)
pdat$class_lab <- factor(pdat$class_lab, 
                         levels = levelorder)

# generate plot
ggplot(pdat, aes(x = source, y = prob, group = class_lab, color = class_lab)) +
  geom_point() +
  geom_line(linewidth=1) +
  theme(
    axis.text.x = element_text(hjust = 0.5, vjust = 0.5, angle = 45),
    axis.text.y = element_text(size = 12),              
    axis.title.x = element_text(size = 14),          
    axis.title.y = element_text(size = 14),            
    plot.title = element_text(hjust = 0.5, size = 16), 
    legend.text = element_text(size = 12),             
    legend.title = element_text(size = 14),
    legend.position = "bottom",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank()) +
  guides(color = guide_legend(nrow = 2), linetype = guide_legend(nrow = 1)) +
  scale_x_discrete("") +
  scale_color_manual("classes (weighted population share)", values = color_6) +
  scale_y_continuous("% Response probability\nwithin each class", 
                     breaks = seq(0, 1, by = 0.1), limits = c(0, 1), labels = scales::percent) 

ggsave(paste0("figures/Fig1_LCA_probabilities_",DATE,".png"), height = 6, width = 12)
ggsave(paste0("figures/Fig1_LCA_probabilities_",DATE,".pdf"), height = 6, width = 12)

rm(pdat, c, levelorder)

## 4.2) Forest plot
#-------------------------------------------------------

pdat <- regout[variable != "(Intercept)",.(class,varlab,riskratio,lower = `2.5 %`,upper = `97.5 %`)]
pdat$class <- factor(pdat$class, levels = rev(levels(data$class)))

pdat$varlab <- factor(pdat$varlab, levels = rev(levels(pdat$varlab)))

pdat[, group := ifelse(varlab %like% "Sex|Age|Edu|depri|DEGURBA", "sociodemographics",
                       ifelse(varlab %like% "frequenc|purpose|prescription|CAST", "cannabis-related", "health-related"))]
pdat$group <- factor(pdat$group, levels = c("sociodemographics","cannabis-related", "health-related"))

ggplot(pdat, aes(x = varlab, y = riskratio, color = class)) + 
  facet_grid(group ~ ., scales = "free", space = "free") +
  geom_hline(yintercept = 1) +
  geom_point(position = position_dodge(0.7), size = 2) +
  geom_errorbar(aes(ymin = lower, ymax = upper), 
                width = 0, position = position_dodge(0.7), linewidth = 0.4) +
  scale_x_discrete("") +
  scale_y_continuous("Risk Ratio", transform = "log10") +
  scale_color_manual("", values = rev(color_6[2:6]),
                     guide = guide_legend(reverse = TRUE)) +
  coord_flip() 

ggsave(paste0("figures/Fig2_Forest plot_",DATE,".png"), height = 8, width = 10)
ggsave(paste0("figures/Fig2_Forest plot_",DATE,".pdf"), height = 8, width = 10)

rm(pdat)


## 4.3) Jitter
#-------------------------------------------------------

pdat <- data[regsample == T, .(class, quant_tot, weights)]
pdat[, weight := sum(weights), by = class]

ggplot(pdat, aes(x = class, y = quant_tot, fill = class)) + 
  geom_boxplot(alpha = 0.8, show.legend = F) + 
  geom_jitter(alpha = 0.3, show.legend = F, width = 0.2) + 
  scale_x_discrete("") + 
  scale_y_continuous("30-day cannabis use quanties\n(logarithmized scale)", trans = "log10") + 
  scale_fill_manual(values = color_6)

ggsave(paste0("figures/Fig3_QUANT JITTER_",DATE,".png"), height = 5, width = 10)
ggsave(paste0("figures/Fig3_QUANT JITTER_",DATE,".pdf"), height = 5, width = 10)

rm(pdat)

# ==================================================================================================================================================================
# ==================================================================================================================================================================
# ==================================================================================================================================================================

## 5) SUPPLEMENTAL FIGURES
# ______________________________________________________________________________________________________________________

## 5.1) weights
#-------------------------------------------------------

ggplot(data[!is.na(edu)], aes(x = sex, y = weights)) +
  facet_grid(edu ~ agegroup) + 
  geom_hline(yintercept = 1) +
  geom_jitter(alpha = 0.1, show.legend = F, width = 0.2)

ggsave(paste0("figures/Supp Fig1_WEIGHTS_",DATE,".png"), height = 5, width = 10)
