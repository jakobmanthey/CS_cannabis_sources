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





# ==================================================================================================================================================================
# ==================================================================================================================================================================
# ==================================================================================================================================================================

# 1) PREPARE DATA
# ______________________________________________________________________________________________________________________

##  run source code
# -------------------------------------------------------

source("0_prepare.R")

input[, min(START)]
input[, max(END)]


##  some helper functions
# -------------------------------------------------------

##  entropy from stack overflow or this ppt: // https://daob.nl/wp-content/uploads/2015/07/ESRA-course-slides.pdf
get_entropy<-function(mod){
  entropy<-function(p)sum(-p*log(p),na.rm = T)
  error_prior <- entropy(mod$P) # Class proportions
  error_post <- mean(apply(mod$posterior, 1, entropy))
  R2_entropy <- (error_prior - error_post) / error_prior
  return(R2_entropy)
}


calculate_lca_metrics <- function(model, num_classes, return_alcpp = TRUE) {
  # Extract posterior probabilities and predicted classes
  posteriors <- as.data.table(model$posterior)
  posteriors[, predclass := model$predclass]
  
  # Ensure correct column order
  class_cols <- paste0("class", 1:num_classes)
  names(posteriors) <- c(class_cols, "predclass")
  setcolorder(posteriors, c(class_cols, "predclass"))
  
  # Calculate the classification table
  classification_table <- posteriors[, lapply(.SD, sum), 
                                     by = predclass, 
                                     .SDcols = class_cols]
  setorder(classification_table, predclass)
  
  # Calculate ALCPP
  alcpp_matrix <- matrix(0, nrow = num_classes, ncol = num_classes)
  for (i in 1:num_classes) {
    class_members <- posteriors[predclass == i]
    alcpp_matrix[i,] <- colMeans(class_members[, ..class_cols])
  }
  
  # Calculate lowest ALCPP diagonal value
  lowest_alcpp <- min(diag(alcpp_matrix))
  
  # Return based on user choice
  if (return_alcpp) {
    return(alcpp_matrix)  # Return full ALCPP matrix
  } else {
    return(lowest_alcpp)  # Return lowest diagonal value of ALCPP
  }
}




##########

## Average posterior class membership probabilities
calculate_classification_metrics <- function(model, num_classes, lowest) {
  # Extract the posterior probabilities
  class_probs <- data.frame(model$posterior)
  
  # Calculate the maximum probability for each participant
  class_probs$max_prob <- apply(class_probs, 1, max)
  
  # Calculate average posterior probability
  avg_posterior_prob <- mean(class_probs$max_prob)
  cat("Average posterior probability: ", avg_posterior_prob, "\n")
  
  # Create a dataframe with posterior probabilities and predicted classes
  posteriors <- data.table(model$posterior, model$predclass)
  
  # Rename the predclass column
  names(posteriors)[num_classes +1] <- "predclass"
  
  # Calculate the classification table by summing the probabilities
  classification_table <- posteriors[, lapply(.SD, function(x) {round(sum(x),3)}), by = predclass, .SDcols = names(posteriors)[names(posteriors) != "predclass"]][order(predclass)]
  
  # Transform the counts to proportions
  classification_table_prop <- copy(classification_table)
  select <- names(classification_table)[names(classification_table) != "predclass"]
  classification_table_prop[, (select) := lapply(.SD, function(x) x/sum(x)), .SDcols = select]
  rm(select)
  
  # Print the proportions table
  print(classification_table_prop)
  
  # Convert to matrix and calculate classification error
  classification_values_only <- as.matrix(classification_table[,2:(num_classes + 1)])
  total_classification_error <- 1 - sum(diag(classification_values_only)) / sum(classification_values_only)
  
  cat("Total classification error: ", total_classification_error, "\n")
  
  # REPORT
  if (lowest == T){
    # report lowest:
    
    select <- names(classification_table_prop)[names(classification_table_prop) != "predclass"]
    temp <- copy(classification_table_prop[,.SD,.SDcols = select])
    return(min(diag(as.matrix(temp))))
    
  } else {
    
    # report results as a list for potential further use
    return(list(avg_posterior_prob = avg_posterior_prob,
                classification_table = classification_table,
                classification_table_prop = classification_table_prop,
                total_classification_error = total_classification_error))
  }
}
  
  

# ==================================================================================================================================================================
# ==================================================================================================================================================================
# ==================================================================================================================================================================

## 2) RESULTS
# ______________________________________________________________________________________________________________________

## 2.1) Descriptives general
#-------------------------------------------------------

data[, table(can12m)] # 1478
data[, table(can30d)] # 1073

## 2.2) Frequency of Sources
#-------------------------------------------------------

data[can12m == T, wtd.mean(source_homegrow, weights)] # 39.7%
data[can12m == T, wtd.mean(source_illegal, weights)] # 38.4%
data[can12m == T, wtd.mean(source_social, weights)] # 34.5%
data[can12m == T, wtd.mean(source_pharma, weights)] # 17.8%
data[can12m == T, wtd.mean(source_other, weights)] # 11.9%

table(rowSums(data[can12m == T, .(source_illegal, source_social, 
                            source_homegrow, source_pharma, source_other)]))
# (378+92+4)/1478 = 32.1%

data[can12m == T, wtd.mean(source_homegrow == T | source_pharma == T, weights)] # 50.7%


# STRICTLY LEGAL
data[can12m == T, mean(only_homegrow == T)]
data[can12m == T, wtd.mean(only_homegrow == T, weights)] # 7.4%
data[can12m == T, mean(only_medical == T)]
data[can12m == T, wtd.mean(only_medical == T, weights)] # 7.2%
data[can12m == T, mean(only_homegrow_or_medical == T)]
data[can12m == T, wtd.mean(only_homegrow_or_medical == T, weights)] # 1.3%

# subsample
data[regsample == T, wtd.mean(source_illegal, weights)] # 41.7%
data[regsample == T, wtd.mean(source_social, weights)] # 35.2%
data[regsample == T, wtd.mean(source_homegrow, weights)] # 49.1%
data[regsample == T, wtd.mean(source_pharma, weights)] # 21.3%
data[regsample == T, wtd.mean(source_other, weights)] # 5.3%


## 2.3) LCA
#-------------------------------------------------------

library(poLCA)

f <- as.formula(paste0("cbind(",paste(
  c("s_pharmacy", "s_friends", "s_unknowns", "s_knowndealer", "s_homegrown_self", "s_homegrown_others", "s_cannabisclub",
    "s_online", "s_socialmedia", "s_no_possess", "s_other"), collapse = ","), ")~1" ))
set.seed(341)
am2 <- poLCA(f, data = data, nclass = 2, maxiter = 100000, nrep = 50)
am3 <- poLCA(f, data = data, nclass = 3, maxiter = 100000, nrep = 50)
am4 <- poLCA(f, data = data, nclass = 4, maxiter = 100000, nrep = 50)
am5 <- poLCA(f, data = data, nclass = 5, maxiter = 100000, nrep = 50)
am6 <- poLCA(f, data = data, nclass = 6, maxiter = 100000, nrep = 50)
am7 <- poLCA(f, data = data, nclass = 7, maxiter = 100000, nrep = 50)
am8 <- poLCA(f, data = data, nclass = 8, maxiter = 100000, nrep = 50)


##  get model fit indices (Supp Table 1)
sumtab = data.frame(Models = stringr::str_c(c(2,3,4,5,6,7,8), " Class"),
                      LL = c(am2$llik, am3$llik, am4$llik, am5$llik, am6$llik, am7$llik, am8$llik),
                      AIC = c(am2$aic, am3$aic, am4$aic, am5$aic, am6$aic, am7$aic, am8$aic),
                      BIC = c(am2$bic, am3$bic, am4$bic, am5$bic, am6$bic, am7$bic, am8$bic), 
                      Entropy = c(get_entropy(am2), 
                                  get_entropy(am3),
                                  get_entropy(am4),
                                  get_entropy(am5),
                                  get_entropy(am6),
                                  get_entropy(am7),
                                  get_entropy(am8)),
                      smallest_n = c(min(table(am2$predclass)),
                                     min(table(am3$predclass)),
                                     min(table(am4$predclass)),
                                     min(table(am5$predclass)),
                                     min(table(am6$predclass)),
                                     min(table(am7$predclass)),
                                     min(table(am8$predclass))),
                      ALCPP_lowest = c(calculate_lca_metrics(am2, 2, return_alcpp = F),
                                       calculate_lca_metrics(am3, 3, return_alcpp = F),
                                       calculate_lca_metrics(am4, 4, return_alcpp = F),
                                       calculate_lca_metrics(am5, 5, return_alcpp = F),
                                       calculate_lca_metrics(am6, 6, return_alcpp = F),
                                       calculate_lca_metrics(am7, 7, return_alcpp = F),
                                       calculate_lca_metrics(am8, 8, return_alcpp = F)))
sumtab %>% knitr::kable()
kable(sumtab) %>%
  kable_styling() %>%
  save_kable(file = "tabs/sup_tab1.html")


##  save plots for all candidate models
for (m in 2:8){
  
  mod <- get(ls()[ls() %like% paste0("^am",m)])
  png(filename = paste0("figures/Supp Fig",m-1,"_LCA_mod_",m,"_classes_",DATE,".png"),
      width = 12, height = 6, unit = "in", res = 300)
  plot(mod)  
  dev.off()
}
 
  
##  add class findings to data
lcamod_all <- am6

##  add classes
#data$class <- NULL
data$class <- NA_character_
data[can12m == T]$class <- lcamod_all$predclass

data$class <- factor(data$class, levels = 1:6, 
                     labels = c(
                       "HOMEGROW",
                       "NO-POSESS",
                       "ILLEGAL",
                       "MIX",
                       "SOCIAL",
                       "MEDICAL"))

data$class <- factor(data$class,
                     levels = c(
                       "MIX",
                       "SOCIAL",
                       "ILLEGAL",
                       "HOMEGROW",
                       "MEDICAL",
                       "NO-POSESS"))

data[, prop.table(table(class))]
class_props <- data[!is.na(class), .(weighted_sum = sum(weights)), by = class][
  , .(class,weighted_prop = weighted_sum / sum(weighted_sum))][order(class)]


## 2.4) Multinomial regression
#-------------------------------------------------------

vars <- names(data)[names(data) %like% "sex|agegr|edu|DEGURBA|gisd|distress|health|freq12m|purpose|prescribed|castrisk"]
f <- as.formula(paste0("class ~ ", paste0(vars, collapse = " + ")))

sub <- copy(data[can12m == T & sex != "other"])
sub$sex <-factor(sub$sex)
sub$edu <-relevel(sub$edu, ref = "mid")

mnmod <- multinom(formula = f,data = sub)
summary(mnmod)

#tab2 <- gtsummary::tbl_regression(mnmod, exponentiate = TRUE)

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

var_levels <- unique(regout$variable)
var_levels <- c(var_levels[var_levels %like% "Intercept"],
                var_levels[var_levels %like% "sex"],
                var_levels[var_levels %like% "agegr"],
                var_levels[var_levels %like% "edu"],
                var_levels[var_levels %like% "gisd"],
                var_levels[var_levels %like% "DEGURBA"],
                var_levels[var_levels %like% "monthly"],
                var_levels[var_levels %like% "weekly"],
                var_levels[var_levels %like% "daily"],
                var_levels[var_levels %like% "purpose"],
                var_levels[var_levels %like% "prescr"],
                var_levels[var_levels %like% "cast"],
                var_levels[var_levels %like% "health"],
                var_levels[var_levels %like% "distress"])

regout$variable <- factor(regout$variable, levels = var_levels)

kable(regout[order(class,variable),.(class,variable,riskratio,ci,p)]) %>%
  kable_styling() %>%
  save_kable(file = "tabs/sup_tab2.html")


## 2.5) Cannabis sourcing and use quantities
#-------------------------------------------------------

# distribution of classes in subsample
#data[regsample == T, prop.table(table(class))]
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

#social+medical+homegrow:
temp[class %like% "SOCIAL|MEDICAL|HOMEGROW", sum(prop)]

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

vars <- names(data)[names(data) %like% "sex|age|edu|gisd|DEGURBA|freq12m|freq30d|quant_pd|earlyonset|purpose|prescribed|castrisk|distress|health"]

table1a <- CreateTableOne(data = data[can12m == T,.SD, .SDcols = vars])
table1a_df <- as.data.frame(print(table1a))
#table1a_df <- as.data.table(print(table1a))

table1b <- CreateTableOne(data = data[regsample == T,.SD, .SDcols = vars])
table1b_df <- as.data.frame(print(table1b))
#table1b_df <- as.data.table(print(table1b))

out <- cbind(table1a_df, table1b_df)

kable(out) %>%
  kable_styling() %>%
  save_kable(file = "tabs/tab1.html")

## BIAS OF REGRESSION SAMPLE (TABLE 1)

### comparison:
summary(glm(age ~ regsample,family = "gaussian", data[can12m == T])) # no diff
summary(glm(sex == "men" ~ regsample,family = "binomial", data[can12m == T])) # more men
summary(glm(edu == "mid" ~ regsample,family = "binomial", data[can12m == T])) # more mid edu
summary(glm(edu == "high" ~ regsample,family = "binomial", data[can12m == T])) # more high edu
summary(glm(gisd_k == "high" ~ regsample,family = "binomial", data[can12m == T])) # no diff
summary(glm(gisd_k == "mid" ~ regsample,family = "binomial", data[can12m == T])) # no diff
summary(glm(DEGURBA == "cities" ~ regsample,family = "binomial", data[can12m == T])) # no diff
summary(glm(DEGURBA == "towns_suburbs" ~ regsample,family = "binomial", data[can12m == T])) # no diff
summary(glm(freq30d ~ regsample,family = "gaussian", data[can12m == T])) # more use days
summary(glm(earlyonset ~ regsample,family = "binomial", data[can12m == T])) # earlier onset
summary(glm(purpose == "only medical" ~ regsample,family = "binomial", data[can12m == T])) # less only medical 
summary(glm(prescribed ~ regsample,family = "binomial", data[can12m == T])) # more prescription
summary(glm(castrisk ~ regsample,family = "binomial", data[can12m == T])) # more castrisk
summary(glm(distress ~ regsample,family = "binomial", data[can12m == T])) # no diff
summary(glm(health_good ~ regsample,family = "binomial", data[can12m == T])) # no diff
summary(glm(health_chronic ~ regsample,family = "binomial", data[can12m == T])) # no diff



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


## 3.1) LCA plot
#-------------------------------------------------------

# prepare data
pdat <- as.data.frame(lcamod_all$probs)
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
                       "HOMEGROW",
                       "NO-POSESS",
                       "ILLEGAL",
                       "MIX",
                       "SOCIAL",
                       "MEDICAL"))

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
  #theme_minimal() +
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
  #scale_linetype_manual(values = linetype_values, labels = label_values) +
  scale_y_continuous("% Response probability\nwithin each class", 
                     breaks = seq(0, 1, by = 0.1), limits = c(0, 1), labels = scales::percent) 

ggsave(paste0("figures/Fig1_LCA_probabilities_",DATE,".png"), height = 6, width = 12)

rm(pdat, c, lab, levelorder)

## 3.2) Forest plot
#-------------------------------------------------------

pdat <- regout[variable != "(Intercept)",.(class,variable,riskratio,lower = `2.5 %`,upper = `97.5 %`)]
pdat$class <- factor(pdat$class, levels = rev(levels(data$class)))

pdat$variable <- factor(pdat$variable, levels = rev(levels(pdat$variable)))

pdat[, group := ifelse(variable %like% "sex|agegr|edu|gisd|DEGURBA", "sociodemographics",
                       ifelse(variable %like% "freq|purpose|prescribed|castrisk", "cannabis-related", "health-related"))]
pdat$group <- factor(pdat$group, levels = c("sociodemographics","cannabis-related", "health-related"))

ggplot(pdat, aes(x = variable, y = riskratio, color = class)) + 
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

rm(pdat)


## 3.3) Jitter
#-------------------------------------------------------

pdat <- data[regsample == T, .(class, quant_tot, weights)]
pdat[, weight := sum(weights), by = class]

ggplot(pdat, aes(x = class, y = quant_tot, fill = class)) + 
  geom_boxplot(alpha = 0.8, show.legend = F) + 
  geom_jitter(alpha = 0.3, show.legend = F, width = 0.2) + 
  #geom_violin(alpha = 0.8, show.legend = F) + 
  scale_x_discrete("") + 
  scale_y_continuous("30-day cannabis use quanties\n(logarithmized scale)", trans = "log10") + 
  scale_fill_manual(values = color_6)

ggsave(paste0("figures/Fig3_QUANT JITTER_",DATE,".png"), height = 5, width = 8)

