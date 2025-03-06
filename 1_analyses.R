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

# ==================================================================================================================================================================
# ==================================================================================================================================================================
# ==================================================================================================================================================================

## 2) RESULTS
# ______________________________________________________________________________________________________________________

## 2.1) Descriptives general
#-------------------------------------------------------

data[, table(can12m)] # 1490
data[, table(can30d)] # 1074


## 2.2) Frequency of Sources
#-------------------------------------------------------

data[can12m == T, mean(source_illegal)] # 35.2%
data[can12m == T, mean(source_social)] # 35.0%
data[can12m == T, mean(source_homegrow)] # 38.5%
data[can12m == T, mean(source_pharma)] # 16.6%
data[can12m == T, mean(source_other)] # 13.1%

data[regsample == T, mean(source_illegal)] # 39.4%
data[regsample == T, mean(source_social)] # 35.7%
data[regsample == T, mean(source_homegrow)] # 48.3%
data[regsample == T, mean(source_pharma)] # 19.8%
data[regsample == T, mean(source_other)] # 5.8%

## 2.3) Regression
#-------------------------------------------------------

mod1 <- glm(source_illegal ~ sex + agegroup + edu + freq30d_group + quant_pd_group + purpose + castrisk, family = "binomial", data = data[regsample==T])
mod2 <- glm(source_social ~ sex + agegroup + edu + freq30d_group + quant_pd_group + purpose + castrisk, family = "binomial", data = data[regsample==T])
mod3 <- glm(source_homegrow ~ sex + agegroup + edu + freq30d_group + quant_pd_group + purpose + castrisk, family = "binomial", data = data[regsample==T])
mod4 <- glm(source_pharma ~ sex + agegroup + edu + freq30d_group + quant_pd_group + purpose + castrisk, family = "binomial", data = data[regsample==T])
mod5 <- glm(source_other ~ sex + agegroup + edu + freq30d_group + quant_pd_group + purpose + castrisk, family = "binomial", data = data[regsample==T])

summary(mod1)
summary(mod2)
summary(mod3)
summary(mod4)
summary(mod5)


##  LCA
library(poLCA)
sub <- data[regsample==T & sex != "other"]
#sub <- data[sex != "other"]
sub$sex <- factor(sub$sex)
sub[, edu2g := factor(ifelse(edu %like% "low|mid", "lowmid", "high"))]
sub[, table(edu,edu2g)]

##OLD
lca_mod2 <- poLCA(cbind(source_illegal_num, source_social_num, source_homegrow_num, source_pharma_num, source_other_num) ~ 
                    #sex + agegroup + edu + freq30d_group + quant_pd_group + purpose + castrisk,
                    sex + age + edu2g + quant_group + purpose + castrisk, 
                  data = sub, 
                  nclass = 2)
lca_mod3 <- poLCA(cbind(source_illegal_num, source_social_num, source_homegrow_num, source_pharma_num, source_other_num) ~ 
                    #sex + agegroup + edu + freq30d_group + quant_pd_group + purpose + castrisk,
                    sex + age + edu2g + quant_group + purpose + castrisk, 
                  data = sub, 
                  nclass = 3)
lca_mod4 <- poLCA(cbind(source_illegal_num, source_social_num, source_homegrow_num, source_pharma_num, source_other_num) ~ 
                    #sex + agegroup + edu + freq30d_group + quant_pd_group + purpose + castrisk,
                    sex + agegroup + edu2g + quant_group + purpose + castrisk, 
                  data = sub, 
                  nclass = 4)

##NEW  no covariates
f <- as.formula(paste0("cbind(",paste(
  c("source_illegal_num", "source_social_num", "source_homegrow_num", "source_pharma_num", "source_other_num")
  , collapse = ","), ")~1" ))
f <- as.formula(paste0("cbind(",paste(
  c("s_apo", "s_friends", "s_unbekannte", "s_dealer", "s_eigenanbau1", "s_eigenanbau2", "s_cav",
    "s_online", "s_socmed", "s_noposs", "s_other"), collapse = ","), ")~1" ))
set.seed(341)
m2 <- poLCA(f, data = sub, nclass = 2, maxiter = 100000, nrep = 50)
m3 <- poLCA(f, data = sub, nclass = 3, maxiter = 100000, nrep = 50)
m4 <- poLCA(f, data = sub, nclass = 4, maxiter = 100000, nrep = 50)
m5 <- poLCA(f, data = sub, nclass = 5, maxiter = 100000, nrep = 50)
m6 <- poLCA(f, data = sub, nclass = 6, maxiter = 100000, nrep = 50)
m7 <- poLCA(f, data = sub, nclass = 7, maxiter = 100000, nrep = 50)
m8 <- poLCA(f, data = sub, nclass = 8, maxiter = 100000, nrep = 50)


##  CHAT GPT
calculate_entropy <- function(lca_model) {
  # Extract posterior probabilities (individual class membership probabilities)
  post_probs <- lca_model$posterior  # This is an NxK matrix (N = number of observations, K = number of classes)
  
  # Calculate entropy for each individual
  entropy_individuals <- apply(post_probs, 1, function(p) {
    -sum(p * log(p), na.rm = TRUE)  # Apply the entropy formula
  })
  
  # Calculate the average entropy across all individuals
  total_entropy <- mean(entropy_individuals)
  
  return(total_entropy)
}

##  stack overflow // https://daob.nl/wp-content/uploads/2015/07/ESRA-course-slides.pdf
get_entropy<-function(mod){
  entropy<-function(p)sum(-p*log(p),na.rm = T)
  error_prior <- entropy(mod$P) # Class proportions
  error_post <- mean(apply(mod$posterior, 1, entropy))
  R2_entropy <- (error_prior - error_post) / error_prior
  return(R2_entropy)
  }



sumtab = data.frame(Models = stringr::str_c(c(2,3,4,5,6,7,8), " Class"),
                    LL = c(m2$llik, m3$llik, m4$llik, m5$llik, m6$llik, m7$llik, m8$llik),
                    AIC = c(m2$aic, m3$aic, m4$aic, m5$aic, m6$aic, m7$aic, m8$aic),
                    BIC = c(m2$bic, m3$bic, m4$bic, m5$bic, m6$bic, m7$bic, m8$bic), 
                    Entropy = c(get_entropy(m2), 
                                get_entropy(m3),
                                get_entropy(m4),
                                get_entropy(m5),
                                get_entropy(m6),
                                get_entropy(m7),
                                get_entropy(m8)))
sumtab %>% knitr::kable()



png("figures/lca_test.png", width = 10, height = 6, units = "in", res = 300)
par(mfrow = c(4, 2))
plot(m2)
plot(m3)
plot(m4)
plot(m5)
plot(m6)
plot(m7)
#plot(m8)
dev.off()

lcamod <- m4



# ==================================================================================================================================================================
# ==================================================================================================================================================================
# ==================================================================================================================================================================

## 3) TABLES
# ______________________________________________________________________________________________________________________

## 3.1) TABLE 1
#-------------------------------------------------------

library(tableone)
library(knitr)
library(kableExtra)

vars <- names(data)[names(data) %like% "sex|age|edu|freq12m|freq30d|quant_pd|purpose|castrisk"]

table1a <- CreateTableOne(data = data[can12m == T,.SD, .SDcols = vars])
table1a_df <- as.data.frame(print(table1a))

table1b <- CreateTableOne(data = data[regsample == T,.SD, .SDcols = vars])
table1b_df <- as.data.frame(print(table1b))

out <- cbind(table1a_df,table1b_df)

kable(out) %>%
  kable_styling() %>%
  save_kable(file = "tabs/tab1.html")

## BIAS OF REGRESSION SAMPLE (TABLE 1)

### comparison:
summary(glm(age ~ regsample,family = "gaussian", data[can12m == T])) # no diff
summary(glm(sex == "men" ~ regsample,family = "binomial", data[can12m == T])) # more men
summary(glm(edu == "mid" ~ regsample,family = "binomial", data[can12m == T])) # more mid edu
summary(glm(edu == "high" ~ regsample,family = "binomial", data[can12m == T])) # more high edu
summary(glm(freq30d ~ regsample,family = "gaussian", data[can12m == T])) # more use days
summary(glm(purpose == "only medical" ~ regsample,family = "binomial", data[can12m == T])) # less only medical 
summary(glm(castrisk ~ regsample,family = "binomial", data[can12m == T])) # more castrisk



# ==================================================================================================================================================================
# ==================================================================================================================================================================
# ==================================================================================================================================================================

## 4) FIGURES
# ______________________________________________________________________________________________________________________

# colors
color_2 <- c("#FF6666", "#66ffff")
color_4 <- c("#9b3721","#489b21","#21859b","#74219b")

## 3.1) LCA plot
#-------------------------------------------------------

# prepare data
pdat <- as.data.frame(lcamod$probs)
pdat <- melt(data.table(class  = rownames(pdat), pdat), id.vars = "class")
pdat <- pdat[stringr::str_sub(variable, -2,-2) == 2]
pdat[, class := stringr::str_sub(class, 7,7)]
pdat[, source := stringr::str_sub(variable, 3,-7)]
pdat[, prob := round(value,10)]
pdat <- pdat[,.(class,source,prob)]
pdat$source <- factor(pdat$source, c("apo", "friends", "unbekannte",
                                      "dealer", "online", "socmed",
                                      "eigenanbau1", "eigenanbau2", "cav",
                                       "noposs", "other"))

class_p <- lcamod$P
pdat$class_lab <- NA_character_
for(c in 1:length(class_p)){
  pdat[class == c, class_lab := paste0(c, " (",round(class_p[c],3)*100,"%)")]    
}


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
  guides(color = guide_legend(nrow = 1), linetype = guide_legend(nrow = 1)) +
  scale_x_discrete("") +
  scale_color_manual("classes (population share)", values = color_4) +
  #scale_linetype_manual(values = linetype_values, labels = label_values) +
  scale_y_continuous("% Response probability within each class", breaks = seq(0, 1, by = 0.1), limits = c(0, 1), labels = scales::percent) 

ggsave(paste0("figures/LCA_testfigure_",DATE,".png"), height = 6, width = 10)

## 3.2) Prevalence use by sex and edu
#-------------------------------------------------------


