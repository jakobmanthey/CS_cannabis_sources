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

table(rowSums(data[can12m == T, .(source_illegal, source_social, 
                            source_homegrow, source_pharma, source_other)]))
# (378+92+4)/1490 = 31.8%

data[regsample == T, mean(source_illegal)] # 39.4%
data[regsample == T, mean(source_social)] # 35.7%
data[regsample == T, mean(source_homegrow)] # 48.3%
data[regsample == T, mean(source_pharma)] # 19.8%
data[regsample == T, mean(source_other)] # 5.8%

## 2.3) LCA
#-------------------------------------------------------

library(poLCA)

f <- as.formula(paste0("cbind(",paste(
  c("s_medical", "s_friends", "s_unknowns", "s_knowndealer", "s_homegrown_self", "s_homegrown_others", "s_cannabisclub",
    "s_online", "s_socialmedia", "s_no_possess", "s_other"), collapse = ","), ")~1" ))
set.seed(341)
am2 <- poLCA(f, data = data, nclass = 2, maxiter = 100000, nrep = 50)
am3 <- poLCA(f, data = data, nclass = 3, maxiter = 100000, nrep = 50)
am4 <- poLCA(f, data = data, nclass = 4, maxiter = 100000, nrep = 50)
am5 <- poLCA(f, data = data, nclass = 5, maxiter = 100000, nrep = 50)
am6 <- poLCA(f, data = data, nclass = 6, maxiter = 100000, nrep = 50)
am7 <- poLCA(f, data = data, nclass = 7, maxiter = 100000, nrep = 50)
am8 <- poLCA(f, data = data, nclass = 8, maxiter = 100000, nrep = 50)

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
                                get_entropy(am8)))
sumtab %>% knitr::kable()

for (m in 2:8){
  
  mod <- get(ls()[ls() %like% paste0("^am",m)])
  png(filename = paste0("figures/LCA_mod_",m,"_classes_",DATE,".png"),
      width = 12, height = 6, unit = "in", res = 300)
  plot(mod)  
  dev.off()
  
}

lcamod_all <- am5


## 2.4) Cannabis sourcing and use quantities
#-------------------------------------------------------

data$class <- NULL
data$class <- NA_character_
data[can12m == T]$class <- lcamod_all$predclass

data$class <- factor(data$class, levels = 1:5, 
                     labels = c(
                       "NO-POSESS",
                       "ILLEGAL",
                       "HOMEGROW",
                       "SOCIAL",
                       "MEDICAL"))

data$class <- factor(data$class,
                     levels = c(
                       "ILLEGAL",
                       "SOCIAL",
                       "HOMEGROW",
                       "MEDICAL",
                       "NO-POSESS"))

data[regsample == T, .(.N,prop = .N/nrow(data[regsample == T])), by = class]
data[regsample == T, mean(quant_tot), by = class][order(class)]
data[regsample == T, median(quant_tot), by = class][order(class)]

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
unique(temp[, .(class,prop = wt_class/wt_tot)])





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
color_6 <- c("#f39c12","#2ecc71","#16a085","#3498db","#8e44ad","#e74c3c")
color_5 <- c("#4e0a39", "#c2154a", "#eb8102", "#f8cb00", "#929d0e")


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
                                     "medical",
                                     "no_possess", "other"))

pdat$class <- factor(pdat$class, levels = 1:5, 
                     labels = c(
                       "NO-POSESS",
                       "ILLEGAL",
                       "HOMEGROW",
                       "SOCIAL",
                       "MEDICAL"))

class_p <- lcamod_all$P
pdat$class_lab <- NA_character_
for(c in 1:length(class_p)){
  pdat[class == levels(pdat$class)[c], 
       ':=' (class_lab = paste0(levels(pdat$class)[c],
                           " (",
                           round(class_p[c],3)*100,"%)"),
             class_prop = class_p[c])]    
  
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
  guides(color = guide_legend(nrow = 1), linetype = guide_legend(nrow = 1)) +
  scale_x_discrete("") +
  scale_color_manual("classes (population share)", values = color_5) +
  #scale_linetype_manual(values = linetype_values, labels = label_values) +
  scale_y_continuous("% Response probability within each class", breaks = seq(0, 1, by = 0.1), limits = c(0, 1), labels = scales::percent) 

ggsave(paste0("figures/LCA_probabilities_",DATE,".png"), height = 6, width = 10)

## 3.2) Jitter
#-------------------------------------------------------

pdat <- data[regsample == T, .(class, quant_tot)]

ggplot(pdat, aes(x = class, y = quant_tot)) + 
  geom_jitter(alpha = 0.5) + 
  geom_boxplot(alpha = 0.5)

ggsave(paste0("figures/QUANT JITTER_",DATE,".png"), height = 6, width = 10)

