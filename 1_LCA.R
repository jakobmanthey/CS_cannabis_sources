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

## 2) LCA
# ______________________________________________________________________________________________________________________

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
  save_kable(file = paste0("tabs/sup_tab4_",DATE,".html"))


##  save plots for all candidate models
for (m in 2:8){
  
  mod <- get(ls()[ls() %like% paste0("^am",m)])
  png(filename = paste0("figures/Supp Fig",m,"_LCA_mod_",m,"_classes_",DATE,".png"),
      width = 12, height = 6, unit = "in", res = 300)
  plot(mod)  
  dev.off()
}
 
  
##  add class findings to data
lcamod <- am6

##  add classes
#data$class <- NULL
data$class <- NA_character_
data[can12m == T]$class <- lcamod$predclass

data$class <- factor(data$class, levels = 1:6, 
                     labels = c(
                       "HOME CULTIVATION",
                       "OTHERS' SUPPLY",
                       "ILLEGAL",
                       "MIX",
                       "SOCIAL",
                       "PHARMACY"))

data$class <- factor(data$class,
                     levels = c(
                       "MIX",
                       "SOCIAL",
                       "ILLEGAL",
                       "HOME CULTIVATION",
                       "PHARMACY",
                       "OTHERS' SUPPLY"))

data[, prop.table(table(class))]

saveRDS(data, paste0("data/data_sources_prepared_after LCA_",DATE,".rds"))
saveRDS(lcamod, paste0("data/LCA model_",DATE,".rds"))


# ==================================================================================================================================================================
# ==================================================================================================================================================================
# ==================================================================================================================================================================