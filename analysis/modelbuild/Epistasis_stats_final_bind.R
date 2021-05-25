library(MASS)
library(gtools)
library(stringi)
library(comprehenr)
library(AICcmodavg)
library(rlist)

#Set source for stepaicc module and data path (need to modify to match repo location)
source("/path/to/repo/EpistasisStats/analysis/stepaicc/stepAICc.r")
setwd("/path/to/repo/EpistasisStats/data/processed/")

#Read in data, define features (non abstracted)
data<-read.table("skempi_bind_processed.txt", header=T, sep='\t')
features <- list("separation", "cplx_type", "intside_int", "size_net")
dep_var <- "ep_abs"

#Generate list of abstractions to build permutated lists for abstracted features
hp_abs <- list("1" = 'hp_ab1', "2" = 'hp_ab2', '3'= 'hp_ab3')
charge_abs <- list("1" = 'charge_ab1', '2' = 'charge_ab2', 
                   '3' = 'charge_ab3')
ss_abs <- list('1' = 'ss_ab1', '2' = 'ss_ab2', '3' = 'ss_ab3')
sasa_abs <- list('1' = 'rel_sasa_net', '2' = 'abs_sasa_net')

#Set the list of features which are continuous or do not have abstractions
staticlist <- paste(features, collapse=" + ")

#Build list of all combinations of abstractions
abstractionslist <- as.list(outer(as.list(outer(as.list(outer(charge_abs, 
                    ss_abs, paste, sep=" + ")), hp_abs, paste, sep=" + ")), sasa_abs, paste, sep=" + "))
allcombos <- as.list(outer(staticlist, abstractionslist, paste, sep=" + "))


#Function to randomize features and find best model using AICc, returns the AICc for the model
randomizedAIC <- function(l) {
  randomized = stri_join_list(as.list(sample(stri_split(l, fixed=" + ")[[1]])), collapse=" + ")
  fullform = paste("ep_abs ~ ", paste(paste("(",randomized), ")^2"))
  rfit <- glm(fullform, data=data)  
  fitStep <- stepAICc(rfit, direction='both',steps=100000)
  return(AICc(fitStep))
}

#Function to randomize features and find best model using corrected AICc, returns the model
randomizedfit <- function(l) {
  randomized = stri_join_list(as.list(sample(stri_split(l, fixed=" + ")[[1]])), collapse=" + ")
  fullform = paste("ep_abs ~ ", paste(paste("(",randomized), ")^2"))
  rfit <- glm(fullform, data=data)  
  fitStep <- stepAICc(rfit, direction='both',steps=100000)
  return(fitStep)
}

#Function to calculate the delta r-squared for all terms in a model. 
delta_r_sq <- function(fullmodel) {
  fmrsq <- summary(fullmodel)$r.squared
  drsq <- list(paste("Full model: ",fmrsq))
  terms <-  stri_split(c(bestmodel$terms[[3]]), fixed=" + ")[[1]]
  singles <- to_list(for (elem in terms) if (grepl(":", elem) != TRUE) c(elem))
  for (elem in singles) {
    sublist <- stri_join_list(to_list(for (t_elem in terms) if (grepl(elem, t_elem) != TRUE) c(t_elem)), 
                              collapse=" + ")
    nrsq <- summary(lm(paste("ep_abs ~ ", sublist), data=data))$r.squared
    drsq <- list.append(drsq, paste(elem,abs(nrsq - fmrsq),sep=": "))
  }
  return(drsq)
}

#Generates two sets with random order to correct for any effects of ordering, errors if there is a mismatch.
run1 <- to_list(for (elem in allcombos) c(randomizedAIC(elem)))
run2 <- to_list(for (elem in allcombos) c(randomizedAIC(elem)))
if (which.min(run1) != which.min(run2)) 
  {print("ERROR: Mismatched models")}

#Save results
bestmodel <- lm(randomizedfit(allcombos[[which.min(run1)]]))
sink('../model_build_results/best_model_bind.txt')
summary(lm(bestmodel))
sink()
sink("../model_build_results/delta_rsq_best_bind.txt")
delta_r_sq(bestmodel)
sink()



