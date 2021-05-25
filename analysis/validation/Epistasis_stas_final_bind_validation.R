library(MASS)
library(gtools)
library(stringi)
library(comprehenr)
library(AICcmodavg)
library(rlist)
'%!in%' <- function(x,y)!('%in%'(x,y))

#Set source for stepaicc module and data path (need to modify to match repo location)
source("/path/to/repo/EpistasisStats/analysis/stepaicc/stepAICc.r")
setwd("/path/to/repo/EpistasisStats/data/processed/")
dep_var <- "ep_abs"

#Read in data, define features (non abstracted)
data<-read.table("skempi_bind_processed.txt", header=T, sep='\t')
features <- list("separation", "cplx_type", "intside_int", "size_net")

#Get list of unique systems in the dataset
unique_systems <- unique(data$cplx)

#Determine number to remove based on desired percentage
number_to_remove <- round(0.1 * length(unique_systems))


#Generate list of abstractions to build permutated lists
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


#Returns the AICc for a given bestfit model using model selection based on AICc
randomizedAIC <- function(l) {
  indata <- noquote(indata)
  randomized = stri_join_list(as.list(sample(stri_split(l, fixed=" + ")[[1]])), collapse=" + ")
  fullform = paste("ep_abs ~ ", paste(paste("(",randomized), ")^2"))
  rfit <- glm(fullform, data=data)  
  fitStep <- stepAICc(rfit, direction='both',steps=100000)
  return(AICc(fitStep))
}

#Returns the bestfit model by model selection using AICc as criteria
randomizedfit <- function(l) {
  randomized = stri_join_list(as.list(sample(stri_split(l, fixed=" + ")[[1]])), collapse=" + ")
  fullform = paste("ep_abs ~ ", paste(paste("(",randomized), ")^2"))
  rfit <- glm(fullform, data=data)  
  fitStep <- stepAICc(rfit, direction='both',steps=100000)
  return(fitStep)
}

#Same as above, just with cplx added as a variable
randomizedAIC_cplx <- function(l) {
  randomized = stri_join_list(as.list(sample(stri_split(l, fixed=" + ")[[1]])), collapse=" + ")
  fullform = paste("ep_abs ~ ", paste(paste("(",randomized), ")^2+cplx"))
  rfit <- glm(fullform, data=data)  
  fitStep <- stepAICc(rfit, direction='both',steps=100000)
  return(AICc(fitStep))
}

#Same as above, just with cplx added as a variable
randomizedfit_cplx <- function(l) {
  randomized = stri_join_list(as.list(sample(stri_split(l, fixed=" + ")[[1]])), collapse=" + ")
  fullform = paste("ep_abs ~ ", paste(paste("(",randomized), ")^2+cplx"))
  rfit <- glm(fullform, data=data)  
  fitStep <- stepAICc(rfit, direction='both',steps=100000)
  return(fitStep)
}

#Same as above, just with a filtered dataset (specified percentage of systems removed)
randomizedAIC_filtered <- function(l) {
  randomized = stri_join_list(as.list(sample(stri_split(l, fixed=" + ")[[1]])), collapse=" + ")
  fullform = paste("ep_abs ~ ", paste(paste("(",randomized), ")^2"))
  rfit <- glm(fullform, data=filtered_data)  
  fitStep <- stepAICc(rfit, direction='both',steps=100000)
  return(AICc(fitStep))
}

#Same as above, just with a filtered dataset (specified percentage of systems removed)
randomizedfit_filtered <- function(l) {
  randomized = stri_join_list(as.list(sample(stri_split(l, fixed=" + ")[[1]])), collapse=" + ")
  fullform = paste("ep_abs ~ ", paste(paste("(",randomized), ")^2"))
  rfit <- glm(fullform, data=filtered_data)  
  fitStep <- stepAICc(rfit, direction='both',steps=100000)
  return(fitStep)
}

#Gets the delta r-squared for all parameters in a given model
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


#Same as above, just for the dataset excluding the percentage of complexes removed
delta_r_sq_filtered <- function(fullmodel) {
  fmrsq <- summary(fullmodel)$r.squared
  drsq <- list(paste("Full model: ",fmrsq))
  terms <-  stri_split(c(bestmodel1$terms[[3]]), fixed=" + ")[[1]]
  singles <- to_list(for (elem in terms) if (grepl(":", elem) != TRUE) c(elem))
  for (elem in singles) {
    sublist <- stri_join_list(to_list(for (t_elem in terms) if (grepl(elem, t_elem) != TRUE) c(t_elem)), 
                              collapse=" + ")
    nrsq <- summary(lm(paste("ep_abs ~ ", sublist), data=filtered_data))$r.squared
    drsq <- list.append(drsq, paste(elem,abs(nrsq - fmrsq),sep=": "))
  }
  return(drsq)
}


#Generates 100 samples, saves data to file
i <- 99
for (iteration in 0:i) {
  removed_subset <- sample(unique_systems, number_to_remove)
  filtered_data <- data[data$cplx %!in% removed_subset,]
  run <- to_list(for (elem in allcombos) c(randomizedAIC_filtered(elem)))
  bestmodel1 <- lm(randomizedfit_filtered(allcombos[[which.min(run)]]))
  sink("../validation_results/leave_10per_out_bind.txt", append=TRUE)
  print("\newentry")
  print(i)
  print("removed complexes: ")
  print(removed_subset)
  print("R_sq ranking: ")
  print(delta_r_sq_filtered(bestmodel1))
  print(summary(lm(bestmodel1)))
  sink()
}