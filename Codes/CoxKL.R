setwd("/Users/Eliuvish/Downloads/SSS/Codes/CoxKL")
library(tidyverse)
library(knitr)
library(ciTools)
library(survival)

source_param_est <- matrix(NA, nrow=300, ncol=10)

for (i in 1:10) {
  print(paste("冲！\n Process...... file no.", i))
  
  y_name = paste("data/source/source_Y", as.character(i), ".csv", sep="")
  Ys = read_csv(y_name)
  colnames(Ys) = c("Time", "Event")
  
  X_name = paste("data/source/source_X", as.character(i), ".csv", sep="")
  Xs = read_csv(X_name)
  
  data = cbind(Ys, Xs)
  
  # Fit the Lognormal model
  model <- survreg(Surv(time = Time, event = Event) ~ . -1, 
                   data = data, dist = "lognormal")
  
  # Print a summary of the model
  # summary(model)
  
  source_param_est[, i] <- model$coefficients
}

write_csv(as.data.frame(source_param_est), "results/source_param_est.csv")
