#install.packages("devtools")

rm(list=ls())
setwd("/Users/sisishao/CoxKL")
library(devtools)
#install_github("SisiShao/CoxKL")
library(CoxKL)
library(survival)

df_source= read.csv("final_data_s_rho0.6_CS0.78_S1000_T100.csv")
df_target= read.csv("final_data_t_rho0.6_CS0.78_S1000_T100.csv")
true_target = read.csv("TrueT_rho0.6_CS0.78_S1000_T100")
true_source = read.csv("TrueS_rho0.6_CS0.78_S1000_T100")

est_source = coxph(Surv(y_s_obs, delta_s) ~ ., data=df_source)
est_source$coefficients
true_source$X0

nrows = dim(df_target)[1]/2

num_covariates = 5
RS_external = as.matrix(df_target[1:nrows, 1:num_covariates]) %*% est_source$coefficients[1:num_covariates] # calculated by external parameters and internal covariates
#RS_external = as.matrix(df_target[1:nrows, 1:num_covariates]) %*% true_source$X0[1:num_covariates] # calculated by external parameters and internal covariates

#Use CoxKL estimate
eta_KL <- function(eta) {
    est = KL_Cox_Estimate(z = df_target[1:nrows,1:num_covariates], 
                          delta = df_target$delta_t[1:nrows],
                    time=df_target$y_t_obs[1:nrows],RS_internal=RS_external, 
                    eta=eta, tol=1.0e-7)
    
    mse = norm(est - true_target$X0[1:num_covariates], "2")
    cindex = glmnet::Cindex(as.matrix(df_target[-c(1:nrows),1:num_covariates])%*%est[1:num_covariates], 
                            Surv(df_target$y_t_obs[-c(1:nrows)], 
                                 df_target$delta_t[-c(1:nrows)]))
    
    cindex
    #mse
    #est
}

library(parallel)

eta_list = seq(0, 10, 0.05)
C_index_good_quality = mclapply(eta_list, eta_KL)
eta_max = eta_list[which.max(C_index_good_quality)]
C_max = C_index_good_quality[which.max(C_index_good_quality)]

plot(eta_list, C_index_good_quality, type="l", col="darkgreen", lwd=2,
     ylab = "C-index", xlab="eta", ylim = c(min(unlist(C_index_good_quality)-0.002),
                                            max(unlist(C_index_good_quality)+0.001)))

points(eta_max, C_max, col="red", cex=2, pch=19)
text(eta_max, unlist(C_max), 
     labels = paste("eta: ", round(eta_max, 2), "\nMSE: ", round(unlist(C_max), 2)),
     pos=1, cex=0.8, offset=.5)


