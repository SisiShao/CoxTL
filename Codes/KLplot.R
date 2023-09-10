plot_eta_cindex <- function(file_source, file_target, file_true_source, file_true_target, num_covariates, mse=FALSE, cindex=TRUE, max_eta=1) {
  # Loading data
  df_source = read.csv(file_source)
  df_target = read.csv(file_target)
  true_target = read.csv(file_true_target)
  true_source = read.csv(file_true_source)
  
  est_source = coxph(Surv(y_s_obs, delta_s) ~ ., data=df_source)
  RS_external = as.matrix(df_target[1:nrow(df_target)/2, 1:num_covariates]) %*% est_source$coefficients[1:num_covariates]
  
  eta_KL <- function(eta) {
    est = KL_Cox_Estimate(z = df_target[1:nrow(df_target)/2,1:num_covariates], 
                          delta = df_target$delta_t[1:nrow(df_target)/2],
                          time=df_target$y_t_obs[1:nrow(df_target)/2],RS_internal=RS_external, 
                          eta=eta, tol=1.0e-7)
    
    mse_val = norm(est - true_target$X0[1:num_covariates], "2")
    cindex_val = glmnet::Cindex(as.matrix(df_target[-c(1:nrow(df_target)/2),1:num_covariates])%*%est[1:num_covariates], 
                                Surv(df_target$y_t_obs[-c(1:nrow(df_target)/2)], 
                                     df_target$delta_t[-c(1:nrow(df_target)/2)]))
    
    if(mse) {
      return(mse_val)
    } else if(cindex) {
      return(cindex_val)
    } else {
      return(list(mse = mse_val, cindex = cindex_val))
    }
  }
  
  eta_list = seq(0, max_eta, 0.05)
  values = mclapply(eta_list, eta_KL)
  
  if(mse) {
    ylab_text <- "MSE"
    eta_optimal = eta_list[which.min(values)]
    optimal_value = values[which.min(values)]
  } else {
    ylab_text <- "C-index"
    eta_optimal = eta_list[which.max(values)]
    optimal_value = values[which.max(values)]
  }
  
  plot(eta_list, values, type="l", col="darkgreen", lwd=2,
       ylab = ylab_text, xlab="eta", ylim = c(min(unlist(values)-0.001),
                                              max(unlist(values)+0.008)))
  points(eta_optimal, optimal_value, col="red", cex=2, pch=19)
  text(eta_optimal, unlist(optimal_value), 
       labels = paste("eta: ", round(eta_optimal, 2), "\nValue: ", round(unlist(optimal_value), 2)),
       pos=4, cex=0.8, offset=.5)
}

# Example call for MSE
plot_eta_cindex("final_data_s_rho0.6_CS0.27_S2000_T500.csv", 
                "final_data_t_rho0.6_CS0.27_S2000_T500.csv",
                "TrueS_rho0.6_CS0.27_S2000_T500",
                "TrueT_rho0.6_CS0.27_S2000_T500",
                15, mse=TRUE, cindex=FALSE, max_eta = 1)

# Example call for C-index
plot_eta_cindex("final_data_s_rho0.6_CS0.27_S2000_T500.csv", 
                "final_data_t_rho0.6_CS0.27_S2000_T500.csv",
                "TrueS_rho0.6_CS0.27_S2000_T500",
                "TrueT_rho0.6_CS0.27_S2000_T500",
                15, mse=FALSE, cindex=TRUE, max_eta =1)