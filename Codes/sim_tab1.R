rm(list=ls())

CoxKL_result <- function(file_source, file_target, file_true_source, 
                            file_true_target, num_covariates, max_eta=1, gap=0.05) {
  # Loading data
  df_source = read.csv(file_source)
  df_target = read.csv(file_target)
  true_target = read.csv(file_true_target)
  true_source = read.csv(file_true_source)
  
  est_source = coxph(Surv(y_s_obs, delta_s) ~ ., data=df_source)
  RS_external = as.matrix(df_target[1:nrow(df_target)/2, 1:num_covariates]) %*%
    est_source$coefficients[1:num_covariates]
  
  eta_KL <- function(eta) {
    est = KL_Cox_Estimate(z = df_target[1:nrow(df_target)/2,1:num_covariates], 
                          delta = df_target$delta_t[1:nrow(df_target)/2],
                          time=df_target$y_t_obs[1:nrow(df_target)/2],
                          RS_internal=RS_external, 
                          eta=eta, tol=1.0e-7)
    
    mse_val = norm(est - true_target$X0[1:num_covariates], "2")
    cindex_val = glmnet::Cindex(as.matrix(df_target[-c(1:nrow(df_target)/2),
                                                    1:num_covariates])%*%est[1:num_covariates], 
                                Surv(df_target$y_t_obs[-c(1:nrow(df_target)/2)], 
                                     df_target$delta_t[-c(1:nrow(df_target)/2)]))
    
    return (data.frame(MSE=mse_val,Cindex=cindex_val)) 
    }
    
    eta_list = seq(0, max_eta, gap)
    results = mclapply(eta_list, eta_KL)
    mse_values <- sapply(results, function(df) df$MSE)
    cindex_values = sapply(results, function(df) df$Cindex)
    mse_eta_optimal = eta_list[which.min(mse_values)]
    
    return (data.frame(eta=eta_list, MSE=mse_values, Cindex=cindex_values))
}

get_optimal_val <- function(df) {
  "df consists of 3 columns: eta, MSE, and Cindex."
  eta_list = df$eta
  mse_values = df$MSE
  cindex_val = df$Cindex
  mse_eta_optimal = eta_list[which.min(mse_values)]
  mse_min = min(mse_values)
  cindex_eta_optimal = eta_list[which.max(cindex_val)]
  cindex_max = max(cindex_val)
  
  return(data.frame(mse_min=mse_min, optimal_mse_eta=mse_eta_optimal,
                    cindex_max = cindex_max, optimal_cindex_eta=cindex_eta_optimal))
}

all_files = list.files("./")

all_results <- list()
hi_results <- list()

for (i in 1:100) {
  print(i)
  if (i == 55) next
  hi = CoxKL_result(all_files[i], 
                    all_files[i+100],
                    all_files[i+200],
                    all_files[i+300], num_covariates = 15, max_eta = 3, gap=0.05)
  
  result = get_optimal_val(hi)
  all_results[[i]] <- result
  hi_results[[i]] <- hi
}

final_df <- do.call(rbind, all_results)
# Create a sequence of numbers from 1 to 99, but remove 55
numbers <- setdiff(1:100, 55)
# Rep each number 61 times
desired_vector <- rep(numbers, each=3/0.05 + 1)
hi_df <- do.call(rbind, hi_results) %>%
  add_column(sim_idx=desired_vector)
head(final_df)
head(hi_df)

#write_csv(final_df, file = "../sim_tab1_result/optimal_vals.csv")
#write_csv(hi_df, file = "../sim_tab1_result/hi_df.csv")

library(viridis)

ggplot(hi_df) +
  geom_line(aes(x=eta, y=MSE, group=sim_idx,  col=as.factor(sim_idx)),
            cex=.5, show.legend = FALSE) +
  theme_bw() #+
  #scale_color_brewer(palette="RdPu") +
  #ylim(0.82, 0.97)

library(gtsummary)
library(dplyr)

table_summary <- final_df %>% 
  tbl_summary(
    statistic = list(all_continuous() ~ "{mean} ({sd})"),
    ) 
table_summary

hi_df %>%
  tbl_summary(
    statistic = list(all_continuous() ~ "{mean} ({sd})")
  )
