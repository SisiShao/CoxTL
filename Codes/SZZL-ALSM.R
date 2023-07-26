rm(list=ls())
library(glmnet)
library(penalized)

# Read Data
true_beta = read_csv("data/beta.csv")
X = read_csv("data/train_test/train_X.csv")
Y = read_csv("data/train_test/train_y.csv"); colnames(Y) = c("Time", "Event")
Omega = read_csv("data/omega.csv")
RS = as.matrix(X) %*% as.matrix(Omega)
df_interval = cbind(X, status=Y$Event, time=Y$Time)

# Run Cox and Cox-KL
mod_KL = KL_Cox_Estimate(X[, -1], delta = Y$Event, time = Y$Time, RS_internal = RS, eta = 0.5)
mod_Cox = coxph(Surv(time, status) ~ . -1, data=df_interval)
mod_KL
mod_Cox$coefficients

# Plot Results
idx=30
plot(mod_KL[1:idx], pch=19, type="b", col="red")
lines(true_beta$Beta[1:idx], pch=19, type="b", col="black")

norm(mod_KL - true_beta$Beta[-1], type="2")
norm(mod_Cox$coefficients[-1] - true_beta$Beta[-1], type="2")

# Shao-Zhou-Zhou-Li's Aggregated Least Square Method (SZZL-ALSM)
X_test = read_csv("data/train_test/test_X.csv")
Y_test = read_csv("data/train_test/test_y.csv"); colnames(Y_test) = c("Time", "Event")
Gang_regressor = cbind(as.matrix(X_test) %*% as.matrix(Omega),
                       as.matrix(X_test) %*% c(0, mod_KL))
colnames(Gang_regressor)[11] = "KL"

# Cox's Multiplicative Intensity Model
X_Gang = Gang_regressor
Y_Gang = Surv(Y_test$Time, Y_test$Event)
cv_fit_Gang = cv.glmnet(X_Gang, Y_Gang, alpha=.0, family="cox")

fit_Gang = glmnet(X_Gang, Y_Gang, alpha=.0, family="cox", lambda = cv_fit_Gang$lambda.min)
fit_Gang$beta

fit_Cox = glmnet(X_Gang, Y_Gang, alpha=1, family="cox", lambda = 0)
fit_Cox$beta

# Penalized Accelerated Failure Time
aft_model <- penalized(Y_Gang ~ ., 
                       data = as.data.frame(X_Gang), 
                       #penalized = ~ ., # specify which variables to penalize
                       lambda1 = 0, # Lasso penalty
                       lambda2 = 1, # Ridge penalty
                       model = "cox")

# Print the model
print(aft_model)
aft_model@fitted
