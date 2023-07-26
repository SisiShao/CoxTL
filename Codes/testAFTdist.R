rm(list=ls())

# Read data and true omegas
i = 9
y_name = paste("data/source/source_Y", as.character(i), ".csv", sep="")
Ys = read_csv(y_name)
colnames(Ys) = c("Time", "Event")
X_name = paste("data/source/source_X", as.character(i), ".csv", sep="")
Xs = read_csv(X_name)

data = cbind(Ys, Xs)

true_hi = read_csv("data/omega.csv")

# Fit the Lognormal model
model1 <- survreg(Surv(time = Time, event = Event) ~ . -1, 
                 data = data, dist = "lognormal")
model2 <- survreg(Surv(time = Time, event = Event) ~ . -1, 
                 data = data, dist = "weibull")
model3 <- survreg(Surv(time = Time, event = Event) ~ . -1, 
                  data = data, dist = "loglogistic")

idx = 30
plot(true_hi$Omega_9[1:idx], pch=19, type="b", col="black")
lines(model1$coefficients[1:idx], col="red", pch=19, type="b")
lines(model2$coefficients[1:idx], col="green", pch=19, type="b")
lines(model3$coefficients[1:idx], col="blue", pch=19, type="b")

norm(true_hi$Omega_9 - model1$coefficients, type="2")
norm(true_hi$Omega_9 - model2$coefficients, type="2")
norm(true_hi$Omega_9 - model3$coefficients, type="2")
