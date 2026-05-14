dat <- read.table("pgalpga2008.txt", header=TRUE,  col.names = c("Distance", "Precision", "FM"))
datF <- subset(dat, FM==1, select=1:2)
plot(datF$Distance, datF$Precision, 
     main = "PGA 2008 Plot", 
     xlab = "Distance", 
     ylab = "Precision", 
     pch = 19,      # Ausgefüllte Punkte
     col = "blue")  # Farbe der Punkte
modell <- lm(Precision ~ Distance, data = datF)
sy = summary(modell)
print(sy)
neue_daten <- data.frame(Distance = 260)
y_hat <- predict(modell, newdata = neue_daten, se.fit = TRUE)
print(y_hat$fit)
sd_hat <- y_hat$se.fit
print(sd_hat)
#vorhersagen konfidenz interval
n <- nrow(datF)
df_modell <- modell$df.residual
se_r <- sy$sigma  # Das ist se_r aus deiner Formel
x_quer <- mean(datF$Distance)
s2_x <- var(datF$Distance)
se_pred <- se_r * sqrt(1 + 1/n + (260 - x_quer)^2 / ((n - 1) * s2_x))
t_wert <- qt(0.975, df = df_modell)
cf_manuell <- y_hat$fit + c(-1, 1) * t_wert * se_pred
print(cf_manuell)
