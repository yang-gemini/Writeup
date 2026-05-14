library(coda)
library("MASS")
library(rjags)
data("OME")
?OME # background on the data
head(OME)

any(is.na(OME)) # check for missing values
dat = subset(OME, OME != "N/A") # manually remove OME missing values identified with "N/A"
dat$OME = factor(dat$OME)
str(dat)

plot(dat$Age, dat$Correct / dat$Trials )
plot(dat$OME, dat$Correct / dat$Trials )
plot(dat$Loud, dat$Correct / dat$Trials )
plot(dat$Noise, dat$Correct / dat$Trials )
mod_glm = glm(Correct/Trials ~ Age + OME + Loud + Noise, data=dat, weights=Trials, family="binomial")
summary(mod_glm)
plot(residuals(mod_glm, type="deviance"))
plot(fitted(mod_glm), dat$Correct/dat$Trials)
X = model.matrix(mod_glm)[,-1] # -1 removes the column of 1s for the intercept
head(X)
mod_string = " model {
	for (i in 1:length(y)) {
		y[i] ~ dbin(phi[i], n[i])
		logit(phi[i]) = b0 + b[1]*Age[i] + b[2]*OMElow[i] + b[3]*Loud[i] + b[4]*Noiseincoherent[i]
	}
	
	b0 ~ dnorm(0.0, 1.0/5.0^2)
	for (j in 1:4) {
		b[j] ~ dnorm(0.0, 1.0/4.0^2)
	}
	
} "

data_jags = as.list(as.data.frame(X))
data_jags$y = dat$Correct # this will not work if there are missing values in dat (because they would be ignored by model.matrix). Always make sure that the data are accurately pre-processed for JAGS.
data_jags$n = dat$Trials
str(data_jags) # make sure that all variables have the same number of observations (712).
mod = jags.model(textConnection(mod_string), data = data_jags, n.chains = 1)

# 3. Einschwingphase (Burn-In): Erste 1.000 Iterationen verwerfen
update(mod, 1000)

# 4. Eine kurze Testkette von z. B. 2.000 Iterationen für die Diagnose ziehen
samples = coda.samples(model = mod, 
                       variable.names = c("b0", "b"), 
                       n.iter = 2000)

mcmc_kette = as.mcmc(samples)
raftery.diag(mcmc_kette, q = 0.025, r = 0.005, s = 0.95)
mod_summary <- summary(samples)
print(mod_summary)

# 1. Posteriore Mittelwerte (Log-Odds) extrahieren
koeffizienten <- summary(samples)$statistics[, "Mean"]

# Falls die Reihenfolge im Vektor abweicht, weisen wir sie explizit zu:
b0 <- koeffizienten["b0"]
b1 <- koeffizienten["b[1]"] # Age
b2 <- koeffizienten["b[2]"] # OMElow
b3 <- koeffizienten["b[3]"] # Loud
b4 <- koeffizienten["b[4]"] # Noiseincoherent

# 2. Linearen Prädiktor (eta) für das spezifische Kind berechnen
# OMElow = 0 und Noiseincoherent = 0, da das Kind hohes OME & kohärenten Stimulus hat
eta <- b0 + (b1 * 60) + (b2 * 0) + (b3 * 50) + (b4 * 0)

# 3. Inverse-Logit-Transformation zur Berechnung der Wahrscheinlichkeit (phi)
wahrscheinlichkeit <- 1 / (1 + exp(-eta))

# Ergebnis anzeigen
print(wahrscheinlichkeit)
# Koeffizienten sauber trennen (Intercept vs. Steigungskoeffizienten)
b0_est <- koeffizienten["b0"]
b_est  <- koeffizienten[c("b[1]", "b[2]", "b[3]", "b[4]")]

# Schritt 2: Linearen Prädiktor (eta) für alle Beobachtungen berechnen
# %*% ist die Matrix-Multiplikation zwischen der Designmatrix X und dem Koeffizientenvektor
eta_all <- b0_est + as.vector(X %*% b_est)

# Schritt 3: Inverse-Logit-Transformation (Zuweisung an das geforderte Objekt 'phat')
phat <- 1 / (1 + exp(-eta_all))

# Überprüfung des Vektors
str(phat)
head(phat)
(tab0.7 = table(phat > 0.7, (dat$Correct / dat$Trials) > 0.7))
sum(diag(tab0.7)) / sum(tab0.7)