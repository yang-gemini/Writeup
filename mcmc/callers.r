# 1. Daten einlesen (wir nutzen ein sauberes Data Frame statt cbind)
dat = read.csv(file="callers.csv", header=TRUE)
head(dat)
pairs(dat)

# Wir behalten die originalen Integer-Calls für dpois und nutzen Spalte 2 als Offset
new_dat <- data.frame(
  calls   = dat[, 1],
  days_active = dat[, 2],
  isgroup2   = dat[, 3],
  age        = dat[, 4]
)
pairs(new_dat)

# 2. JAGS Modell mit Offset-Formel definieren
library("rjags")

mod_string = " model {
    for (i in 1:length(calls)) {
      	calls[i] ~ dpois( days_active[i] * lam[i] )
		    log(lam[i]) = b0 + b[1]*age[i] + b[2]*isgroup2[i]
    } 
    
    # Uninformative Priors (Präzision = 0.01)
    b0 ~ dnorm(0.0, 0.01) 
    for(j in 1:2) {
        b[j] ~ dnorm(0.0, 0.01)
    }
} "

set.seed(102)

# Daten für JAGS als Liste vorbereiten
data_jags = as.list(new_dat)

params = c("b0", "b")

# Modell kompiliere und initialisieren
mod = jags.model(textConnection(mod_string), data=data_jags, n.chains=3)

# 3. Burn-in & eigentliche MCMC-Simulation durchführen (Schließt die Lücke zu mod_sim!)
update(mod, n.iter=1000) # 1000 Iterationen zum "Einlaufen" verwerfen
mod_sim = coda.samples(model=mod, variable.names=params, n.iter=5000)

## 4. Konvergenzdiagnostik (Funktioniert jetzt, da mod_sim existiert)
plot(mod_sim)

gelman.diag(mod_sim)
autocorr.diag(mod_sim)
autocorr.plot(mod_sim)
effectiveSize(mod_sim)

## 5. DIC berechnen
# Wichtig: dic.samples benötigt das ursprüngliche 'mod'-Objekt, nicht die Samples!
dic = dic.samples(mod, n.iter=1000)
print(dic)
mod_csim <- as.matrix(mod_sim)
mean(mod_csim[, "b[2]"] > 0)

## Vorhersagen (Alter = 29, Gruppe = 2, Tage aktiv = 30)
v <- c(29, 1)
new_x <- matrix(v, nrow = 2)

# Koeffizienten-Mittelwerte extrahieren
p_coef <- apply(mod_csim, 2, mean)

# 1. Zeige die extrahierten Koeffizienten zur Kontrolle
print(p_coef[c("b[1]", "b[2]")])
print(p_coef["b0"])

# 2. Berechne den linearen Prädiktor auf der Log-Skala
# Das Skalarprodukt %*% muss komplett in die Klammer von exp()
log_lam <- p_coef["b0"] + (p_coef[c("b[1]", "b[2]")] %*% new_x)

# 3. Transformieren mit exp() und mit dem Offset (30 Tage) multiplizieren
lam_hat <- 30 * exp(log_lam)
print(lam_hat)

# 4. Wahrscheinlichkeit für höchstens 3 Calls berechnen
res <- 1 - ppois(q = 2, lambda = lam_hat)
print(res)
