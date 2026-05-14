library(rjags)
data("warpbreaks")
?warpbreaks
head(warpbreaks)

mod1_string = " model {
    for(i in 1:length(y)) {
        # 1. Fehlende schließende Klammer bei 'prec' behoben
        y[i] ~ dnorm(mu[tensGrp[i]], prec[tensGrp[i]])
    }
    
    # 2. Schleife angepasst: warpbreaks$tension hat 3 Stufen, 
    # aber Sie möchten separate Varianzen/Präzisionen für jede Gruppe!
    for (j in 1:3) {
        mu[j] ~ dnorm(0.0, 1.0/1.0e6)
        
        # 3. Inverse-Gamma(1/2, 1/2) Prior für die Varianz (gemäß Aufgabenstellung)
        # In JAGS wird dgamma(shape, rate) für die Präzision genutzt.
        # Ein Inv-Gamma(a, b) für die Varianz entspricht einem Gamma(a, b) für die Präzision.
        prec[j] ~ dgamma(0.5, 0.5)
        sig[j] <- sqrt(1.0 / prec[j]) # Pfeil-Zuweisung für deterministische Knoten
    }
} "

set.seed(83)

# Datenvorbereitung
data1_jags = list(y = log(warpbreaks$breaks), tensGrp = as.numeric(warpbreaks$tension))

params1 = c("mu", "sig")

# Modell initialisieren
mod1 = jags.model(textConnection(mod1_string), data=data1_jags, n.chains=3)
update(mod1, 1e3)

# MCMC-Stichproben ziehen
mod1_sim = coda.samples(model=mod1,
                        variable.names=params1,
                        n.iter=5e3)

## Konvergenzdiagnostik
plot(mod1_sim)
gelman.diag(mod1_sim)
autocorr.diag(mod1_sim)
effectiveSize(mod1_sim)
dic_mod1 <- dic.samples(model = mod1, 
                        n.iter = 5000, 
                        type = "pD")

# Ergebnis anzeigen
print(dic_mod1)
