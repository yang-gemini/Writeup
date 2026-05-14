data("PlantGrowth")
?PlantGrowth
head(PlantGrowth)

mod_string = " model {
    for (i in 1:length(y)) {
        # Präzision hängt jetzt von der Gruppe des Datenpunkts ab
        y[i] ~ dnorm(mu[grp[i]], prec[grp[i]])
    }

    for (j in 1:3) {
        mu[j] ~ dnorm(0.0, 1.0/1.0e6)
        # Jede Gruppe bekommt eine eigene Prior-Verteilung für die Präzision
        prec[j] ~ dgamma(5/2.0, 5*1.0/2.0)
        # Berechnung der separaten Standardabweichungen
        sig[j] = sqrt( 1.0 / prec[j] )
    }
} "

set.seed(82)
str(PlantGrowth)
data_jags = list(y=PlantGrowth$weight, 
                 grp=as.numeric(PlantGrowth$group))

params = c("mu", "sig")

inits = function() {
  inits = list("mu"=rnorm(3,0.0,100.0), "prec"=rgamma(3,1.0,1.0))
}

mod = jags.model(textConnection(mod_string), data=data_jags, inits=inits, n.chains=3)
update(mod, 1e3)

mod_sim = coda.samples(model=mod,
                       variable.names=params,
                       n.iter=5e3)
mod_csim = as.mcmc(do.call(rbind, mod_sim)) # combined chains
par(mfrow = c(2, 2)) # Teilt das Fenster in ein 2x2 Raster
pdf("JAGS_Plots.pdf")
plot(mod_sim)
dev.off() # Schlie
gelman.diag(mod_sim)
autocorr.diag(mod_sim)
effectiveSize(mod_sim)
# DIC Berechnung (muss auf 'mod' zugreifen, nicht auf 'mod_csim')
dic_result2 = dic.samples(mod, n.iter=5e3)
print(dic_result2)

print(summary(mod_csim))