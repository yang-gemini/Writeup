library("rjags")

dat = read.csv(file="pctgrowth.csv", header=TRUE)
head(dat)

mod_string = " model {
for (i in 1:length(y)) {
  y[i] ~ dnorm(mu[grp[i]], prec1)
}

for (j in 1:5) {
  mu[j] ~ dnorm(theta,  prec2)
}
theta ~ dnorm(0, 1/1e6)
prec1 ~ dgamma(1.0/2.0, 3.0/2.0)
prec2 ~ dgamma(2.0/2.0, 2.0/2.0) 
sig2 <- 1.0/prec1
another_sig2 <- 1.0 / prec2
} "

set.seed(113)

# 1. Sicherstellen, dass die Gruppe als 1, 2, 3, 4, 5 formatiert ist
dat$grp <- as.integer(as.factor(dat$grp))

# 2. Daten explizit für JAGS als benannte Liste vorbereiten
data_jags = list(
  y   = dat$y,   # falls die Spalte in der CSV 'y' heißt
  grp = dat$grp  # falls die Spalte in der CSV 'grp' heißt
)

params = c("mu", "sig2")

mod = jags.model(textConnection(mod_string), data=data_jags, n.chains=3)
update(mod, 1e3)

mod_sim = coda.samples(model=mod,
                       variable.names=params,
                       n.iter=10e3)
mod_csim = as.mcmc(do.call(rbind, mod_sim))

## convergence diagnostics
plot(mod_sim)

gelman.diag(mod_sim)
autocorr.diag(mod_sim)
autocorr.plot(mod_sim)
effectiveSize(mod_sim)
print(summary(mod_csim))
## compute DIC
dic = dic.samples(mod, n.iter=1e3)

################## ANOVA +++++++++++++++++++
anova_mod_string = " model {
    for (i in 1:length(y)) {
        y[i] ~ dnorm(mu[grp[i]], prec)
    }
    
    for (j in 1:5) {
        mu[j] ~ dnorm(0.0, 1.0/1.0e6)
    }
    
    prec ~ dgamma(1.0/2.0, 3.0/2.0)
    sig = sqrt( 1.0 / prec )
} "

params = c("mu", "sig")

anova_inits = function() {
  inits = list("mu"=rnorm(5,0.0,100.0), "prec"=rgamma(1,1.0,1.0))
}

anova_mod = jags.model(textConnection(anova_mod_string), data=data_jags, inits=anova_inits, n.chains=3)
update(mod, 1e3)

anova_mod_sim = coda.samples(model=anova_mod,
                       variable.names=params,
                       n.iter=10e3)
anova_mod_csim = as.mcmc(do.call(rbind, mod_sim))

## convergence diagnostics
plot(anova_mod_sim)

gelman.diag(anova_mod_sim)
autocorr.diag(anova_mod_sim)
autocorr.plot(anova_mod_sim)
effectiveSize(anova_mod_sim)
print(summary(anova_mod_csim))
## compute DIC
anova_dic = dic.samples(anova_mod, n.iter=1e3)
print(anova_dic)
