library("COUNT")
data("badhealth")
head(badhealth)
any(is.na(badhealth))
hist(badhealth$numvisit, breaks=20)

plot(jitter(log(numvisit)) ~ jitter(age), data=badhealth, subset=badh==0, xlab="age", ylab="log(visits)")
points(jitter(log(numvisit)) ~ jitter(age), data=badhealth, subset=badh==1, col="red")

library("rjags")
mod_string = " model {
    for (i in 1:length(numvisit)) {
        numvisit[i] ~ dpois(lam[i])
        log(lam[i]) = int + b_badh*badh[i] + b_age*age[i] + b_intx*age[i]*badh[i]
    }
    
    int ~ dnorm(0.0, 1.0/1e6)
    b_badh ~ dnorm(0.0, 1.0/1e4)
    b_age ~ dnorm(0.0, 1.0/1e4)
    b_intx ~ dnorm(0.0, 1.0/1e4)
} "

set.seed(102)

data_jags = as.list(badhealth)

params = c("int", "b_badh", "b_age", "b_intx")

mod = jags.model(textConnection(mod_string), data=data_jags, n.chains=3)
update(mod, 1e3)

mod_sim = coda.samples(model=mod,
                       variable.names=params,
                       n.iter=5e3)
mod_csim = as.mcmc(do.call(rbind, mod_sim))

## convergence diagnostics
plot(mod_sim)

gelman.diag(mod_sim)
autocorr.diag(mod_sim)
autocorr.plot(mod_sim)
effectiveSize(mod_sim)

## compute DIC
dic = dic.samples(mod, n.iter=1e3)

X = as.matrix(badhealth[,-1])
X = cbind(X, with(badhealth, badh*age))
head(X)

(pmed_coef = apply(mod_csim, 2, median))

llam_hat = pmed_coef["int"] + X %*% pmed_coef[c("b_badh", "b_age", "b_intx")]
lam_hat = exp(llam_hat)

hist(lam_hat)

resid = badhealth$numvisit - lam_hat
plot(resid) # the data were ordered
plot(lam_hat, badhealth$numvisit)
abline(0.0, 1.0)
plot(lam_hat[which(badhealth$badh==0)], resid[which(badhealth$badh==0)], xlim=c(0, 8), ylab="residuals", xlab=expression(hat(lambda)), ylim=range(resid))
points(lam_hat[which(badhealth$badh==1)], resid[which(badhealth$badh==1)], col="red")
var(resid[which(badhealth$badh==0)])

var(resid[which(badhealth$badh==1)])
summary(mod_sim)

x1 = c(0, 35, 0) # good health
x2 = c(1, 35, 35) # bad health

head(mod_csim)

loglam1 = mod_csim[,"int"] + mod_csim[,c(2,1,3)] %*% x1
loglam2 = mod_csim[,"int"] + mod_csim[,c(2,1,3)] %*% x2

lam1 = exp(loglam1)
lam2 = exp(loglam2)

(n_sim = length(lam1))

y1 = rpois(n=n_sim, lambda=lam1)
y2 = rpois(n=n_sim, lambda=lam2)

plot(table(factor(y1, levels=0:18))/n_sim, pch=2, ylab="posterior prob.", xlab="visits")
points(table(y2+0.1)/n_sim, col="red")

mean(y2 > y1)
dic_result1 = dic.samples(mod, n.iter=00000)
print(dic_result1)

############ alternative model ohne interactive terms #####################
mod_string_alt = " model {
    for (i in 1:length(numvisit)) {
        numvisit[i] ~ dpois(lam[i])
        log(lam[i]) = int + b_badh*badh[i] + b_age*age[i]
    }
    
    int ~ dnorm(0.0, 1.0/1e6)
    b_badh ~ dnorm(0.0, 1.0/1e4)
    b_age ~ dnorm(0.0, 1.0/1e4)
    b_intx ~ dnorm(0.0, 1.0/1e4)
} "

data_jags = as.list(badhealth)

params = c("int", "b_badh", "b_age", "b_intx")

mod1 = jags.model(textConnection(mod_string_alt), data=data_jags, n.chains=3)
update(mod1, 1e3)

mod1_sim = coda.samples(model=mod1,
                       variable.names=params,
                       n.iter=5e3)
dic_result2 = dic.samples(mod1, n.iter=100000)
print(dic_result2)




