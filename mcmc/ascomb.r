library("rjags")

mod_string = " model {
    for (i in 1:length(education)) {
        education[i] ~ dnorm(mu[i], prec)
        mu[i] = b0 + b[1]*income[i] + b[2]*young[i] + b[3]*urban[i]
    }
    
    b0 ~ dnorm(0.0, 1.0/1.0e6)
    for (i in 1:3) {
        b[i] ~ dnorm(0.0, 1.0/100)
    }
    
    prec ~ dgamma(1.0/2.0, 1.0*1500.0/2.0)
    	## Initial guess of variance based on overall
    	## variance of education variable. Uses low prior
    	## effective sample size. Technically, this is not
    	## a true 'prior', but it is not very informative.
    sig2 = 1.0 / prec
    sig = sqrt(sig2)
} "

data_jags = as.list(Anscombe)
set.seed(72)

params1 = c("b", "sig")

inits1 = function() {
  inits = list("b"=rnorm(3, 0.0, 100.0), "prec"=rgamma(1, 1.0, 1.0))
}

mod1 = jags.model(textConnection(mod_string), data=data_jags, inits=inits1, n.chains=3)
update(mod1, 1000) # burn-in

mod1_sim = coda.samples(model=mod1,
                        variable.names=params1,
                        n.iter=100000)

mod1_csim = do.call(rbind, mod1_sim) # combine multiple chains
#plot(mod1_csim)
plot(mod)
#autocorr.plot(mod1_sim)
#plot(resid(mod))
dic_result = dic.samples(mod1, n.iter=100000)
print(dic_result)
print(summary(mod1_csim))
print(mean(mod1_csim[, "b[1]"] > 0))
quantile(mod1_csim[, "b[1]"], probs = c(0.025, 0.975))