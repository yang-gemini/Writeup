update_mu = function(n, ybar, sig2, mu_0, sig2_0) {
  sig2_1 = 1.0 / (n / sig2 + 1.0 / sig2_0)
  mu_1 = sig2_1 * (n * ybar / sig2 + mu_0 / sig2_0)
  rnorm(n=1, mean=mu_1, sd=sqrt(sig2_1))
}

gibbs = function(d, n_coeff, n_iter, init, prior) {
  
  mu_out <- matrix(0, nrow = n_coeff, ncol = n_iter)

  mu_now = init$mu
  
  ## Gibbs sampler
  for (i in 1:n_iter) {
    for (j in 1:n_coeff) {
        mu_now[j] = update_mu(n=n, ybar=ybar, sig2=sig2_now, mu_0=prior$mu_0, sig2_0=prior$sig2_0)
    }
    mu_out[i] = mu_now
  }
}

# y = c(1.2, 1.4, -0.5, 0.3, 0.9, 2.3, 1.0, 0.1, 1.3, 1.9)
y = c(-0.2, -1.5, -5.3, 0.3, -0.8, -2.2)
ybar = mean(y)
n = length(y)

## prior
prior = list()
prior$mu_0 = 1.0
prior$sig2_0 = 1.0
prior$n_0 = 6.0 # prior effective sample size for sig2
prior$s2_0 = sd(y) # prior point estimate for sig2
prior$nu_0 = prior$n_0 / 2.0 # prior parameter for inverse-gamma
prior$beta_0 = prior$n_0 * prior$s2_0 / 2.0 # prior parameter for inverse-gamma

hist(y, freq=FALSE, xlim=c(-1.0, 3.0)) # histogram of the data
curve(dnorm(x=x, mean=prior$mu_0, sd=sqrt(prior$sig2_0)), lty=2, add=TRUE) # prior for mu
points(y, rep(0,n), pch=1) # individual data points
points(ybar, 0, pch=19) # sample mean

set.seed(53)

init = list()
init$mu = 0.0

post = gibbs(y=y, n_iter=5000, init=init, prior=prior)
print(mean(post[,1]))





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