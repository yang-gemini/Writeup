yt = scan("eeg.txt")
plot.ts(yt, 
        main = "Zeitreihen-Verlauf", 
        col = "blue", 
        ylab = "Messwert", 
        xlab = "Zeit (t)")

set.seed(2021)
T=length(yt) # number of time points

## Case 1: Conditional likelihood
p=8
y=rev(yt[(p+1):T]) # response
X=t(matrix(yt[rev(rep((1:p),T-p)+rep((0:(T-p-1)),rep(p,T-p)))],p,T-p));
XtX=t(X)%*%X
XtX_inv=solve(XtX)
phi_MLE=XtX_inv%*%t(X)%*%y # MLE for phi
s2=sum((y - X%*%phi_MLE)^2)/(length(y) - p) #unbiased estimate for v

cat("\n MLE of conditional likelihood for phi: ", phi_MLE, "\n",
    "Estimate for v: ", s2, "\n")

#####################################################################################
##  AR(2) case 
### Posterior inference, conditional likelihood + reference prior via 
### direct sampling                 
#####################################################################################
library(MASS)
n_sample=500 # posterior sample size

## step 1: sample v from inverse gamma distribution
v_sample=1/rgamma(n_sample, (T-2*p)/2, sum((y-X%*%phi_MLE)^2)/2)

## step 2: sample phi conditional on v from normal distribution
phi_sample=matrix(0, nrow = n_sample, ncol = p)
for(i in 1:n_sample){
  phi_sample[i, ]=mvrnorm(1,phi_MLE,Sigma=v_sample[i]*XtX_inv)
}

## plot histogram of posterior samples of phi and nu
par(mfrow = c(3, 3), cex.lab = 1.3)
for(i in 1:8){
  hist(phi_sample[, i], xlab = bquote(phi), 
       main = bquote("Histogram of "~phi[.(i)]))
}

hist(v_sample, xlab = bquote(nu), main = bquote("Histogram of "~v))

### Frequency und Preiodic
poly_coefs <- c(1, -phi_MLE[1], -phi_MLE[2], -phi_MLE[3], -phi_MLE[4], -phi_MLE[5], -phi_MLE[6], -phi_MLE[7], -phi_MLE[8] )

# 2. Find the roots of the polynomial
roots <- polyroot(poly_coefs)

# 3. Calculate the reciprocal roots (1 / z)
reciprocal_roots <- 1 / roots

# 4. Calculate the Moduli (should be < 1 for a stationary process)
moduli <- Mod(reciprocal_roots)

# 5. Calculate the Periods (lambda) for complex roots
# Using the argument (angle) of the complex number: period = 2 * pi / |argument|
angles <- Arg(reciprocal_roots)
periods <- (2 * pi) / abs(angles)

# --- Output the results ---
cat("Reziproke Wurzeln:\n")
print(reciprocal_roots)

cat("\nModuli (Beträge):\n")
print(moduli)

cat("\nPerioden (in Zeitschritten):\n")
print(periods)