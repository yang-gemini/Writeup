library(mvtnorm)

dat <- read.delim("earthquakes.txt")
y.sample = as.numeric(dat$Quakes)

## set up
p = 3                  ## order of AR process
n.all = length(y.sample) ## T, total number of data
n = n.all - p          ## Anzahl der effektiven Beobachtungen in Y

Y = matrix(y.sample[4:n.all], ncol=1)
Fmtx = matrix(c(y.sample[3:(n.all-1)], y.sample[2:(n.all-2)], y.sample[1:(n.all-3)]), nrow=p, byrow=TRUE)

## posterior inference
## set the prior
m0 = matrix(rep(0, p), ncol=1)
C0 = 10 * diag(p)
n0 = 0.02
d0 = 0.02

## calculate parameters that will be reused in the loop
e = Y - t(Fmtx) %*% m0
Q = t(Fmtx) %*% C0 %*% Fmtx + diag(n) # n war vorher nicht definiert!
Q.inv = chol2inv(chol(Q))
A = C0 %*% Fmtx %*% Q.inv
m = m0 + A %*% e
C = C0 - A %*% Q %*% t(A)
n.star = n + n0
d.star = as.numeric(t(e) %*% Q.inv %*% e) + d0 # Fehler behoben

n.sample = 5000
nu.sample = rep(0, n.sample)
phi.sample = matrix(0, nrow=n.sample, ncol=p)

# set.seed() AUSSERHALB der Schleife setzen!
for (i in 1:n.sample) {
  set.seed(123)
  nu.new = 1 / rgamma(1, shape=n.star/2, rate=d.star/2)
  nu.sample[i] = nu.new
  phi.new = rmvnorm(1, mean=m, sigma=nu.new * C)
  phi.sample[i, ] = phi.new
}

phi_mean = colMeans(phi.sample) 

# Wenn Sie statt eines einzelnen Punktes (y_hat) echte Konfidenzintervalle für die 
# Vorhersage wollen, müssen Sie y_hat für JEDE gezogene phi-Stichprobe berechnen:
y_hat_matrix = matrix(0, nrow=n.sample, ncol=n.all)
for(s in 1:n.sample) {
  y_hat_matrix[s, 1:3] = y.sample[1:3]
  for(i in 4:n.all) {
    y_hat_matrix[s, i] = y.sample[i-1] * phi.sample[s, 1] + y.sample[i-2] * phi.sample[s, 2] + y.sample[i-3] * phi.sample[s, 3]
  }
}

summary.vec95 = function(vec){
  c(unname(quantile(vec, 0.025)), mean(vec), unname(quantile(vec, 0.975)))
}

summary.y = apply(y_hat_matrix, MARGIN=2, summary.vec95)

# Plotten der Ergebnisse
png(filename = "./images/conjugate-prediction.png", 
    width = 1920,           # Breite in Pixeln
    height = 1080,           # Höhe in Pixeln
    res = 150)              # Auflösung erhöhen für scharfen Text (DPI)

plot(y.sample, type='b', xlab='Time', ylab='Quakes', pch=16, ylim=c(min(y.sample)-5, max(y.sample)+5))
lines(summary.y[2,], type='b', col='blue', lty=2, pch=4)
lines(summary.y[1,], type='l', col='purple', lty=3)
lines(summary.y[3,], type='l', col='purple', lty=3)
legend("topright", legend=c('Truth', 'Mean', '95% C.I.'), lty=1:3,
       col=c('black', 'blue', 'purple'), horiz = TRUE, pch=c(16, 4, NA))
dev.off()

## 1. Residuen berechnen (Wahre Werte minus geschätzte Werte)
y_fitted = summary.y[2, (p+1):n.all] 
y_true   = y.sample[(p+1):n.all]
residuals = y_true - y_fitted

## 2. Grafikfenster für 3 Diagnose-Plots aufteilen
png(filename = "./images/conjugate-residual.png", 
    width = 1920,           # Breite in Pixeln
    height = 1080,           # Höhe in Pixeln
    res = 150)              # Auflösung erhöhen für scharfen Text (DPI)

plot(residuals, type='h', col='darkblue', xlab='Zeit', ylab='Residuen',
     main='Residual over the time')
abline(h=0, col='red', lty=2)
dev.off()
# Plot: Autokorrelation der Residuen (Sollte ab Lag 1 innerhalb der blauen Linien sein)
# WICHTIG: Wenn hier Balken weit herausragen, reicht ein AR(2) Prozess nicht aus!
png(filename = "./images/conjugate-acf.png", 
    width = 1920,           # Breite in Pixeln
    height = 1080,           # Höhe in Pixeln
    res = 150)              # Auflösung erhöhen für scharfen Text (DPI)

acf(residuals, main='ACF of Residual', lag.max=20)

# Grafikfenster zurücksetzen
dev.off()
