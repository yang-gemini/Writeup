library(MCMCpack)
library(mvtnorm)

model_comp_mix = function(tot_num_comp, y, p_order) {
  
  ## Hyperparameters und Dimensionen zuweisen
  p = p_order
  n.all = length(y)
  K = tot_num_comp
  
  m0 = matrix(rep(0, p), ncol = 1)
  C0 = 10 * diag(p)
  C0.inv = 0.1 * diag(p)
  n0 = 0.02
  d0 = 0.02
  a = rep(1, K)
  
  Y = matrix(y[(p + 1):n.all], ncol = 1)
  
  # Dynamischer Aufbau der Designmatrix Fmtx passend zur AR-Ordnung p
  Fmtx_list = lapply(1:p, function(i) y[(p + 1 - i):(n.all - i)])
  Fmtx = matrix(unlist(Fmtx_list), nrow = p, byrow = TRUE)
  n = length(Y)
  
  #### Gibbs-Sampler Teilfunktionen
  sample_omega = function(L.cur) {
    n.vec = sapply(1:K, function(k) { sum(L.cur == k) })
    rdirichlet(1, a + n.vec)
  }
  
  sample_L_one = function(beta.cur, omega.cur, nu.cur, y.cur, Fmtx.cur) {
    prob_k = function(k) {
      beta.use = beta.cur[((k - 1) * p + 1):(k * p)]
      omega.cur[k] * dnorm(y.cur, mean = sum(beta.use * Fmtx.cur), sd = sqrt(nu.cur))
    }
    prob.vec = sapply(1:K, prob_k)
    if(sum(prob.vec) == 0) prob.vec = rep(1, K)
    L.sample = sample(1:K, 1, prob = prob.vec / sum(prob.vec))
    return(L.sample)
  }
  
  sample_L = function(y, x, beta.cur, omega.cur, nu.cur) {
    L.new = sapply(1:n, function(j) { sample_L_one(beta.cur, omega.cur, nu.cur, y.cur = y[j, ], Fmtx.cur = x[, j]) })
    return(L.new)
  }
  
  sample_nu = function(L.cur, beta.cur) {
    n.star = n0 + n + p * K
    err.y = function(idx) {
      L.use = L.cur[idx]
      beta.use = beta.cur[((L.use - 1) * p + 1):(L.use * p)]
      err = Y[idx, ] - sum(Fmtx[, idx] * beta.use)
      return(err^2)
    }
    err.beta = function(k.cur) {
      beta.use = beta.cur[((k.cur - 1) * p + 1):(k.cur * p)]
      beta.use.minus.m0 = matrix(beta.use, ncol = 1) - m0
      t(beta.use.minus.m0) %*% C0.inv %*% beta.use.minus.m0
    }
    d.star = d0 + sum(sapply(1:n, err.y)) + sum(sapply(1:K, err.beta))
    1 / rgamma(1, shape = n.star / 2, rate = d.star / 2)
  }
  
  sample_beta = function(k, L.cur, nu.cur) {
    idx.select = (L.cur == k)
    n.k = sum(idx.select)
    if (n.k == 0) {
      m.k = m0
      C.k = C0
    } else {
      y.tilde.k = matrix(Y[idx.select, ], ncol = 1)
      Fmtx.tilde.k = Fmtx[, idx.select, drop = FALSE] # drop=FALSE fixiert Fehler bei n.k == 1
      e.k = y.tilde.k - t(Fmtx.tilde.k) %*% m0
      Q.k = t(Fmtx.tilde.k) %*% C0 %*% Fmtx.tilde.k + diag(n.k)
      Q.k.inv = chol2inv(chol(Q.k))
      A.k = C0 %*% Fmtx.tilde.k %*% Q.k.inv
      m.k = m0 + A.k %*% e.k
      C.k = C0 - A.k %*% Q.k %*% t(A.k)
    }
    rmvnorm(1, m.k, nu.cur * C.k)
  }
  
  nsim = 10000
  beta.mtx = matrix(0, nrow = p * K, ncol = nsim)
  L.mtx = matrix(0, nrow = n, ncol = nsim)
  omega.mtx = matrix(0, nrow = K, ncol = nsim)
  nu.vec = rep(0, nsim)
  
  # Sinnvolle Zufallsinitialisierung statt reiner Nullen
  beta.cur = rnorm(p * K, 0, 0.5)
  L.cur = sample(1:K, n, replace = TRUE)
  omega.cur = rep(1 / K, K)
  nu.cur = var(Y)
  
  set.seed(42) # Einmaliger globaler Seed für Stabilität
  for (i in 1:nsim) {
    omega.cur = sample_omega(L.cur)
    omega.mtx[, i] = omega.cur
    
    L.cur = sample_L(Y, Fmtx, beta.cur, omega.cur, nu.cur)
    L.mtx[, i] = L.cur
    
    nu.cur = sample_nu(L.cur, beta.cur)
    nu.vec[i] = nu.cur
    
    beta.cur = as.vector(sapply(1:K, function(k) { sample_beta(k, L.cur, nu.cur) }))
    beta.mtx[, i] = beta.cur
    
    if (i %% 2000 == 0) print(paste("K =", K, "| Iteration:", i))
  }
  
  #### DIC Berechnung
  cal_log_likelihood_mix_one = function(idx, beta, nu, omega) {
    norm.lik = rep(0, K)
    for (k.cur in 1:K) {
      mean.norm = sum(Fmtx[, idx] * beta[((k.cur - 1) * p + 1):(k.cur * p)])
      norm.lik[k.cur] = dnorm(Y[idx, 1], mean.norm, sqrt(nu), log = FALSE)
    }
    log.lik = log(sum(norm.lik * omega))
    return(log.lik)
  }
  
  cal_log_likelihood_mix = function(beta, nu, omega) {
    sum(sapply(1:n, function(idx) { cal_log_likelihood_mix_one(idx = idx, beta = beta, nu = nu, omega = omega) }))
  }
  
  sample.select.idx = seq(5001, 10000, by = 1)
  
  beta.mix = rowMeans(beta.mtx[, sample.select.idx, drop=FALSE])
  nu.mix = mean(nu.vec[sample.select.idx])
  omega.mix = rowMeans(omega.mtx[, sample.select.idx, drop=FALSE])
  
  log.lik.bayes.mix = cal_log_likelihood_mix(beta.mix, nu.mix, omega.mix)
  post.log.lik.mix = sapply(sample.select.idx, function(k) { cal_log_likelihood_mix(beta.mtx[, k], nu.vec[k], omega.mtx[, k]) })
  E.post.log.lik.mix = mean(post.log.lik.mix)
  
  p_DIC.mix = 2 * (log.lik.bayes.mix - E.post.log.lik.mix)
  DIC.mix = -2 * log.lik.bayes.mix + 2 * p_DIC.mix
  
  return(DIC.mix)
}

# --- DATEN UND LOOPS ---
# Annahme: 'earthquakes.txt' existiert im Ordner
dat <- read.delim("earthquakes.txt")
y = as.numeric(dat$Quakes)

# Wir definieren p = 3 als AR-Ordnung für das Modell
ar_p = 3 

# Schleife über 2 bis 5 Komponenten ausführen
components <- 2:5
mix.model.all = sapply(components, function(k) model_comp_mix(tot_num_comp = k, y = y, p_order = ar_p))

# --- PLOT DER ERGEBNISSE ---
plot(components, mix.model.all, type = "b", pch = 19, col = "darkblue",
     xlab = "Anzahl der Mischungskomponenten (K)", 
     ylab = "Deviance Information Criterion (DIC)",
     main = "Modellauswahl via DIC für Mixture AR-Modell",
     xaxt = "n") # Standard-Achse deaktivieren für saubere Beschriftung

axis(1, at = components) # Saubere Ticks bei 2, 3, 4, 5 einfügen
grid(col = "gray", lty = "dotted") # Hintergrundgitter für bessere Lesbarkeit

# Markiert die beste Komponente (niedrigster DIC ist am besten)
best_k <- components[which.min(mix.model.all)]
points(best_k, min(mix.model.all), col = "red", cex = 2, lwd = 2)
text(best_k, min(mix.model.all), labels = paste("Minimum K =", best_k), pos = 3, col = "red")

