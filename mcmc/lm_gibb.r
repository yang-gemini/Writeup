update_mu_sig2 <- function(beta_list, X, sig_init) {
  mean_b <- mean(beta_list)
  sum_X <- sum(X^2)
  return (rnorm(1, mean = mean_b, sd = sqrt(sig_init/sum_X)))
}

n_iter <- 5000
# Leere Vektoren für die gezogenen Werte (Erinnerung: leere Matrix/Vektoren erstellen)
b0_proben <- numeric(n_iter)
b1_proben <- numeric(n_iter)
b2_proben <- numeric(n_iter)
b3_proben <- numeric(n_iter)
sig2_proben <- numeric(n_iter)

# Startwerte setzen
b0 <- 0
b1 <- 0
b2 <- 0
b3 <- 0
sig2 <- 1
Y = amoda$education
X0 = amoda$v1 
X1 = amoda$v2
X2 = amoda$v2
# 3. Der Gibbs-Sampling Loop
for (i in 1:n_iter) {
  
  # Schritt 1: Bedingtes Update für beta0
  res_b0 <- Y - b1 * X1 - b2 * X2 - b3 * X3
  mean_b0 <- mean(res_b0)
  var_b0 <- sig2 / n
  b0 <- rnorm(1, mean = mean_b0, sd = sqrt(var_b0))
  
  # Schritt 2: Bedingtes Update für beta1
  res_b1 <- Y - b2 * X2 - b3 * X3 - b0
  b1 <- update_mu_sig2(res_b1, X1, sig_init)

  # Schritt 3: Bedingtes Update für beta2
  res_b2 <- Y - b1 * X1 - b3 * X3 - b0
  b2 <- update_mu_sig2(res_b2, X2, sig_init)

  # Schritt 4: Bedingtes Update für beta3
  res_b3 <- Y - b1 * X1 - b2 * X2 - b0
  b3 <- update_mu_sig2(res_b3, X3, sig_init)
 
  # Schritt 5: Bedingtes Update für die Varianz (sigma^2)
  residual <- Y - (b0 + b1 * X1 + b2 * X2 + b3 * X3)
  sum_quad_error <- sum(residual^2)
  # Ziehen aus der Inversen Gamma-Verteilung via Inverser Chi-Quadrat-Logik
  sig2 <- sum_quad_error / rchisq(1, df = length(Y) - 2)
  
  # Werte speichern
  b0_proben[i] <- b0
  b1_proben[i] <- b1
  b2_proben[i] <- b2
  b3_proben[i] <- b3
  sig2_proben[i] <- sig2
}

# 4. Ergebnisse auswerten (nach einer "Burn-in"-Phase von 1000 Schritten)
burn_in <- 1001:n_iter
cat("Geschätztes Beta0:", mean(b0_proben[burn_in]), "\n")
cat("Geschätztes Beta1:", mean(b1_proben[burn_in]), "\n")
cat("Geschätztes Beta2:", mean(b2_proben[burn_in]), "\n")
cat("Geschätztes Beta3:", mean(b3_proben[burn_in]), "\n")
cat("Geschätzte Varianz:", mean(sig2_proben[burn_in]), "\n")