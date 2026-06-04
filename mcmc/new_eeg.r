# 1. Load data and define parameters
yt = scan('eeg.txt')
T = length(yt)
p = 8

# 2. Construct the Response vector (from t = 9 to T)
y = yt[(p+1):T] 

# 3. Construct the Design Matrix X using embed()
# embed(yt, p+1) creates a matrix where each row is c(y_t, y_{t-1}, ..., y_{t-p})
# We drop the first column (which is y_t) to get just the lags.
X = embed(yt, p+1)[, -1]

# 4. Ordinary Least Squares / Conditional MLE Estimation
XtX = t(X) %*% X
XtX_inv = solve(XtX)
phi_MLE = XtX_inv %*% t(X) %*% y 

# 5. Unbiased estimate for residual variance (v)
# Degrees of freedom = length(y) - p (number of coefficients estimated)
residuals = y - X %*% phi_MLE
s2 = sum(residuals^2) / (length(y) - p) 

# 6. Print Results
cat("\n MLE of conditional likelihood for phi:\n")
print(phi_MLE)
cat("\n Estimate for v:", s2, "\n")