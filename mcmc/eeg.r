yt=scan('eeg.txt')
T=length(yt)
p=8
y=rev(yt[(p+1):T]) # response
X=t(matrix(yt[rev(rep((1:p),T-p)+rep((0:(T-p-1)),rep(p,T-p)))],p,T-p));
XtX=t(X)%*%X
XtX_inv=solve(XtX)
phi_MLE=XtX_inv%*%t(X)%*%y # MLE for phi
s2=sum((y - X%*%phi_MLE)^2)/(length(y) - p) #unbiased estimate for v

cat("\n MLE of conditional likelihood for phi: ", phi_MLE, "\n",
    "Estimate for v: ", s2, "\n")

roots=1/polyroot(c(1, -phi_MLE)) # compute reciprocal characteristic roots
r=Mod(roots)
# compute moduli of reciprocal roots
lambda=2*pi/Arg(roots) # compute periods of reciprocal roots
# print results modulus and frequency by decreasing order
res <- cbind(r, abs(lambda))[order(r, decreasing=TRUE), ]
print(res[, 1])
print(res[, 2])

#[c(1,3,5,7)])
#print(cbind(r, abs(lambda))[order(r, decreasing=TRUE), ][c(2,4,6,8)])