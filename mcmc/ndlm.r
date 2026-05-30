#################################################
##### Univariate DLM: Known, constant variances
#################################################
set_up_dlm_matrices <- function(FF, GG, VV, WW){
  return(list(FF=FF, GG=GG, VV=VV, WW=WW))
}

set_up_initial_states <- function(m0, C0){
  return(list(m0=m0, C0=C0))
}

### forward update equations ###
forward_filter <- function(data, matrices, initial_states){
  ## retrieve dataset
  yt <- data$yt
  T <- length(yt)
  
  ## retrieve a set of quadruples 
  # FF, GG, VV, WW are scalar
  FF <- matrices$FF  
  GG <- matrices$GG
  VV <- matrices$VV
  WW <- matrices$WW
  
  ## retrieve initial states
  m0 <- initial_states$m0
  C0 <- initial_states$C0
  
  ## create placeholder for results
  d <- dim(GG)[1]
  at <- matrix(NA, nrow=T, ncol=d)
  Rt <- array(NA, dim=c(d, d, T))
  ft <- numeric(T)
  Qt <- numeric(T)
  mt <- matrix(NA, nrow=T, ncol=d)
  Ct <- array(NA, dim=c(d, d, T))
  et <- numeric(T)
  
  
  for(i in 1:T){
    # moments of priors at t
    if(i == 1){
      at[i, ] <- GG %*% t(m0)
      Rt[, , i] <- GG %*% C0 %*% t(GG) + WW
      Rt[,,i] <- 0.5*Rt[,,i]+0.5*t(Rt[,,i]) 
    }else{
      at[i, ] <- GG %*% t(mt[i-1, , drop=FALSE])
      Rt[, , i] <- GG %*% Ct[, , i-1] %*% t(GG) + WW
      Rt[,,i] <- 0.5*Rt[,,i]+0.5*t(Rt[,,i]) 
    }
    
    # moments of one-step forecast:
    ft[i] <- t(FF) %*% (at[i, ]) 
    Qt[i] <- t(FF) %*% Rt[, , i] %*% FF + VV
    
    # moments of posterior at t:
    At <- Rt[, , i] %*% FF / Qt[i]
    et[i] <- yt[i] - ft[i]
    mt[i, ] <- at[i, ] + t(At) * et[i]
    Ct[, , i] <- Rt[, , i] - Qt[i] * At %*% t(At)
    Ct[,,i] <- 0.5*Ct[,,i] + 0.5*t(Ct[,,i]) 
  }
  cat("Forward filtering is completed!") # indicator of completion
  return(list(mt = mt, Ct = Ct, at = at, Rt = Rt, 
              ft = ft, Qt = Qt))
}


forecast_function <- function(posterior_states, k, matrices){
  
  ## retrieve matrices
  FF <- matrices$FF
  GG <- matrices$GG
  WW <- matrices$WW
  VV <- matrices$VV
  mt <- posterior_states$mt
  Ct <- posterior_states$Ct
  
  ## set up matrices
  T <- dim(mt)[1] # time points
  d <- dim(mt)[2] # dimension of state parameter vector
  
  ## placeholder for results
  at <- matrix(NA, nrow = k, ncol = d)
  Rt <- array(NA, dim=c(d, d, k))
  ft <- numeric(k)
  Qt <- numeric(k)
  
  
  for(i in 1:k){
    ## moments of state distribution
    if(i == 1){
      at[i, ] <- GG %*% t(mt[T, , drop=FALSE])
      Rt[, , i] <- GG %*% Ct[, , T] %*% t(GG) + WW
      Rt[,,i] <- 0.5*Rt[,,i]+0.5*t(Rt[,,i]) 
    }else{
      at[i, ] <- GG %*% t(at[i-1, , drop=FALSE])
      Rt[, , i] <- GG %*% Rt[, , i-1] %*% t(GG) + WW
      Rt[,,i] <- 0.5*Rt[,,i]+0.5*t(Rt[,,i]) 
    }
    
    ## moments of forecast distribution
    ft[i] <- t(FF) %*% t(at[i, , drop=FALSE])
    Qt[i] <- t(FF) %*% Rt[, , i] %*% FF + VV
  }
  cat("Forecasting is completed!") # indicator of completion
  return(list(at=at, Rt=Rt, ft=ft, Qt=Qt))
}

## obtain 95% credible interval
get_credible_interval <- function(mu, sigma2, 
                                  quantile = c(0.025, 0.975)){
  z_quantile <- qnorm(quantile)
  bound <- matrix(0, nrow=length(mu), ncol=2)
  bound[, 1] <- mu + z_quantile[1]*sqrt(as.numeric(sigma2)) # lower bound
  bound[, 2] <- mu + z_quantile[2]*sqrt(as.numeric(sigma2)) # upper bound
  return(bound)
}

### smoothing equations ###
backward_smoothing <- function(data, matrices, 
                               posterior_states){
  ## retrieve data 
  yt <- data$yt
  T <- length(yt) 
  
  ## retrieve matrices
  FF <- matrices$FF
  GG <- matrices$GG
  
  ## retrieve matrices
  mt <- posterior_states$mt
  Ct <- posterior_states$Ct
  at <- posterior_states$at
  Rt <- posterior_states$Rt
  
  ## create placeholder for posterior moments 
  mnt <- matrix(NA, nrow = dim(mt)[1], ncol = dim(mt)[2])
  Cnt <- array(NA, dim = dim(Ct))
  fnt <- numeric(T)
  Qnt <- numeric(T)
  for(i in T:1){
    # moments for the distributions of the state vector given D_T
    if(i == T){
      mnt[i, ] <- mt[i, ]
      Cnt[, , i] <- Ct[, , i]
      Cnt[, , i] <- 0.5*Cnt[, , i] + 0.5*t(Cnt[, , i]) 
    }else{
      inv_Rtp1<-solve(Rt[,,i+1])
      Bt <- Ct[, , i] %*% t(GG) %*% inv_Rtp1
      mnt[i, ] <- mt[i, ] + Bt %*% (mnt[i+1, ] - at[i+1, ])
      Cnt[, , i] <- Ct[, , i] + Bt %*% (Cnt[, , i + 1] - Rt[, , i+1]) %*% t(Bt)
      Cnt[,,i] <- 0.5*Cnt[,,i] + 0.5*t(Cnt[,,i]) 
    }
    # moments for the smoothed distribution of the mean response of the series
    fnt[i] <- t(FF) %*% t(mnt[i, , drop=FALSE])
    Qnt[i] <- t(FF) %*% t(Cnt[, , i]) %*% FF
  }
  cat("Backward smoothing is completed!")
  return(list(mnt = mnt, Cnt = Cnt, fnt=fnt, Qnt=Qnt))
}

####################### Example: Lake Huron Data ######################
plot(LakeHuron) # 98 observations total 
k=4
T=length(LakeHuron)-k # We take the first 
# 94 observations only as our data
ts_data=LakeHuron[1:T]
ts_validation_data <- LakeHuron[(T+1):98]

data <- list(yt = ts_data)

## set up dlm matrices
GG <- as.matrix(1)
FF <- as.matrix(1)
VV <- as.matrix(1)
WW <- as.matrix(1)
m0 <- as.matrix(570)
C0 <- as.matrix(1e4)

## wrap up all matrices and initial values
matrices <- set_up_dlm_matrices(FF, GG, VV, WW)
initial_states <- set_up_initial_states(m0, C0)

## filtering and smoothing 
results_filtered <- forward_filter(data, matrices, 
                                   initial_states)
results_smoothed <- backward_smoothing(data, matrices, 
                                       results_filtered)

index=seq(1875, 1972, length.out = length(LakeHuron))
index_filt=index[1:T]


par(mfrow=c(2,1))
plot(index, LakeHuron, main = "Lake Huron Level ",type='l',
     xlab="time",ylab="feet",lty=3,ylim=c(575,583))
points(index,LakeHuron,pch=20)
lines(index_filt, results_filtered$mt, type='l', 
      col='red',lwd=2)
lines(index_filt, results_smoothed$mnt, type='l', 
      col='blue',lwd=2)


# Now let's look at the DLM package 
library(dlm)
model=dlmModPoly(order=1,dV=1,dW=1,m0=570,C0=1e4)
results_filtered_dlm=dlmFilter(LakeHuron[1:T],model)
results_smoothed_dlm=dlmSmooth(results_filtered_dlm)

plot(index_filt, LakeHuron[1:T], ylab = "level", 
     main = "Lake Huron Level",
     type='l', xlab="time",lty=3,ylim=c(575,583))
points(index_filt,LakeHuron[1:T],pch=20)
lines(index_filt,results_filtered_dlm$m[-1],col='red',lwd=2)
lines(index_filt,results_smoothed_dlm$s[-1],col='blue',lwd=2)

# Similarly, for the second order polynomial and the co2 data:
T=length(co2)
data=list(yt = co2)

FF <- (as.matrix(c(1,0)))
GG <- matrix(c(1,1,0,1),ncol=2,byrow=T)
VV <- as.matrix(200)
WW <- 0.01*diag(2)
m0 <- t(as.matrix(c(320,0)))
C0 <- 10*diag(2)

## wrap up all matrices and initial values
matrices <- set_up_dlm_matrices(FF,GG, VV, WW)
initial_states <- set_up_initial_states(m0, C0)

## filtering and smoothing 
results_filtered <- forward_filter(data, matrices, 
                                   initial_states)
results_smoothed <- backward_smoothing(data, matrices, 
                                       results_filtered)

#### Now, using the DLM package: 
model=dlmModPoly(order=2,dV=200,dW=0.01*rep(1,2),
                 m0=c(320,0),C0=10*diag(2))
# filtering and smoothing 
results_filtered_dlm=dlmFilter(data$yt,model)
results_smoothed_dlm=dlmSmooth(results_filtered_dlm)

par(mfrow=c(2,1))
plot(as.vector(time(co2)),co2,type='l',xlab="time",
     ylim=c(300,380))
lines(as.vector(time(co2)),results_filtered$mt[,1],
      col='red',lwd=2)
lines(as.vector(time(co2)),results_smoothed$mnt[,1],
      col='blue',lwd=2)

plot(as.vector(time(co2)),co2,type='l',xlab="time",
     ylim=c(300,380))
lines(as.vector(time(co2)),results_filtered_dlm$m[-1,1],
      col='red',lwd=2)
lines(as.vector(time(co2)),results_smoothed_dlm$s[-1,1],
      col='blue',lwd=2)


