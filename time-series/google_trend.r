# download data 
library(gtrendsR)
timeseries_data<-gtrends("time series",time="all")
plot(timeseries_data)
names(timeseries_data)

timeseries_data=timeseries_data$interest_over_time
data=list(yt=timeseries_data$hits)

library(dlm)
#model_seasonal=dlmModTrig(s=12,q=4,dV=0,dW=1)
#model_trend=dlmModPoly(order=2,dV=10,dW=rep(1,2),m0=c(40,0))
model_seasonal=dlmModTrig(s=12,q=3,dV=0,dW=1)
model_trend=dlmModPoly(order=1,dV=10,dW=rep(1,1), m0=c(40))
model=model_trend+model_seasonal
model$C0=10*diag(7)
n0=1
S0=10
k=length(model$m0)
T=length(data$yt)

Ft=array(0,c(1,k,T))
Gt=array(0,c(k,k,T))
for(t in 1:T){
   Ft[,,t]=model$FF
   Gt[,,t]=model$GG
}

source('all_dlm_functions_unknown_v.R')
source('discountfactor_selection_functions.R')

matrices=set_up_dlm_matrices_unknown_v(Ft=Ft,Gt=Gt)
initial_states=set_up_initial_states_unknown_v(model$m0,
                                               model$C0,n0,S0)

df_range=seq(0.9,1,by=0.005)

## fit discount DLM
## MSE
results_MSE <- adaptive_dlm(data, matrices, 
               initial_states, df_range,"MSE",forecast=FALSE)

## print selected discount factor
print(paste("The selected discount factor:",results_MSE$df_opt))

## retrieve filtered results
results_filtered <- results_MSE$results_filtered
ci_filtered <- get_credible_interval_unknown_v(
  results_filtered$ft,results_filtered$Qt,results_filtered$nt)

## retrieve smoothed results
results_smoothed <- results_MSE$results_smoothed
ci_smoothed <- get_credible_interval_unknown_v(
  results_smoothed$fnt, results_smoothed$Qnt, 
  results_filtered$nt[length(results_smoothed$fnt)])

## plot smoothing results 
par(mfrow=c(1,1))
index <- timeseries_data$date
plot(index, data$yt, ylab='Google hits',
     main = "Google Trends: time series", type = 'l',
     xlab = 'time', lty=3,ylim=c(0,100))
lines(index, results_smoothed$fnt, type = 'l', col='blue', 
      lwd=2)
lines(index, ci_smoothed[, 1], type='l', col='blue', lty=2)
lines(index, ci_smoothed[, 2], type='l', col='blue', lty=2)

# Plot trend and rate of change 
par(mfrow=c(2,1))
plot(index,data$yt,pch=19,cex=0.3,col='lightgray',xlab="time",
     ylab="Google hits",main="trend")
lines(index,results_smoothed$mnt[,1],lwd=2,col='magenta')
plot(index,results_smoothed$mnt[,2],col='darkblue',lwd=2,
     type='l', ylim=c(-0.6,0.6), xlab="time",
     ylab="rate of change")
abline(h=0,col='red',lty=2)

# Plot seasonal components 
par(mfrow=c(2,2))
plot(index,results_smoothed$mnt[,3],lwd=2,col="darkgreen",
     type='l', xlab="time",ylab="",main="period=12",
     ylim=c(-12,12))
plot(index,results_smoothed$mnt[,5],lwd=2,col="darkgreen",
     type='l',xlab="time",ylab="",main="period=6",
     ylim=c(-12,12))
plot(index,results_smoothed$mnt[,7],lwd=2,col="darkgreen",
     type='l', xlab="time",ylab="",main="period=4",
     ylim=c(-12,12))
plot(index,results_smoothed$mnt[,9],lwd=2,col="darkgreen",
     type='l', xlab="time",ylab="",main="period=3",
     ylim=c(-12,12))

#Estimate for the observational variance: St[T]
results_filtered$St[T]
