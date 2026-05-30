#######################
##### DLM package #####
#######################

library(dlm)
k=4
T=length(LakeHuron)-k # We take the first 
# 94 observations only as our data
index=seq(1875, 1972, length.out = length(LakeHuron))
index_filt=index[1:T]

model=dlmModPoly(order=1,dV=1,dW=0.01,m0=570,C0=10)
results_filtered_dlm=dlmFilter(LakeHuron[1:T],model)
results_smoothed_dlm=dlmSmooth(results_filtered_dlm)

plot(index_filt, LakeHuron[1:T], ylab = "level", 
     main = "Lake Huron Level",
     type='l', xlab="time",lty=3,ylim=c(575,583))
points(index_filt,LakeHuron[1:T],pch=20)
lines(index_filt,results_filtered_dlm$m[-1],col='red',lwd=2)
lines(index_filt,results_smoothed_dlm$s[-1],col='blue',lwd=2)
legend(1880,577, legend=c("filtered", "smoothed"),
       col=c("red", "blue"), lty=1, cex=0.8)
