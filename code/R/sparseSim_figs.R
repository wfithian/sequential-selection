## plots p=40

simulation.data <- read.csv("../snr_5_alpha_05.csv")
names(simulation.data)

p=40
s=7


pdf("../../figs/simulation_snr_5_alpha_05_signal_var.pdf",
    width=8,height=2.5)
par(xaxs="i",mar=c(2.1,3.1,3.1,0.7), mfrow=c(1,4))

for(k in 4:7) {
  which.signal <- which(simulation.data[[paste0("variable_selected_",k)]] <=s)
  main <- paste("Step ",as.character(k),sep="")
  plot(ecdf(simulation.data[[paste0("nominal_pvalue_",k)]][which.signal]),xlim=c(0,1),
       main=main,xlab="",ylab="")
del=runif(length(which.signal),0,.01)
  plot(ecdf(simulation.data[[paste0("maxT_pvalue_",k)]][which.signal]+del),col="red",add=T)
    plot(ecdf(simulation.data[[paste0("maxT_unknown_pvalue_",k)]][which.signal]+del),col="green",add=T)
  plot(ecdf(simulation.data[[paste0("saturated_pvalue_",k)]][which.signal]),col="blue",add=T)
  abline(0,1,lty=3,col="gray")
if(k==4) legend (list(x=c(.3,1),y=c(.25,0)),c("Nominal","MaxZ","MaxT","Saturated"),col=c("black","red","green","blue"),lty=1,cex=.7)
}
dev.off()

#null is true


##
pdf("../../figs/simulation_snr_5_alpha_05_null_true.pdf",
    width=8,height=2.5)
par(xaxs="i",mar=c(2.1,3.1,3.1,0.7), mfrow=c(1,4))


for(k in (s+1):(s+4)) {

which.null <- which(k> simulation.data$completion_idx+1)  # correction for zero-based
  main <- paste("Step ",as.character(k),sep="")
  plot(ecdf(simulation.data[[paste0("nominal_pvalue_",k)]][which.null]),xlim=c(0,1),
       main=main,xlab="",ylab="")
plot(ecdf(simulation.data[[paste0("maxT_pvalue_",k)]][which.null]),col="red",add=T)
  plot(ecdf(simulation.data[[paste0("maxT_unknown_pvalue_",k)]][which.null]),col="green",add=T)
  plot(ecdf(simulation.data[[paste0("saturated_pvalue_",k)]][which.null]),col="blue",add=T)
  abline(0,1,lty=3,col="gray")
if(k==(s+1)) legend (list(x=c(.3,1),y=c(.25,0)),c("Nominal","MaxZ","MaxT","Saturated"),col=c("black","red","green","blue"),lty=1,cex=.7)
}
dev.off()


pdf("../../figs/simulation_snr_5_alpha_05_null_false.pdf",
    width=8,height=2.5)
par(xaxs="i",mar=c(2.1,3.1,3.1,0.7),
    mfrow=c(1,4))
for(k in 6:9) {
  which.nonnull <- which(simulation.data$completion_idx >= k)
   main <- paste("Step ",as.character(k),sep="")
  del=runif(length(which.nonnull),0,.01)
  plot(ecdf(simulation.data[[paste0("nominal_pvalue_",k)]][which.nonnull]+del),xlim=c(0,1),
       main=main,xlab="",ylab="")
plot(ecdf(simulation.data[[paste0("maxT_pvalue_",k)]][which.nonnull]+del),col="red",add=T)
  plot(ecdf(simulation.data[[paste0("maxT_unknown_pvalue_",k)]][which.nonnull]+del),col="green",add=T)
  plot(ecdf(simulation.data[[paste0("saturated_pvalue_",k)]][which.nonnull]),col="blue",add=T)
  abline(0,1,lty=3,col="gray")
if(k==6) legend (list(x=c(.3,1),y=c(.25,0)),c("Nominal","MaxZ","MaxT","Saturated"),col=c("black","red","green","blue"),lty=1,cex=.7)
}
dev.off()


pdf("../../figs/simulation_snr_5_alpha_05_noise_var.pdf",
    width=8,height=2.5)
par(xaxs="i",mar=c(2.1,3.1,3.1,0.7),
    mfrow=c(1,4))
for(k in 1:4) {
  which.noise <- which(simulation.data[[paste0("variable_selected_",k)]] >= s)
  main <- paste("Step ",as.character(k),sep="")
  plot(ecdf(simulation.data[[paste0("nominal_pvalue_",k)]][which.noise]),xlim=c(0,1),
       main=main,xlab="",ylab="")
  del=runif(length(which.noise),0,.01)
  plot(ecdf(simulation.data[[paste0("maxT_pvalue_",k)]][which.noise]+del),col="red",add=T)
  plot(ecdf(simulation.data[[paste0("maxT_unknown_pvalue_",k)]][which.noise]+del),col="green",add=T)
  plot(ecdf(simulation.data[[paste0("saturated_pvalue_",k)]][which.noise]),col="blue",add=T)
  abline(0,1,lty=3,col="gray")
  if(k==1) legend (list(x=c(.3,1),y=c(.25,0)),c("Nominal","MaxZ","MaxT","Saturated"),col=c("black","red","green","blue"),lty=1,cex=.7)
}
dev.off()

