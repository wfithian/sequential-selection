## plots p=40

simulation.data <- read.csv("../snr_5_alpha_05.csv")
names(simulation.data)

p=40
s=7


pdf("../../figs/simulation_snr_5_alpha_05_signal_var.pdf",
    width=8,height=5)
par(xaxs="i",mar=c(2.1,3.1,3.1,0.7), mfrow=c(2,5))
for(k in 1:10) {
  which.signal <- which(simulation.data[[paste0("variable_selected_",k)]] <=s)
  main <- paste("Step ",as.character(k),sep="")
  plot(ecdf(simulation.data[[paste0("nominal_pvalue_",k)]][which.signal]),xlim=c(0,1),
       main=main,xlab="",ylab="")
del=runif(length(which.signal),0,.01)
  plot(ecdf(simulation.data[[paste0("maxT_pvalue_",k)]][which.signal]+del),col="red",add=T)
    plot(ecdf(simulation.data[[paste0("maxT_unknown_pvalue_",k)]][which.signal]+del),col="green",add=T)
  plot(ecdf(simulation.data[[paste0("saturated_pvalue_",k)]][which.signal]),col="blue",add=T)
  abline(0,1,lty=3,col="gray")
if(k==1) legend (list(x=c(.3,1),y=c(.25,0)),c("Nominal","MaxT","Maxt/Unkn","Saturated"),col=c("black","red","green","blue"),lty=1,cex=.7)
}
dev.off()

#null is true


##
pdf("../../figs/simulation_snr_5_alpha_05_null_true.pdf",
    width=8,height=5)
par(xaxs="i",mar=c(2.1,3.1,3.1,0.7), mfrow=c(2,5))


for(k in (s+1):(s+10)) {

which.null <- which(k> simulation.data$completion_idx+1)  # correction for zero-based
  main <- paste("Step ",as.character(k),sep="")
  plot(ecdf(simulation.data[[paste0("nominal_pvalue_",k)]][which.null]),xlim=c(0,1),
       main=main,xlab="",ylab="")
plot(ecdf(simulation.data[[paste0("maxT_pvalue_",k)]][which.null]),col="red",add=T)
  plot(ecdf(simulation.data[[paste0("maxT_unknown_pvalue_",k)]][which.null]),col="green",add=T)
  plot(ecdf(simulation.data[[paste0("saturated_pvalue_",k)]][which.null]),col="blue",add=T)
  abline(0,1,lty=3,col="gray")
if(k==(s+1)) legend (list(x=c(.3,1),y=c(.25,0)),c("Nominal","MaxT","Maxt/Unkn","Saturated"),col=c("black","red","green","blue"),lty=1,cex=.7)
}
dev.off()


pdf("../../figs/simulation_snr_5_alpha_05_null_false.pdf",
    width=8,height=5)
par(xaxs="i",mar=c(2.1,3.1,3.1,0.7),
    mfrow=c(2,5))
for(k in 1:10) {
  which.nonnull <- which(simulation.data$completion_idx >= k)
   main <- paste("Step ",as.character(k),sep="")
  del=runif(length(which.nonnull),0,.01)
  plot(ecdf(simulation.data[[paste0("nominal_pvalue_",k)]][which.nonnull]+del),xlim=c(0,1),
       main=main,xlab="",ylab="")
plot(ecdf(simulation.data[[paste0("maxT_pvalue_",k)]][which.nonnull]+del),col="red",add=T)
  plot(ecdf(simulation.data[[paste0("maxT_unknown_pvalue_",k)]][which.nonnull]+del),col="green",add=T)
  plot(ecdf(simulation.data[[paste0("saturated_pvalue_",k)]][which.nonnull]),col="blue",add=T)
  abline(0,1,lty=3,col="gray")
if(k==1) legend (list(x=c(.3,1),y=c(.25,0)),c("Nominal","MaxT","Maxt/Unkn","Saturated"),col=c("black","red","green","blue"),lty=1,cex=.7)
}
dev.off()


pdf("../../figs/simulation_snr_5_alpha_05_noise_var.pdf",
    width=8,height=5)
par(xaxs="i",mar=c(2.1,3.1,3.1,0.7),
    mfrow=c(2,5))
for(k in 1:10) {
  which.noise <- which(simulation.data[[paste0("variable_selected_",k)]] >= s)
  main <- paste("Step ",as.character(k),sep="")
  plot(ecdf(simulation.data[[paste0("nominal_pvalue_",k)]][which.noise]),xlim=c(0,1),
       main=main,xlab="",ylab="")
  del=runif(length(which.noise),0,.01)
  plot(ecdf(simulation.data[[paste0("maxT_pvalue_",k)]][which.noise]+del),col="red",add=T)
  plot(ecdf(simulation.data[[paste0("maxT_unknown_pvalue_",k)]][which.noise]+del),col="green",add=T)
  plot(ecdf(simulation.data[[paste0("saturated_pvalue_",k)]][which.noise]),col="blue",add=T)
  abline(0,1,lty=3,col="gray")
  if(k==1) legend (list(x=c(.3,1),y=c(.25,0)),c("Nominal","MaxT","Maxt/Unkn","Saturated"),col=c("black","red","green","blue"),lty=1,cex=.7)
}
dev.off()

#plots for p=200

simulation.data <- read.csv("../snr_5_alpha_05_p200.csv")
names(simulation.data)

p=200
s=7

pdf("../../figs/simulation_snr_5_alpha_05_signal_var_p200.pdf",
    width=8,height=5)
par(xaxs="i",mar=c(2.1,3.1,3.1,0.7), mfrow=c(2,5))
for(k in 1:10) {
  which.signal <- which(simulation.data[[paste0("variable_selected_",k)]] <= s)
 main <- paste("Step ",as.character(k),sep="")
  plot(ecdf(simulation.data[[paste0("nominal_pvalue_",k)]][which.signal]),xlim=c(0,1),
       main=main,xlab="",ylab="")
del=runif(length(which.signal),0,.01)
  plot(ecdf(simulation.data[[paste0("maxT_pvalue_",k)]][which.signal]+del),col="red",add=T)
   plot(ecdf(simulation.data[[paste0("maxT_unknown_pvalue_",k)]][which.signal]+del),col="green",add=T)
  plot(ecdf(simulation.data[[paste0("saturated_pvalue_",k)]][which.signal]),col="blue",add=T)
  abline(0,1,lty=3,col="gray")
  if(k==1) legend (list(x=c(.3,1),y=c(.25,0)),c("Nominal","MaxT","Maxt/Unkn","Saturated"),col=c("black","red","green","blue"),lty=1,cex=.7)
}

dev.off()



# model selection stuff

library(xtable)

simulation.data <- read.csv("../snr_5_alpha_05.csv")
names(simulation.data)


#methods <- c(paste(rep(c("simple", "forward", "strong"),2),
#                   rep(c("maxT_forward","saturated"),each=3),
#                  sep="_"),
 #            "knockoffs")
methods=c("nominal_simple","nominal_forward","nominal_strong",
         "maxT_simple","maxT_forward","maxT_strong",
         "maxT_identify_simple", "maxT_identify_forward", "maxT_identify_strong",
         "saturated_simple" ,"saturated_forward", "saturated_strong",
            "knockoff")

error.types <- c("pscreen","FWER.mod","FDR.mod","FDR.var")
error.rates <- matrix(NA,13,4,
                      dimnames=list(methods,error.types))
for(i in 1:12) {
  R <- simulation.data[[paste0(methods[i],"_R")]]
  k0 <- simulation.data[["completion_idx"]]+1
  V.mod <- pmax(0,R - k0)
  V.var <- simulation.data[[paste0(methods[i],"_V_var")]]
  error.rates[i,"pscreen"] <- mean(R >= k0)
  error.rates[i,"FWER.mod"] <- mean(V.mod > 1)
  error.rates[i,"FDR.mod"] <- mean(V.mod / pmax(R,1))
  error.rates[i,"FDR.var"] <- mean(V.var / pmax(R,1))
}

# What's wrong here?
with(simulation.data, table(knockoff_R - knockoff_V))
with(simulation.data, summary(simple_saturated_screen))
with(simulation.data, summary(simple_saturated_V_model))
with(simulation.data, which.max(simple_saturated_V_model))
simulation.data[46, c("completion_idx", "simple_saturated_R", "simple_saturated_V_model")]
simulation.data$completion_index[46]
simulation.data$simple_saturated_R[46]
simulation.data$simple_saturated_V_model[46]


error.rates["knockoff", "pscreen"] <- with(simulation.data, mean(knockoff_R - knockoff_V == 7))
error.rates["knockoff", "FDR.var"] <- with(simulation.data, mean(knockoff_V / pmax(knockoff_R,1)))
                                           

colnames(error.rates) <- c("$p_{\\text{screen}}", 
                           "$\\text{FWER}_{\\text{model}}$",
                           "$\\text{FDR}_{\\text{model}}$",
                           "$\\text{FDR}_{\\text{variable}}$")
rownames(error.rates) <- c(paste(rep(c("Nominal","MaxT","Maxt-Identify","Saturated"),each=3),
                                 rep(c("Simple", "Forward", "Strong"),3),
                                 sep=", "),
                           "Knockoffs")
xtable(error.rates, digits=3)


# variable-wise FDR for selective inference methods
FDR.var <- matrix(NA,3,3,
                  dimnames=list(c("saturated", "selected", "nominal"),
                                c("simple", "forward", "strong")))
for(method in rownames(FDR.var)) {
  for(rule in colnames(FDR.var)) {
    V.var <- simulation.data[[paste(rule,method,"V_var",sep="_")]]
    R <- simulation.data[[paste(rule,method,"R",sep="_")]]
    FDR.var[method,rule] <- mean(V.var/pmax(R,1))
  }
}

# pretty LaTeX table
xtable(FDR.var,digits=3)

# Knockoff FDR
mean(simulation.data$knockoff_V/pmax(simulation.data$knockoff_R,1))


# model-wise FDR for selective inference methods
FDR.mod <- matrix(NA,3,3,
                  dimnames=list(c("saturated", "selected", "nominal"),
                                c("simple", "forward", "strong")))
for(method in rownames(FDR.var)) {
  for(rule in colnames(FDR.var)) {
    V.mod <- simulation.data[[paste(rule,method,"V_model",sep="_")]]
    R <- simulation.data[[paste(rule,method,"R",sep="_")]]
    FDR.mod[method,rule] <- mean(V.mod/pmax(R,1))
  }
}

# pretty LaTeX table
xtable(FDR.mod,digits=3)
