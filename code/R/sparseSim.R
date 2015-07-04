library(xtable)

simulation.data <- read.csv("../snr_5_alpha_05.csv")
names(simulation.data)


methods <- c(paste(rep(c("simple", "forward", "strong"),2),
                   rep(c("MaxT","saturated"),each=3),
                  sep="_"),
             "knockoffs")
error.types <- c("pscreen","FWER.mod","FDR.mod","FDR.var")
error.rates <- matrix(NA,7,4,
                      dimnames=list(methods,error.types))
for(i in 1:6) {
  R <- simulation.data[[paste0(methods[i],"_R")]]
  k0 <- simulation.data[["completion_idx"]]
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
simulation.data[46, c("completion_index", "simple_saturated_R", "simple_saturated_V_model")]
simulation.data$completion_index[46]
simulation.data$simple_saturated_R[46]
simulation.data$simple_saturated_V_model[46]


error.rates["knockoffs", "pscreen"] <- with(simulation.data, mean(knockoff_R - knockoff_V == 7))
#error.rates["knockoffs", "FDR.var"] <- with(simulation.data, mean(knockoff_V / pmax(knockoff_R,1)))
o=simulation
error.rates["knockoffs", "FDR.var"] <- with(simulation.data, mean(knockoff_V / pmax(knockoff_R,1)))

## Fix this later...
colnames(error.rates) <- c("$p_{\\text{screen}}", 
                           "$\\text{FWER}_{\\text{model}}$",
                           "$\\text{FDR}_{\\text{model}}$",
                           "$\\text{FDR}_{\\text{variable}}$")
rownames(error.rates) <- c(paste(rep(c("Selected","Saturated"),each=3),
                                 rep(c("Simple", "Forward", "Strong"),2),
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


######## PLOTS

pdf("../talks/figs/simulation_snr_5_alpha_05_all.pdf",
    width=8,height=5)
par(xaxs="i",mar=c(2.1,3.1,3.1,0.7),
    mfrow=c(2,5))
for(k in 1:10) {
  main <- bquote(p*.(k))
  plot(ecdf(simulation.data[[paste0("saturated_",k)]]),xlim=c(0,1),
       main=main,xlab="",ylab="")
  plot(ecdf(simulation.data[[paste0("select_",k)]]),col="red",add=T)
  plot(ecdf(simulation.data[[paste0("nominal_",k)]]),col="blue",add=T)
  abline(0,1,lty=3,col="gray")
}
dev.off()

pdf("../talks/figs/simulation_snr_5_alpha_05_noise_var.pdf",
    width=8,height=5)
par(xaxs="i",mar=c(2.1,3.1,3.1,0.7),
    mfrow=c(2,5))
for(k in 1:10) {
  which.noise <- which(simulation.data[[paste0("var_",k)]] >= 7)
  main <- bquote(p*.(k))
  plot(ecdf(simulation.data[[paste0("saturated_pvalue_",k)]][which.noise]),xlim=c(0,1),
       main=main,xlab="",ylab="")
  plot(ecdf(simulation.data[[paste0("select_pvalue_",k)]][which.noise]),col="red",add=T)
  plot(ecdf(simulation.data[[paste0("nominal_pvalue_",k)]][which.noise]),col="blue",add=T)
  abline(0,1,lty=3,col="gray")
}
dev.off()

pdf("../talks/figs/simulation_snr_5_alpha_05_signal_var.pdf",
    width=8,height=5)
par(xaxs="i",mar=c(2.1,3.1,3.1,0.7),
    mfrow=c(2,5))
for(k in 1:10) {
  which.signal <- which(simulation.data[[paste0("var_",k)]] < 7)
  main <- bquote(p*.(k))
  plot(ecdf(simulation.data[[paste0("saturated_",k)]][which.signal]),xlim=c(0,1),
       main=main,xlab="",ylab="")
  plot(ecdf(simulation.data[[paste0("select_",k)]][which.signal]),col="red",add=T)
  plot(ecdf(simulation.data[[paste0("nominal_",k)]][which.signal]),col="blue",add=T)
  abline(0,1,lty=3,col="gray")
}
dev.off()

pdf("../talks/figs/simulation_snr_5_alpha_05_null_true.pdf",
    width=8,height=5)
par(xaxs="i",mar=c(2.1,3.1,3.1,0.7),
    mfrow=c(2,5))
for(k in 8:17) {
  which.null <- which(simulation.data$completion_index < k)
  main <- bquote(p*.(k))
  plot(ecdf(simulation.data[[paste0("saturated_",k)]][which.null]),xlim=c(0,1),
       main=main,xlab="",ylab="")
  plot(ecdf(simulation.data[[paste0("select_",k)]][which.null]),col="red",add=T)
  plot(ecdf(simulation.data[[paste0("nominal_",k)]][which.null]),col="blue",add=T)
  abline(0,1,lty=3,col="gray")
}
dev.off()

pdf("../talks/figs/simulation_snr_5_alpha_05_null_false.pdf",
    width=8,height=5)
par(xaxs="i",mar=c(2.1,3.1,3.1,0.7),
    mfrow=c(2,5))
for(k in 1:10) {
  which.nonnull <- which(simulation.data$completion_index >= k)
  main <- bquote(p*.(k))
  plot(ecdf(simulation.data[[paste0("saturated_",k)]][which.nonnull]),xlim=c(0,1),
       main=main,xlab="",ylab="")
  plot(ecdf(simulation.data[[paste0("select_",k)]][which.nonnull]),col="red",add=T)
  plot(ecdf(simulation.data[[paste0("nominal_",k)]][which.nonnull]),col="blue",add=T)
  abline(0,1,lty=3,col="gray")
}
dev.off()

