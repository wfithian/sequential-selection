B <- 10000
mu <- c(4,4)
pvals <- NULL
for(b in 1:B) {
  y <- mu + rnorm(2)
  if(abs(y[1]) > abs(y[2])) {
    pvals <- rbind(pvals, c(
      integrate(function(x) abs(x)*dnorm(x)/integrate(function(u) u*dnorm(u), 0, 10)$value, 
                abs(y[1]), 10)$value,
      pnorm(-abs(y[1]))/pnorm(-abs(y[2]))
    ))
  }
}
mean(pvals[,1]<.05) # power of selected-model test
mean(pvals[,2]<.05) # power of saturated-model test

pdf("../talks/figs/bivariateSelVSat.pdf",width=4,height=3)
par(xaxs="i",mar=c(4.1,4.1,0.2,0.7))
plot(ecdf(pvals[,2]),xlim=c(0,1),main="",xlab=expression(p[1]),ylab="CDF")
plot(ecdf(pvals[,1]),col="red",add=T)
legend("bottomright",legend=c("Saturated","Selected"),col=1:2,lty=1)
dev.off()
