B <- 1E5
rho <- .99
set.seed(1)
Z1 <- qnorm(runif(B, pnorm(4), 1)) #Gaussian conditional on being bigger than 4
Z2.indep <- rnorm(B) # second Gaussian, independent of Z1
Z2.corr <- Z1*rho + sqrt(1-rho^2)*Z2.indep # second Gaussian, correlated with Z1

pdf("../../figs/SidakNo.pdf",height=6,width=10)
par(mfrow=c(1,2))
plot(rep(Z1[1:1E3],2), c(Z2.indep[1:1E3], Z2.corr[1:1E3]), col=rep(1:2,each=1E3), cex=.5, asp=1,
     main="Truncated Gaussian", xlab="Z1", ylab="Z2")
abline(v=4,lty=2)
legend("bottomright", col=1:2, legend=c("Independent", "Correlated"), pch=1)

plot(ecdf(pmax(Z1,Z2.indep)),main="ECDF of max(|Z1|, |Z2|)",xlim=c(4,5.3))
plot(ecdf(pmax(Z1,Z2.corr)),add=T,col="red")
legend("bottomright", col=1:2, legend=c("Independent", "Correlated"), lty=1, bg="white")
dev.off()

hist(,breaks=seq(4,7,.02),border="black",xlim=c(4,5.3))
hist(pmax(Z1,Z2.corr),breaks=seq(4,7,.02),border="red",xlim=c(4,5.3),add=T)
