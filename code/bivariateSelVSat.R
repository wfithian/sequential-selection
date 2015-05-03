selected.p <- function(y) {
  total.mass <- integrate(function(x) (1-2*pnorm(-x))*dnorm(x), 
                          0, abs(y[1])+10)$value
  tail.mass  <- integrate(function(x) (1-2*pnorm(-x))*dnorm(x), 
                          abs(y[1]), abs(y[1])+10)$value
  tail.mass / total.mass
}

saturated.p <- function(y) {
  pnorm(-abs(y[1]))/pnorm(-abs(y[2]))  
}

right.wedge.mass <- function(mu) {
  # probability that y1 > y2
  p1 <- pnorm((mu[1]-mu[2])/sqrt(2))    
  # probability that y1 > -y2
  p2 <- pnorm((mu[1]+mu[2])/sqrt(2))
  p1 * p2
}

rtruncnorm.right <- function(n, mu) {
  # probability that y1 > y2
  p1 <- pnorm((mu[1]-mu[2])/sqrt(2))
  # sample y1 - y2
  ydiff <- -qnorm(p1*runif(n))*sqrt(2) + mu[1] - mu[2]
  # probability that y1 > -y2
  p2 <- pnorm((mu[1]+mu[2])/sqrt(2))
  # sample y1 + y2
  ysum <- -qnorm(p2*runif(n))*sqrt(2) + mu[1] + mu[2]
  y <- cbind((ysum+ydiff) / 2, (ysum - ydiff) / 2)
}

rtruncnorm <- function(n, mu) {
  p.tot <- right.wedge.mass(mu) + right.wedge.mass(-mu) 
  p.right <- right.wedge.mass(mu) / p.tot
  is.right.lobe <- runif(n) < p.right 
  y <- matrix(NA,n,2)
  y[is.right.lobe,] <- rtruncnorm.right(sum(is.right.lobe),mu)
  y[!is.right.lobe,] <- -rtruncnorm.right(sum(!is.right.lobe),-mu)
  y
}

plot(rtruncnorm(1000,c(.5,3)),xlim=c(-5,5),ylim=c(-5,5),asp=1)


set.seed(1)
B <- 10000
mu <- c(0,0)
pvals <- matrix(NA,B,2)
y <- rtruncnorm(B,mu)
for(b in 1:B) {
  pvals[b,] <- c(selected.p(y[b,]), saturated.p(y[b,]))
}
mean(pvals[,1]<.05) # power of selected-model test
mean(pvals[,2]<.05) # power of saturated-model test


pdf("../talks/figs/bivariateSelVSat_rocCurve.pdf",width=4,height=3)
par(xaxs="i",mar=c(4.1,4.1,0.2,0.7))
plot(ecdf(pvals[,2]),xlim=c(0,1),main="",xlab=expression(p[1]),ylab="CDF")
plot(ecdf(pvals[,1]),col="red",add=T)
legend("bottomright",legend=c("Saturated","Selected"),col=1:2,lty=1)
dev.off()

# Now, try for a grid of mu values
mu1.vals <- seq(0,5,.1)
mu2.vals <- seq(0,5,1)
mu.vals <- expand.grid(mu1=mu1.vals,mu2=mu2.vals)

set.seed(1)
B <- 10000
pvals.sat <- matrix(NA,B,nrow(mu.vals))
pvals.sel <- matrix(NA,B,nrow(mu.vals))
for(i in 1:nrow(mu.vals)) {
  y <- rtruncnorm(B,as.numeric(mu.vals[i,]))
  for(b in 1:B) {
    pvals.sel[b,i] <- selected.p(y[b,])
    pvals.sat[b,i] <- saturated.p(y[b,])
  }  
}

pow.sat <- matrix(colMeans(pvals.sat < .05),
                  length(mu1.vals),length(mu2.vals))
pow.sel <- matrix(colMeans(pvals.sel < .05),
                  length(mu1.vals),length(mu2.vals))
image(pow.sat)
image(pow.sel)

## Smooth power function
library(splines)
resp.sat <- B*(cbind(pow.sat[,5],1-pow.sat[,5]))
resp.sat <- rbind(resp.sat[nrow(resp.sat):1,], resp.sat)
resp.sel <- B*(cbind(pow.sel[,5],1-pow.sel[,5]))
resp.sel <- rbind(resp.sel[nrow(resp.sel):1,], resp.sel)
pred <- c(-rev(mu1.vals),mu1.vals)

sm.pow.sat <- glm(resp.sat ~ ns(pred,30), family=binomial)
sm.pow.sel <- glm(resp.sel ~ ns(pred,30), family=binomial)

pdf("../talks/figs/bivariateSelVSat_powCurves.pdf",width=6,height=5)
par(xaxs="i",mar=c(4.1,4.1,3.1,0.7))
plot(pred,predict(sm.pow.sat,type="response"),type="l",ylim=c(0,1),
     xlab=expression(mu[1]),ylab="Power",
     main=expression(paste("Power, ", mu[2]==4)))
lines(pred,predict(sm.pow.sel,type="response"),type="l",col="red")
abline(h=.05,lty=2)
abline(h=0:1,lty=2,col="gray")
legend("bottomright",lty=1,col=1:2,legend=c("Saturated","Selected"),
       bg="white")
dev.off()

