#########################################
# Make plots of the conditioning sets and null distributions
#########################################

setwd("~/Dropbox/git/sequential-selection/code")

# Conditioning Sets
pdf("../figs/bivariateSelVSat_condSets.pdf",height=4.3,width=4)
par(mar=c(4.1,4.1,3.1,0.1))
y <- t(c(2.9,2.5))
plot(y,xlim=c(-5,5),ylim=c(-5,5),xlab=expression(Y[1]),ylab=expression(Y[2]),
     main="Conditioning Sets")#"Full vs. Reduced Model: First Step",asp=1)
polygon(c(0,10,10),c(0,10,-10),lty=2,col="#F4E918")
polygon(c(0,-10,-10),c(0,10,-10),lty=2,col="#F4E918")
abline(h=0)
abline(v=0)
text(2,.5,expression(A[1]))
text(y+c(.3,-.4),labels="Y")
lines(c(y[2],10),c(y[2],y[2]),lwd=2,col="brown")
lines(c(-y[2],-10),c(y[2],y[2]),lwd=2,col="brown")
points(y,pch=16)
dev.off()

# Null Distributions
pdf("../figs/bivariateSelVSat_nullDists.pdf",height=4.3,width=4)
par(mar=c(4.1,4.1,3.1,0.1),yaxs="i")
x <- seq(-6,6,.01)
plot(x,(abs(x)>2.5)*dnorm(x)/2/pnorm(-2.5),ylim=c(0,1.4),lty=1,
     col="brown",type="l",
     main="Conditional Null Distributions",
     ylab="Density",xlab=expression(Y[1]))
selected.model.density <- (1-2*pnorm(-abs(x)))*dnorm(x) * 2
polygon(c(x,0),c(selected.model.density,0),lty=2,col="#F4E918")
lines(x,(abs(x)>2.5)*dnorm(x)/2/pnorm(-2.5),col="brown")
legend("topleft",legend=c("Saturated Model","Selected Model","Observed Value"),lty=1:3,bg="white", col=c("brown","black","black"))
abline(v=2.9,lty=3)
dev.off()


####################################
##  Code to produce p-values
####################################

# y1 is a (vector of) values of y1
selected.p <- function(y1) {
  ## total probability of right wedge, under the null
  total.mass <- integrate(function(x) (1-2*pnorm(-x))*dnorm(x), 
                          0, max(abs(y1))+10)$value

  tail.mass <- numeric(length(y1))
  for(i in 1:length(y1)) {
    ## probability of right wedge *and* Y > y1, under the null
    tail.mass[i]  <- integrate(function(x) (1-2*pnorm(-x))*dnorm(x), 
                          abs(y1[i]), max(abs(y1))+10)$value
  }
  tail.mass / total.mass
}

## Rejection threshold for selected test at alpha=0.05
selected.threshold.05 <- optimize(function(y) (selected.p(y)-0.05)^2, interval=c(0,4))$minimum

## Test at level 0.05 under selected model
## Returns TRUE if reject, FALSE o/w
selected.test <- function(y1) {
  return(abs(y1) > selected.threshold.05)
}

saturated.p <- function(y1, y2) {
  pnorm(-abs(y1))/pnorm(-abs(y2))  
}

saturated.test <- function(y1, y2) {
  return(saturated.p(y1,y2) < 0.05)
}


## p-Values when Y=(2.9,2.5)
selected.p(2.9)
saturated.p(2.9,2.5)



####################################
##  Code to simulate from Y|A 
####################################

## Utility functions
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





####################################
##  Simulations
####################################


set.seed(1)
B <- 1E4
mu <- c(4,4)
pvals <- matrix(NA,B,2)
y <- rtruncnorm(B,mu)
pvals <- cbind(selected.p(y[,1]), saturated.p(y[,1], y[,2]))
mean(pvals[,1]<.05) # power of selected-model test
mean(pvals[,2]<.05) # power of saturated-model test

pdf("../figs/bivariateSelVSat_rocCurve.pdf",width=4,height=3)
par(xaxs="i",mar=c(4.1,4.1,0.2,0.7))
plot(ecdf(pvals[,2]),xlim=c(0,1),main="",xlab=expression(p[1]),ylab="CDF")
plot(ecdf(pvals[,1]),col="red",add=T)
lines(c(max(pvals[,1]),1),c(1,1),col="red")
abline(h=1,col="gray",lty=2)
legend("bottomright",legend=c("Saturated","Selected"),col=1:2,lty=1,bg="white")
dev.off()

# Now, try for a grid of mu values
mu1.vals <- c(seq(0,.1,.02),seq(.2,8,.2))
mu2.vals <- c(0,4)#seq(0,5,1)
set.seed(1)
B <- 1E5
nrej.sel <- nrej.sat <- 
    matrix(NA,length(mu1.vals),length(mu2.vals),
           dimnames=list(as.character(mu1.vals),
                         as.character(mu2.vals)))
for(mu2 in mu2.vals) {
  for(mu1 in mu1.vals) {
    y <- rtruncnorm(B,c(mu1,mu2))
    nrej.sel[as.character(mu1), as.character(mu2)] <-
      sum(selected.test(y[,1]))
    nrej.sat[as.character(mu1), as.character(mu2)] <-
      sum(saturated.test(y[,1], y[,2]))
  }
  cat("finished for mu2 = ", mu2, "\n")
}


library(splines)

#mu2 = 4

## Smooth power function using logistic regression with splines
resp.sat <- cbind(nrej.sat[,"4"],B-nrej.sat[,"4"])
resp.sat <- rbind(resp.sat[nrow(resp.sat):1,], resp.sat)
resp.sel <- cbind(nrej.sel[,"4"],B-nrej.sel[,"4"])
resp.sel <- rbind(resp.sel[nrow(resp.sel):1,], resp.sel)
pred <- c(-rev(mu1.vals),mu1.vals)

sm.pow.sat <- glm(resp.sat ~ ns(pred,40), family=binomial)
sm.pow.sel <- glm(resp.sel ~ ns(pred,40), family=binomial)

pdf("../figs/bivariateSelVSat_powCurves_4.pdf",width=4.5,height=4.3)
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

# mu2 = 0

resp.sat <- cbind(nrej.sat[,"0"],B-nrej.sat[,"0"])
resp.sat <- rbind(resp.sat[nrow(resp.sat):1,], resp.sat)
resp.sel <- cbind(nrej.sel[,"0"],B-nrej.sel[,"0"])
resp.sel <- rbind(resp.sel[nrow(resp.sel):1,], resp.sel)
pred <- c(-rev(mu1.vals),mu1.vals)

sm.pow.sat <- glm(resp.sat ~ ns(pred,40), family=binomial)
sm.pow.sel <- glm(resp.sel ~ ns(pred,40), family=binomial)


pdf("../figs/bivariateSelVSat_powCurves_0.pdf",width=4.5,height=4.3)
par(xaxs="i",mar=c(4.1,4.1,3.1,0.7))
plot(pred,predict(sm.pow.sat,type="response"),type="l",ylim=c(0,1),
     xlab=expression(mu[1]),ylab="Power",
     main=expression(paste("Power, ", mu[2]==0)))
lines(pred,predict(sm.pow.sel,type="response"),type="l",col="red")
abline(h=.05,lty=2)
abline(h=0:1,lty=2,col="gray")
legend("bottomright",lty=1,col=1:2,legend=c("Saturated","Selected"),
       bg="white")
dev.off()



###########################################
### Dependence of Saturated-Model p-Values
###########################################

# Second-step p-value for bivariate normal. Note that the saturated and selected models are the same in the second step.
p2 <- function(y1, y2) {
  1 - (1-2*pnorm(-abs(y2))) / (1-2*pnorm(-abs(y1)))
}

set.seed(1)
B <- 1E5
y <- matrix(abs(rnorm(2*B)),B,2) # Signs don't matter for the p-values
p.vals <- matrix(NA,B,2)
ymax <- apply(y,1,max)
ymin <- apply(y,1,min)
saturated.p1 <- saturated.p(ymax,ymin)
#selected.p1 <- selected.p(ymax)
both.p2 <- p2(ymax,ymin)
cor(saturated.p1, both.p2)
#cor(selected.p1, both.p2)

bin.2 <- function(p.vals) cut(p.vals,breaks=seq(0,1,.2))

rounded <- data.frame(p1 = bin.2(saturated.p1), p2 = bin.2(both.p2))
tbl <- xtabs(~p1 + p2, data = rounded)
tbl <- cbind(tbl, Total=rowSums(tbl))
tbl <- rbind(tbl, Total=colSums(tbl))
xtable(100*tbl/B, digits=1, align="l|ccccc|c")



#rounded <- data.frame(p1 = bin.2(selected.p1), p2 = bin.2(both.p2))
#tbl <- xtabs(~p1 + p2, data = rounded)
#tbl <- cbind(tbl, Total=rowSums(tbl))
#tbl <- rbind(tbl, Total=colSums(tbl))
#xtable(100*tbl/B, digits=1, align="l|ccccc|c")


image(xtabs(~rounded[,1] + rounded[,2]))
