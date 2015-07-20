
library(lars)
library(selectiveInference,lib.loc="/Users/tibs/dropbox/git/R/mylib")


library(truncnorm)
 library(xtable)


source("exact.R")
source("hitAndRUn.R")

set.seed(635)
options(error=dump.frames)
BIG=10e9


x=read.table("data64.txt")
x=as.matrix(x)
x=scale(x,T,T)
nams=scan("data64.names",what="")
y=scan("diab.y")
y=y-mean(y)
set.seed(44)
n=nrow(x)
p=ncol(x)

maxp=p

# table for paper
a=lsfit(x,y)
s=sqrt(sum(a$res^2)/(n-p))


a=forwardStep(x,y,sigma=s)
pv=rep(NA,p)
for(j in 1:p){
    temp=lsfit(x[,a$pred[1:j]],y)
    temp2=ls.diag(temp)
    z=temp$coef/temp2$std.err
    pv[j]=2*(1-pnorm(abs(z[length(z)])))
}
      
aa=forwardStepInf(a,x,y,nsteps=p)  #saturated

d=lar(x,y,norm=F,maxsteps=25)

dd=larInf(x,y,d,sigma=s)  

#maxt
b0=forwardStepMaxtInf(a,x,y,nsteps=25,initialStep=1,trace=T,maxboot=50000,sigma=s,AbMeth="will")
save(b0,file="diab.RData")  #known sigma

b=forwardStepMaxtInf(a,x,y,nsteps=25,initialStep=1,trace=T,maxboot=50000,fix.normy=T,AbMeth="will")
save(b,file="diab.unk.RData")  #unknown sigma

#make table
load("diab.RData")
load("diab.unk.RData")

b0$pv=c(b0$pv,rep(NA,p-length(b0$pv)))
b$pv=c(b$pv,rep(NA,p-length(b$pv)))
    

tab=cbind(a$pred,round(pv,2),round(aa$pv,2),round(b$pv,2))

xtable(tab)
