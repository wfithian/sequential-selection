
library(lars)
library(selectiveInference,lib.loc="/Users/tibs/dropbox/git/R/mylib")


library(truncnorm)
 library(xtable)

source("lassoinf.R")
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
res2=lassoinf(x,y,type="step",maxp=maxp)  #naive

a=forwardStep(x,y,sigma=s)
aa=forwardStepInf(a,x,y,nsteps=p)  #saturated

d=lar(x,y,norm=F,maxsteps=25)

dd=larInf(x,y,d,sigma=s)  

#maxt 
b=forwardStepMaxtInf(a,x,y,nsteps=25,initialStep=1,trace=T,maxboot=50000,fix.normy=T,AbMeth="will")

save(b,file="diab.RData")  #known sigma
save(b,file="diab.unk.RData")  #unknown sigma

#make table
load("diab.unk.RData")

b$pv=c(b$pv,rep(NA,p-length(b$pv)))
    
tab=cbind(res2$res[,1],round(res2$res[,3],4),round(aa$pv,4),round(b$pv,4))
xtable(tab)
