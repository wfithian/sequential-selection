# functions for exact selective max computation




forwardStepMaxtInf=function(fsfit,x,y,nsteps,AbMeth="will",initialStep=1,sigma=1,maxboot=100000,minboot=300,nonpar=F,projres=T,fix.normy=F,scaleres=F,return.ystar=F,ntry=NULL,try.hitAndRun=TRUE,nburn=2000,ntotal=8000,trace=T){

# we do acc/ref as long we can get minboot samples in maxboot tries. o/w we switch to hit and run
 #   and use  it for the rest of the steps
#note- i don;t think idea behind  fix.normy works in general (trying to cond on ||y||)
  
    n=nrow(x)
    p=ncol(x)

    sigma.call=sigma

pv=rep(NA,nsteps)
    pvless=NULL
    if(!is.null(ntry)){
         pvless=matrix(NA,nsteps,length(ntry))
     }
sel=fsfit$pred
    
    if(AbMeth=="big"){
  fsinf=forwardStepInf(fsfit,x,y,sigma=sigma,nsteps=nsteps)
    gam=-fsinf$A
    stepind=fsinf$stepind
}
    
    xr=x;r=y;fit1=rep(0,n)
 totb=goodb=rep(0,nsteps)

    go.accrej=TRUE
    for(k in initialStep:nsteps){
        cat(k,fill=T)
        if(is.null(sigma.call)) sigma=sqrt(sum(y^2)/n)
        if(k>1){ sel1=sel[1:(k-1)]
        xx=x[,sel1]
                 proj=xx%*%solve(t(xx)%*%xx)%*%t(xx)
            bhat1=lsfit(x[,sel1,drop=F],y)$coef[-1]
                 fit1=x[,sel1,drop=F]%*%bhat1  
     
             xr=lsfit(xx,x)$res
            r=lsfit(xx,y)$res
                  if(is.null(sigma.call)) sigma=sqrt(sum(r^2)/(n-k))
             }
ip=t(xr)%*%r/sqrt(diag(t(xr)%*%xr))
  if(k>1)  ip[sel1]=0
    tt=max(abs(ip))
if(k>1){
         if(AbMeth=="big"){
       gam0=gam[stepind<=(k-1),]
        u0=rep(0,nrow(gam0))
        stepind0=stepind[stepind<=(k-1)]
         }
         if(AbMeth=="will"){
           willAb=compute.willAb(x,y,sel1)
           gam0=-willAb$A;u0=-willAb$b
       }
     }
ystarall=matrix(NA,minboot,n)
sy=sum(y*y)
 facres=1      
 if(scaleres) facres= sigma/sqrt(var( (y-fit1)))
      
 if(k==1) {
        if(!nonpar) {e= sigma*matrix(rnorm(minboot*n),ncol=n);fac=1}
         if(nonpar){ fac=sqrt(n/(n-k-1))
                     e=matrix(sample(r,size=n*minboot,replace=T),ncol=n)
                 }
         y1=matrix(fit1,nrow=minboot,ncol=n,byrow=T)+fac*e
          goodb[1]=totb[1]=minboot
      }
if(k>1) {
     y1=NULL
     if(go.accrej)  {
          go.accrej=FALSE
         junk=genystar5(gam0,u0,k,fit1,diag(n)-proj,r,sigma,sy,facres,minboot=minboot,maxboot=maxboot,nonpar=nonpar,projres=projres,fix.normy=fix.normy)
         y1=junk$ystar;goodb[k]=junk$goodb;totb[k]=junk$totb
                      if(goodb[k]==minboot) go.accrej=TRUE
                  }
        if(try.hitAndRun & goodb[k]<minboot){
            cat(c("step=",k,"; using hit and run"),fill=T)
            if(AbMeth=="big"){
            A=fsinf$A[fsinf$stepind<nsteps,]
            b=rep(0,nrow(A))
        }
            if(AbMeth=="will") {A=willAb$A;b=willAb$b}

            if(!fix.normy){
                junk=hitAndRun(x,y,fsfit,A,b,sigma=sigma,nsteps,nburn=nburn,ntotal=ntotal,return.ystar=T,
                    trace=trace,fix.normy=fix.normy)
                 ystar=junk$ystar
            }

            if(fix.normy) {
                cat("hitAndRun on sphere",fill=T)
                mu=as.vector(y-lsfit(xx,y)$res)
              junk=conditional.cond(A,b,t(xx),t(xx)%*%y) 
           con=list(A=A,b=b,mu=junk$mu,Sigma=junk$Sigma)
            out=sample_from_sphere(con,y,burnin=nburn,ndraw=ntotal,white=F,trace=F)
         ystar=out$Z
            }
           
           y1=rbind(y1,ystar)
            
       
        }

        }
        #note: i center y1 below, since original y was centered
        # better to call estpv  fun below
        y1=t(scale(t(y1),T,F))
        rstar=t(y1)
      if(k>1) rstar=(diag(n)-proj)%*%t(y1)
       ipp=t(xr)%*%rstar/sqrt(diag(t(xr)%*%xr))
       if(k>1) for(kk in 1:(k-1)){
          ipp[sel1[kk],]=0
      }
      
        ttstar=apply(abs(ipp),2,max)
     
        
pv[k]=mean(ttstar>tt)
     if(!is.null(ntry)){
         for(ii in 1:length(ntry)){
              pvless[k,ii]=mean(ttstar[1:ntry[ii]]>tt)
    }}

    }
  out2=NULL
    if(return.ystar) out2=y1
 return(list(pv=pv,totb=totb,goodb=goodb,sel=sel,pvless=pvless,ntry=ntry,facres=facres,ystar=out2))

}






maxselpv=function(x,y,nsteps,sigma=1,maxboot=50000,minboot=500, projres=F,fix.normy=F,return.ttstar=F,ntry=NULL){
    # R version
a=myfs(x,y,stand=F,nsteps=nsteps)
pv=rep(NA,nsteps)
sel=a$pred
    junk=myfs.pval(a,x,y,sigma,nsteps=nsteps)
gam=-junk$A
    stepind=junk$stepind
    xr=x;r=y;fit1=rep(0,n)
 totb=goodb=rep(0,nsteps)
proj=matrix(0,nrow=n,ncol=n)
    for(k in 1:nsteps){
        cat(k)
        if(k>1){ sel1=sel[1:(k-1)]
        xx=x[,sel1]
        proj=xx%*%solve(t(xx)%*%xx)%*%t(xx)
      fit1=proj%*%y
     xr=lsfit(x[,sel1,drop=F],x)$res
     r=y-fit1
             }
ip=t(xr)%*%r/sqrt(diag(t(xr)%*%xr))
  if(k>1)  ip[sel1]=0
    tt=max(abs(ip))
if(k>1) gam0=gam[stepind<=(k-1),]
ystarall=matrix(NA,minboot,n)
bb=0
 b=0
while(bb<minboot & b<maxboot){
b=b+1
if(!projres)  ystar=fit1+sigma*rnorm(n)
if(projres) {e=sigma*rnorm(n);ystar=fit1+(diag(n)-proj)%*%e}
nnn=T
if(fix.normy) {
         u2=(diag(n)-proj)%*%rnorm(n)
             fac=sqrt((sum(y^2)-sum(fit1^2))/sum(u2^2))
             ystar=fit1+fac*u2
     }
  
if(k>1) nnn= min(gam0%*%ystar)>0
      if(nnn){bb=bb+1
      ystarall[bb,]=ystar
    
    }
}
goodb[k]=bb
 totb[k]=b
        
if(goodb[k]<minboot) cat(c("step=",k,": need more boots"),fill=T)
y1=ystarall
        #note: i center y1 below, since original y was centered
        y1=t(scale(t(y1),T,F))
        rstar=t(y1)
if(k>1) rstar=(diag(n)-proj)%*%t(y1)
ipp=t(xr)%*%rstar/sqrt(diag(t(xr)%*%xr))
   if(k>1) for(kk in 1:(k-1)){
          ipp[sel1[kk],]=0
      }
ttstar=apply(abs(ipp),2,max)
pv[k]=mean(ttstar>tt,na.rm=T)
    }
if(!return.ttstar){return(list(pv=pv,totb=totb,goodb=goodb,sel=sel))}
if(return.ttstar){return(list(pv=pv,totb=totb,goodb=goodb,sel=sel,tt=tt,ttstar=ttstar))}
}




maxselpv22=function(x,y,nsteps,sigma=1,maxboot=50000,minboot=500,ord=F){

    #fortran version with extra screen
a=myfs(x,y,stand=F,nsteps=nsteps)
pv=rep(NA,nsteps)
sel=a$pred
    b=myfs.pval(a,x,y,sigma,nsteps=nsteps)
gam=-b$A
    stepind=b$stepind
    xr=x;r=y;fit1=rep(0,n)
 totb=goodb=rep(0,nsteps)

    for(k in 1:nsteps){
        cat(k)
        if(k>1){ sel1=sel[1:(k-1)]
        xx=x[,sel1]
        proj=xx%*%solve(t(xx)%*%xx)%*%t(xx)
      fit1=proj%*%y
     xr=lsfit(x[,sel1,drop=F],x)$res
     r=y-fit1
             }
ip=t(xr)%*%r/sqrt(diag(t(xr)%*%xr))
  if(k>1)  ip[sel1]=0
    tt=max(abs(ip))
if(k>1) {gam0=gam[stepind<=(k-1),]
         stepind0=stepind[stepind<=(k-1)]
     }
ystarall=matrix(NA,minboot,n)


 if(k==1) {
          y1=matrix(fit1,nrow=minboot,ncol=n,byrow=T)+sigma*matrix(rnorm(minboot*n),ncol=n)
          goodb[1]=totb[1]=minboot
      }
if(k>1) {
        junk=genystar2(gam0,stepind0,fit1,sigma,minboot=minboot,maxboot=maxboot,ord=ord,y=y)
         y1=junk$ystar;goodb[k]=junk$goodb;totb[k]=junk$totb
        
        }
        #note: i center y1 below, since original y was centered
        y1=t(scale(t(y1),T,F))
        rstar=t(y1)
if(k>1) rstar=(diag(n)-proj)%*%t(y1)
ipp=t(xr)%*%rstar/sqrt(diag(t(xr)%*%xr))
   if(k>1) for(kk in 1:(k-1)){
          ipp[sel1[kk],]=0
      }
ttstar=apply(abs(ipp),2,max)
pv[k]=mean(ttstar>tt)
     
    }
return(list(pv=pv,totb=totb,goodb=goodb,sel=sel))
}



comp.beta.IS=function(x,y,fsfit,nsteps){
    # compute bhat's for imp sampling- increase small coefs to sqrt(2*log(p))
    p=ncol(x)
b=lsfit(x[,fsfit$pred[1:nsteps]],y)
bhat=b$coef[-1]
bb=ls.diag(b)
se=bb$std.err[-1]
z=(bhat/se)
o=abs(z)<sqrt(2*log(p-nsteps))
del=rep(0,nsteps)
del[o]=sign(bhat[o])*se[o]*sqrt(2*log(p-nsteps))-bhat[o]

bhat2=bhat+del
    return(bhat2)
}


  genystar5=function(gam,u,nvar,fit,resproj,res,sigma,sy,facres,seed=1234,minboot=200,maxboot=5000,nonpar=F,
      projres=F,fix.normy=F,scaleres=F){
        #generate samples y* with Gam%*%y* ge 0
      #nvar is # of predictors entered (for adjusting estimate of sigma^2)
       n=ncol(gam)
       kk=nrow(gam)
      mode(n)="integer"
       mode(kk)="integer"
         mode(nvar)="integer"
     mode(minboot)="integer"
       mode(maxboot)="integer"
      mode(gam)="single"
       mode(u)="single"
        mode(sigma)="single"
        mode(sy)="single"
         mode(facres)="single"
        mode(fit)="single"
        mode(resproj)="single"
       mode(seed)="integer"
nonpar=(1*nonpar)
      mode(nonpar)="integer"
       projres=1*projres
         mode(projres)="integer"
       fixnormy=1*(fix.normy)
       mode(fixnormy)="integer"
#scale up residuals to make var unbiassed
       res=res*sqrt(n/(n-nvar-1))
 mode(res)="single"
     
 
 if(!is.loaded("genystar5")) dyn.load("/Users/tibs/dropbox/PAPERS/FourOfUs/genystar5.so")
  
     junk=.Fortran("genystar5",
     n,
     kk,
     t(gam),
     u,
     fit,
     resproj,
     res,
     sigma,
     minboot,
         maxboot,
          nonpar,
           projres,
           fixnormy,
           sy,
           facres,
           seed,
          goodb=integer(1),
          totb=integer(1),
       ystar=single(n*minboot),
      scr=single(n),
       scr2=single(kk),
         scr3=single(n)
)
 
ystar=matrix(junk$ystar,nrow=minboot,ncol=n,byrow=T)
       return(list(ystar=ystar,goodb=junk$goodb,totb=junk$totb))
   }

compute.willAb=function(x,y,mod){
#use will's trick to compute A,b for FS
#mod is  current model
p=ncol(x)
k=length(mod)

cc=jhat=norms=NULL

# compute c and norms from page 18
norms0=sqrt(colSums(x^2))
cc0=max(abs(t(x)%*%y)/norms0)
jhat0=which.max(abs(t(x)%*%y)/norms0)


if(k>1){
    cc=jhat=rep(NA,k-1)
    norms=vector("list",k-1)
for(i in 1:(k-1)){

 mod2=mod[1:i]
 junk=lsfit(x[,mod2,drop=F],y)
 junk2=lsfit(x[,mod2,drop=F],x)
 norms[[i]]=sqrt(colSums(junk2$res^2))

 fit=y-junk$res
 ip=abs(t(x)%*%junk$res/norms[[i]])
 ip[mod2]=0
   cc[i]=max(ip)
 jhat[i]=which.max(ip)
}}

# now compute vm and vp from 35,36

modall=mod[1:k]
vmm=vpp=matrix(NA,p,k-1)
vm0=vp0=rep(NA,p)



for(j in (1:p)[-modall]){
    vm0[j]=-cc0*norms0[j]
    vp0[j]=+cc0*norms0[j]
 if(k>1){
 for(i in 1:(k-1)){
  mod2=mod[1:i]
 fit=y-lsfit(x[,mod2,drop=F],y)$res
 vmm[j,i]=sum(x[,j]*fit)-cc[i]*norms[[i]][j]
 vpp[j,i]=sum(x[,j]*fit)+cc[i]*norms[[i]][j]

}}}


vm=apply(cbind(vm0,vmm),1,max)
vp=apply(cbind(vp0,vpp),1,min)

xx=x[,-modall,drop=F]

A=rbind(-t(xx),t(xx))
b=c(-vm[!is.na(vm)], vp[!is.na(vp)])

# upper and lower bounds for p-val approx
L=rep(-Inf,p)
U=rep(Inf,p)
if(k>1){
for(j in (1:p)[-modall]){
  mod2=mod[1:(k-1)]
 fit=y-lsfit(x[,mod2,drop=F],y)$res
  L[j]=(vm[j]-sum(x[,j]*fit))/norms[[k-1]][j]
    U[j]=(vp[j]-sum(x[,j]*fit))/norms[[k-1]][j]
      }
}
 return(list(A=A,b=b,cc0=cc0,cc=cc,norms0=norms0,norms=norms,jhat0=jhat0,jhat=jhat,vm=vm,vp=vp,L=L,U=U))
}



exactLassoStepInf=function(larfit,x,y,nsteps,initialStep=1,sigma=1,maxboot=100000,minboot=300,nonpar=F,projres=T,fix.normy=F,scaleres=F,return.more=F,ntry=NULL,try.hitAndRun=TRUE,nburn=2000,ntotal=8000,trace=T){

# we do acc/ref as long we can get minboot samples in maxboot tries. o/w we switch to hit and run
 #   and use  it for the rest of the steps
#note- i don;t think idea behind  fix.normy works in general (trying to cond on ||y||)
  
    n=nrow(x)
    p=ncol(x)
sel=unlist(lar$act)
    sigma.call=sigma

pv=rep(NA,nsteps)
 
    
  
    
    xr=x;r=y;fit1=rep(0,n)
 totb=goodb=rep(0,nsteps)

    go.accrej=TRUE
    for(k in initialStep:nsteps){
        cat(k,fill=T)
        if(is.null(sigma.call)) sigma=sqrt(sum(y^2)/n)
        if(k>1){ sel1=sel[1:(k-1)]
        xx=x[,sel1]
                 proj=xx%*%solve(t(xx)%*%xx)%*%t(xx)
            bhat1=lsfit(x[,sel1,drop=F],y)$coef[-1]
                 fit1=x[,sel1,drop=F]%*%bhat1  
     
             xr=lsfit(xx,x)$res
            r=lsfit(xx,y)$res
                  if(is.null(sigma.call)) sigma=sqrt(sum(r^2)/(n-k))
             }
ip=t(xr)%*%r/sqrt(diag(t(xr)%*%xr))
  if(k>1)  ip[sel1]=0
    tt=max(abs(ip))
if(k>1){
         if(AbMeth=="big"){
       gam0=gam[stepind<=(k-1),]
        u0=rep(0,nrow(gam0))
        stepind0=stepind[stepind<=(k-1)]
         }
         if(AbMeth=="will"){
           willAb=compute.willAb(x,y,sel1)
           gam0=-willAb$A;u0=-willAb$b
       }
     }
ystarall=matrix(NA,minboot,n)
sy=sum(y*y)
 facres=1      
 if(scaleres) facres= sigma/sqrt(var( (y-fit1)))
      
 if(k==1) {
        if(!nonpar) {e= sigma*matrix(rnorm(minboot*n),ncol=n);fac=1}
         if(nonpar){ fac=sqrt(n/(n-k-1))
                     e=matrix(sample(r,size=n*minboot,replace=T),ncol=n)
                 }
         y1=matrix(fit1,nrow=minboot,ncol=n,byrow=T)+fac*e
          goodb[1]=totb[1]=minboot
      }
if(k>1) {
     y1=NULL
     if(go.accrej)  {
          go.accrej=FALSE
         junk=genystar5(gam0,u0,k,fit1,diag(n)-proj,r,sigma,sy,facres,minboot=minboot,maxboot=maxboot,nonpar=nonpar,projres=projres,fix.normy=fix.normy)
         y1=junk$ystar;goodb[k]=junk$goodb;totb[k]=junk$totb
                      if(goodb[k]==minboot) go.accrej=TRUE
                  }
        if(try.hitAndRun & goodb[k]<minboot){
            cat(c("step=",k,"; using hit and run"),fill=T)
            if(AbMeth=="big"){
            A=fsinf$A[fsinf$stepind<nsteps,]
            b=rep(0,nrow(A))
        }
            if(AbMeth=="will") {A=willAb$A;b=willAb$b}
            junk=hitAndRun(x,y,fsfit,A,b,sigma=sigma,nsteps,nburn=nburn,ntotal=ntotal,return.ystar=T,trace=trace,fix.normy=fix.normy)
           # ystar=sample_truncnorm(fsinf$A, fsinf$b, y, rnorm(n),
           #                         Sigma = diag(ncol(fsinf$A)),
           #                         mu = rep(0,ncol(fsinf$A)),
           #                         how_often=1000,
            #                        burnin=2000,
           #                         ndraw=8000,
           #                         thinning=1,
           #                         use_A=TRUE)
            #  y1=rbind(y1,ystar)
           y1=rbind(y1,junk$ystar)
            
       
        }

        }
        #note: i center y1 below, since original y was centered
        y1=t(scale(t(y1),T,F))
        rstar=t(y1)
      if(k>1) rstar=(diag(n)-proj)%*%t(y1)
       ipp=t(xr)%*%rstar/sqrt(diag(t(xr)%*%xr))
       if(k>1) for(kk in 1:(k-1)){
          ipp[sel1[kk],]=0
      }
      
        ttstar=apply(abs(ipp),2,max)
     
        
pv[k]=mean(ttstar>tt)
     if(!is.null(ntry)){
         for(ii in 1:length(ntry)){
              pvless[k,ii]=mean(ttstar[1:ntry[ii]]>tt)
    }}

    }
    
 return(list(pv=pv,totb=totb,goodb=goodb,sel=sel,pvless=pvless,ntry=ntry,facres=facres))

}


foolam=function(x,y,larfit,nsteps,sigma,initialStep=1,minboot=300,maxboot=500000,nonpar=F,fix.normy=F){
    # selective pvalues for lar   using next-lam statistic
lam=larfit$lam
nk=larfit$nk
sel=unlist(larfit$act)
sy=sum(y*y)
    xr=x;r=y;fit1=rep(0,n)
pv= totb=goodb=rep(0,nsteps)
nextlamstar=rep(NA,minboot)



 
    for(k in initialStep:nsteps){
        cat(k,fill=T)
        
        if(k>1){ sel1=sel[1:(k-1)]
        xx=x[,sel1]
                 proj=xx%*%solve(t(xx)%*%xx)%*%t(xx)
            bhat1=lsfit(x[,sel1,drop=F],y)$coef[-1]
                 fit1=x[,sel1,drop=F]%*%bhat1  
     
             xr=lsfit(xx,x)$res
            r=lsfit(xx,y)$res
     
             }
        
          if(k==1) nextlam=lam[1]
         if(k>1)  nextlam=   lam[k]-lam[k-1]

         if(k>1){
           gam0=larfit$Gamma[1:nk[k-1],]
           u0=rep(0,nrow(gam0))
           }
          ystarall=matrix(NA,minboot,n)
         sy=sum(y*y)
            facres=1
      
          if(k==1) {
          e= sigma*matrix(rnorm(minboot*n),ncol=n);fac=1
           y1=matrix(fit1,nrow=minboot,ncol=n,byrow=T)+fac*e
           goodb[1]=totb[1]=minboot
             }
        if(k>1) {
  
   
         junk=genystar5(gam0,u0,k,fit1,diag(n)-proj,r,sigma,sy,facres,minboot=minboot,maxboot=maxboot,nonpar=nonpar,projres=T,fix.normy=fix.normy)
         y1=junk$ystar;goodb[k]=junk$goodb;totb[k]=junk$totb
       
    
        #note: i center y1 below, since original y was centered
        y1=t(scale(t(y1),T,F))
 }
         for(iii in 1:nrow(y1)){
           #   cat(iii)
            aaa=lar(x,y1[iii,],normalize=F,maxsteps=k)
            if(k==1) nextlamstar[iii]=aaa$lam[1]
            if(k>1)  nextlamstar[iii]=aaa$lam[k]-aaa$lam[k-1]
        }
       
pv[k]=mean(nextlamstar>nextlam)
     

    }
    
 return(list(pv=pv,totb=totb,goodb=goodb,sel=sel))

}


