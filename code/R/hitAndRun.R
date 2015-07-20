


hitAndRun=function(x,y,fsfit,A,b,nsteps,sigma,nburn=1000,ntotal=2000,gap.reset=trunc(ntotal/40),ntry=NULL,return.ystar=F,fix.normy=F,trace=F){

    # not yet working for fix.normy=T!!!!
 
if(nsteps==1) stop("not yet implemented for nsteps=1")
    n=length(y)

pp=nrow(A)

mod=fsfit$pred[1:nsteps]
ystar=matrix(NA,ntotal,n)
xx=x[,mod[1:(nsteps-1)]]
mu=as.vector(y-lsfit(xx,y)$res)
newb=b-A%*%mu
ystar[1,]=y-mu

ii=1
go=T

while(ii<ntotal){
    z=rnorm(n)   #random direction
    if(ii%%gap.reset==0) z=sample(c(1,rep(0,n-1)))
    z=lsfit(xx,z,int=F)$res
    z=z/sqrt(sum(z*z))
    junk=tf.jonvs(ystar[ii,],A,newb, z)   #find truncation limits
    if(junk$vm>junk$vp){
        temp=junk$vp;junk$vp= junk$vm;junk$vm=temp  # not sure why this is needed(for fix.normy=T only)
    }
    
    if(ii%%100==0 & trace)  cat(c(ii,junk$vm,junk$vp),fill=T)
    cc=rtruncnorm(1,a=junk$vm,b=junk$vp,mean=0,sd=sigma*sqrt(sum(z*z)))
    yy=lsfit(z,ystar[ii,],int=F)$res+cc*z
 
    if(fix.normy) {
         go=F;
         syy=yy
         yy=sqrt(sum(y*y))*yy/sqrt(sum(yy*yy))
         cat(c(max(A%*%(syy+mu)-b),max(A%*%(yy+mu)-b)),fill=T)
         if (max(A%*%(yy+mu)-b)<0) {go=T}
         
     }
    if(!fix.normy | go)  {
    ystar[ii+1,]=yy;
    val=max(A%*%(ystar[ii+1,]+mu)-b)
    if(val>0) c("Warning max>0",ii+1,val,fill=T)
    ii=ii+1
}
}

ystar=scale(ystar,-mu,F)

ystar=ystar[-(1:nburn),]

 
pv=estpv(x,y,mod,nsteps,ystar)
pvless=NULL
if(!is.null(ntry)){
    pvless=rep(NA,length(ntry))
   
    for(ii in 1:length(ntry)){
           nn=nrow(ystar)-ntry[ii]
        pvless[ii]=estpv(x,y,mod,nsteps,ystar[-(1:nn),])}
}

if(!return.ystar) return(list(pv=pv,pvless=pvless,ntry=ntry))
if(return.ystar) return(list(pv=pv,ystar=ystar,pvless=pvless,ntry=ntry))
}

#from our package but modified

tf.jonvs = function(y,A,b,eta) {
    
  g = A %*% eta/sum(eta^2)
  f = b - A%*%y + g*sum(eta*y)
  temp=f/g
  vm = suppressWarnings(max((temp)[g<0]))
  vp = suppressWarnings(min((temp)[g>0]))
 # vz = suppressWarnings(min(f[g==0]))
  return(list(vm=vm,vp=vp))
}


estpv=function(x,y,mod,nsteps,ystar){
xx=x[,mod[1:(nsteps-1)]]

xr=lsfit(xx,x)$res
r=lsfit(xx,y)$res

ip=abs(t(xr)%*%r/sqrt(diag(t(xr)%*%xr)))

ip[mod[1:(nsteps-1)]]=0  #what's going on ROB
tt=max(ip)

rstar=lsfit(xx,t(ystar))$res

ipp=t(xr)%*%rstar/sqrt(diag(t(xr)%*%xr))
ipp[mod[1:(nsteps-1)],]=0


ttstar=apply(abs(ipp),2,max)
pv=mean(ttstar>tt,na.rm=T)
return(pv)
}


compvs=function(y,eta, A, b){
     if(!is.loaded("compvs")) dyn.load("/Users/tibs/dropbox/PAPERS/FourOfUs/compvs.so")
     n=length(y)
     k=nrow(A)
     mode(y)="single"
      mode(A)="single"
      mode(b)="single"
      mode(n)="integer"
      mode(k)="integer"
      mode(eta)="single"
     junk=.Fortran("compvs",
     y,
     as.matrix(A),
     b,
     n,
     k,
     eta,
     vm=single(1),
     vp=single(1),
     scr=single(k),
     scr2=single(k)
         )
     vm = suppressWarnings(max((junk$scr2)[junk$scr<0]))
  vp = suppressWarnings(min((junk$scr2)[junk$scr>0]))
     return(list(vm=vm,vp=vp))
 }
    

          

hitAndRun.yuval=
function(x,y,a,aa,nsteps,sigma,nburn=1000,ntotal=2000,return.ystar=F){
n=length(y)
A=aa$A[aa$stepind<nsteps,]
mod=a$pred
eta=rnorm(n)
xx=x[,mod[1:(nsteps-1)]]
r=lsfit(xx,y)$res
mu=as.vector(y-r)
ystar=sample_truncnorm(A, b=rep(0,nrow(A)),mu=mu, Sigma=diag(rep(sigma^2,n)),initial=y, eta=eta, burnin=nburn,ndraw=ntotal)

pv=estpv(x,y,mod,nsteps,ystar)
if(!return.ystar) out=pv
if(return.ystar) out=list(pv=pv,ystar=ystar)
return(out)
}
       

hitAndRun.boot=function(x,y,a,aa,nsteps,sigma,nburn=1000,ntotal=2000,gap.reset=trunc(ntotal/40),ntry=NULL,return.ystar=F,trace=F){
    # version that rescales resduals so var = sigma^2
if(nsteps==1) stop("not yet implemented for nsteps=1")
    n=length(y)

A=aa$A[aa$stepind<nsteps,]
pp=nrow(A)
b=rep(0,pp)
mod=a$pred[1:nsteps]
ystar=matrix(NA,ntotal,n)
xx=x[,mod[1:(nsteps-1)]]
mu=as.vector(y-lsfit(xx,y)$res)
newb=b-A%*%mu
ystar[1,]=y-mu
y0=y-mu
fac=sigma/sqrt(var(y0))
ii=1


while(ii<ntotal){
 z=rnorm(n)
     if(ii%%gap.reset==0) z=sample(c(1,rep(0,n-1)))
    z=lsfit(xx,z,int=F)$res
    z=z/sqrt(sum(z*z))
    junk=tf.jonvs(ystar[ii,],A,newb, z)
    temp=fac*sample(y0,size=n,replace=T)
    cc=sum(z*temp)
    if((junk$vm< cc) & (junk$vp> cc)){
      yy=lsfit(z,ystar[ii,],int=F)$res+cc*z
      ystar[ii+1,]=yy;
      val=max(A%*%(ystar[ii+1,]+mu)-b)
     if(val>0) c("Warning max>0",ii+1,val,fill=T)
     ii=ii+1
  }
}

ystar=scale(ystar,-mu,F)

ystar=ystar[-(1:nburn),]

 
pv=estpv(x,y,mod,nsteps,ystar)
pvless=NULL
if(!is.null(ntry)){
    pvless=rep(NA,length(ntry))
   
    for(ii in 1:length(ntry)){
           nn=nrow(ystar)-ntry[ii]
        pvless[ii]=estpv(x,y,mod,nsteps,ystar[-(1:nn),])}
}

if(!return.ystar) return(list(pv=pv,pvless=pvless,ntry=ntry))
if(return.ystar) return(list(pv=pv,ystar=ystar,pvless=pvless,ntry=ntry))
}

#from our package but modified

tf.jonvs = function(y,A,b,eta) {
    
  g = A %*% eta/sum(eta^2)
  f = b - A%*%y + g*sum(eta*y)
  temp=f/g
  vm = suppressWarnings(max((temp)[g<0]))
  vp = suppressWarnings(min((temp)[g>0]))
 # vz = suppressWarnings(min(f[g==0]))
  return(list(vm=vm,vp=vp))
}

                          
