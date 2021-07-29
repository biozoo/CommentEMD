setwd('C:\\D\\data\\Comment') # Please change the route before run this code
library(rEDM)
library(Kendall)
library(dplyr)
library(deSolve)
seed=15486
set.seed(seed)

###Function for making block
make_block <- function(data,cols,delays,lib=c(1,NROW(data))){
  lib <- matrix(lib,ncol = 2)
  data <- as.matrix(data)
  ncol <- length(cols)
  nrow <- dim(data)[1]
  block <- array(NA,dim = c(nrow,ncol))
  colnames(block) <- 1:ncol
  for (i in 1:ncol){
    I <- 1:nrow
    I_delay <- intersect(I,I+delays[i])
    block[I_delay-delays[i],i] <- data[I_delay,cols[i]]
    if (delays[i] < 0){
      # remove data points that fall at start of lib segments
      block[lib[,1] - (0:(delays[i]+1)),i] <- NA
      colnames(block)[i] <- paste(colnames(data)[cols[i]],'_t-',abs(delays[i]),sep="")
    } else if (delays[i] > 0) {
      # remove data points that fall at end of lib segments
      block[lib[,2] - (0:(delays[i]+1)),i] <- NA
      colnames(block)[i] <- paste(colnames(data)[cols[i]],'_t+',abs(delays[i]),sep="")
    } else {
      colnames(block)[i] <- paste(colnames(data)[cols[i]],'_t',sep="")
    }
  }
  return(block)
}



####### Function for generating pairwise lag series 
laf=function(x,y,lagf){
  n <- NROW(x)
  x.t=x;y.t=y
  if(lagf<=0){x.t=x.t[(1-lagf):n];y.t=y.t[1:(n+lagf)]} # if lagf<0, y is leading
  if(lagf>0){x.t=x.t[1:(n-lagf)];y.t=y.t[(1+lagf):n]}  # if lagf>0, x is leading           
  return(cbind(x.t,y.t))
}



####Fig. 1 Moran effect (CCM)
nsim=10000
eT=R1=R2=N1=N2=rep(0,nsim)
r1=3.4;r2=2.9
ph1=0.5;ph2=0.6
s1=0.4;s2=0.35
D1=3;D2=3
R1[1]=1;R2[1]=1;N1[1]=0.5;N1[1]=0.5
#set.seed(2301)
set.seed(15486)
#set.seed(3486)
for(t in 1:nsim){
  eT[t]=rnorm(1)
  R1[t+1] = N1[t]*(r1*(1-N1[t]))*exp(-ph1*eT[t])
  N1[t+1] = s1*N1[t] + max(R1[t-D1], 0)
  
  R2[t+1] = N2[t]*(r2*(1-N2[t]))*exp(-ph2*eT[t])
  N2[t+1] = s2*N2[t] + max(R2[t-D2], 0)
  
}

n=200
dam=data.frame(Time=1:n,cbind(R1,R2,N1,N2)[(nsim-n+1):nsim,])
cor(dam[,'N1'],dam[,'N2'])
write.table(dam,"Data_moran.txt",sep='\t',row.names=F)
write.table(dam[,'N1'],"N1.txt",sep='\t',row.names=F,col.names=F)
write.table(dam[,'N2'],"N2.txt",sep='\t',row.names=F,col.names=F)
write.table(dam[,'Time'],"time_moran.txt",sep='\t',row.names=F,col.names=F)

dam=read.table("Data_moran.txt",sep='\t',header=T)
n=nrow(dam)
# File path
dam.n=scale(dam[,-1], center = TRUE, scale = TRUE)

E.test.n1=NULL
for(E.t in 2:8){
  cmxy.t <- ccm(dam.n, E = E.t,
                lib_column = "N1", target_column = "N2",
                lib_sizes = n, tp=-1,random_libs = F)
  E.test.n1=rbind(E.test.n1,cmxy.t)
}
(E_n1=E.test.n1$E[which.max(E.test.n1$rho)[1]])

sim_n2=simplex(dam.n[,'N2'],lib=c(1,n),pred=c(1,n),E=c(2:8))
(E_n2=sim_n2$E[which.max(sim_n2$rho)[1]])

E.test.n2=NULL
for(E.t in 2:8){
  cmxy.t <- ccm(dam.n, E = E.t,
                lib_column = "N2", target_column = "N1",
                lib_sizes = n,tp=-1,random_libs = F)
  E.test.n2=rbind(E.test.n2,cmxy.t)
}
(E_n2=E.test.n2$E[which.max(E.test.n2$rho)[1]])

# N1 cross-map N2 (test N2 causes N1)
libs=c(seq(20,80,20),seq(100,n,50))
n1_xmap_n2 <- ccm(dam.n, E = E_n1,
                  lib_column = "N1", target_column = "N2",
                  lib_sizes = libs, 
                  num_samples = 100,replace=T,RNGseed=2301)
#plot(rho~lib_size,x_xmap_y)
# N2 cross-map N1 (test N1 causes N2)
n2_xmap_n1 <- ccm(dam.n, E = E_n2,
                  lib_column = "N2", target_column = "N1",
                  lib_sizes = libs, 
                  num_samples = 100,replace=T,RNGseed=2301)

n12q=ccm_means(n1_xmap_n2)
n21q=ccm_means(n2_xmap_n1)

# Fig. 1 Plot CCM forecast skill vs library size 
windows()
plot(n12q[,'rho']~libs,type="l",col="red",ylim=c(0,1),lwd=2,
     main="Convergent cross mapping CCM",xlab="Library size",ylab=expression(rho))

lines(n21q[,'rho']~libs,col="blue",lwd=2)

legend(0.6*n,1,c("N1 xmap N2","N2 xmap N1"),lty=c(1,1),col=c("red","blue"))
abline(h=cor(dam[,'N1'],dam[,'N2']),lty=3)



####################################################
###Fig. 2  White noise time series
# White noise generation very time consuming

if(F){
  white_noise=NULL
  for(i in 1:10000){
    white_noise=rbind(white_noise,arima.sim(model = list(order = c(0, 0, 0)), n = 100))
  }
  white_noise=write.csv('white_noise.csv',row.names=F)
}
# CCM computation (1000 times)
white_noise=read.csv('white_noise.csv',header=F)
sapn=1000;resap=NULL
for(i in 1:sapn) resap=rbind(resap,sample(10000,2))
Edim=3
tlmax=10
lib_siz=seq(2,10,2) 
ccmda=NULL
for(k in 1:sapn){
  da.t=t(white_noise[resap[k,],])
  colnames(da.t)=c('X','Y')
  da.t=apply(apply(da.t,2,as.numeric),2,scale)
  
  # a sequence of library size
  for(j in 1:(length(lib_siz)-4)){
    da.j=da.t
    da.j=da.j[1:lib_siz[length(lib_siz)],]
    x_xmap_y <- ccm(da.j, E = Edim, # The embedding dimension E for each link were determined in previous step
                    lib_column = 'X', target_column = 'Y',
                    lib_sizes = lib_siz, tp=0,RNGseed = seed,
                    random_libs = F)
    y_xmap_x <- ccm(da.j, E = Edim, # The embedding dimension E for each link were determined in previous step
                    lib_column = 'Y', target_column = 'X',
                    lib_sizes = lib_siz, tp=0,RNGseed = seed,
                    random_libs = F)
    
    # Take average for the predictive skill under each library size
    aveg.xy=cbind(unique(x_xmap_y$lib_size),aggregate(x_xmap_y[,c('rho')], by=list(as.factor(x_xmap_y$lib_size)), mean)[,'x'],
                  aggregate(x_xmap_y[,c('mae')], by=list(as.factor(x_xmap_y$lib_size)), mean)[,'x'],
                  aggregate(x_xmap_y[,c('rmse')], by=list(as.factor(x_xmap_y$lib_size)), mean)[,'x'])
    
    
    ccm_mean.xy=data.frame(x_xmap_y[1:nrow(aveg.xy),-c(1:4)]);
    ccm_mean.xy[,c('lib_size','rho','mae','rmse')]=aveg.xy
    ccm_mean.xy[ccm_mean.xy[,'rho']<0,'rho']=0
    
    aveg.yx=cbind(unique(x_xmap_y$lib_size),aggregate(x_xmap_y[,c('rho')], by=list(as.factor(x_xmap_y$lib_size)), mean)[,'x'],
                  aggregate(x_xmap_y[,c('mae')], by=list(as.factor(x_xmap_y$lib_size)), mean)[,'x'],
                  aggregate(x_xmap_y[,c('rmse')], by=list(as.factor(x_xmap_y$lib_size)), mean)[,'x'])
    
    
    ccm_mean.yx=data.frame(y_xmap_x[1:nrow(aveg.yx),-c(1:4)]);
    ccm_mean.yx[,c('lib_size','rho','mae','rmse')]=aveg.yx
    ccm_mean.yx[ccm_mean.yx[,'rho']<0,'rho']=0
    ccmda=rbind(ccmda,data.frame(w=rep(k),tl=rep(lib_siz[(j+4)]),rbind(ccm_mean.xy,ccm_mean.yx)))
  }
}

# The computation of convergence statistics
ind=unique(ccmda[,'w'])
tl=unique(ccmda[,'tl'])
ccmda.t=ccmda[ccmda[,'tl']==tlmax,]
drho=NULL
for(i in 1:length(ind)){
  ccmda.i=ccmda.t[ccmda.t[,'w']==ind[i],]
  ccmda.i[ccmda.i[,'rho']<0,'rho']=0
  ccm.i.x=ccmda.i[ccmda.i[,'lib_column']=='X',]
  ccm.i.y=ccmda.i[ccmda.i[,'lib_column']=='Y',]
  
  rho.Lmax.x=ccm.i.x[nrow(ccm.i.x),'rho']
  rho.Lmin.x=ccm.i.x[1,'rho']
  rho.Lmax.y=ccm.i.y[nrow(ccm.i.y),'rho']
  rho.Lmin.y=ccm.i.y[1,'rho']
  delta_rho.x=rho.Lmax.x-rho.Lmin.x
  delta_rho.y=rho.Lmax.y-rho.Lmin.y
  z.x=abs(0.5*(log((1+rho.Lmax.x)/(1-rho.Lmax.x))-log((1+rho.Lmin.x)/(1-rho.Lmin.x)))*(2/(tlmax-3))^-0.5)
  z.p.x=(1-pnorm(z.x))
  z.y=abs(0.5*(log((1+rho.Lmax.y)/(1-rho.Lmax.y))-log((1+rho.Lmin.y)/(1-rho.Lmin.y)))*(2/(tlmax-3))^-0.5)
  z.p.y=(1-pnorm(z.y))
  
  # Kendall's tau test
  kend.x=MannKendall(ccm.i.x[,'rho'])
  kend.tau.x=kend.x$tau[1]
  kend.p.x=kend.x$sl[[1]]
  kend.y=MannKendall(ccm.i.y[,'rho'])
  kend.tau.y=kend.y$tau[1]
  kend.p.y=kend.y$sl[[1]]
  
  drho=rbind(drho,rbind(c(ind[i],tlmax,1,delta_rho.x,z.p.x,kend.tau.x,kend.p.x),
                        c(ind[i],tlmax,2,delta_rho.x,z.p.y,kend.tau.y,kend.p.y)))
}

colnames(drho)=c('ind','tl','type','rho','p_z','tau','p_tau')
windows()
par(mfcol=c(3,1),mar=c(3,4,3,1))
drho[drho[,'rho']<0,'rho']=0
hist(drho[,'rho'],main=paste('Distribution of delta rho at L=',tlmax), xlab='delta rho',freq =F,breaks=seq(0,1,0.05),col='blue')
hist(drho[,'p_z'],main=paste('Distribution of p-value (Z test) at L=',tlmax), xlab='delta rho',freq =F,breaks=seq(0,1,0.05),col=2)
abline(v=0.05,col=2)

ptau=drho[,'p_tau'];ptau[drho[,'tau']<0]=1
hist(ptau,main=paste('Distribution of p-value (tau test) at L=',tlmax), xlab='delta rho',freq =F,breaks=seq(0,1,0.05),col=2)
abline(v=0.05,col=2)



