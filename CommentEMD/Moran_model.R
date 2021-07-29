rm(list = ls())
seed=15486
set.seed(seed)

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

if(F){
  write.table(dam[,'N1'],"N1.txt",sep='\t',row.names=F,col.names=F)
  write.table(dam[,'N2'],"N2.txt",sep='\t',row.names=F,col.names=F)
  write.table(dam[,'Time'],"time_moran.txt",sep='\t',row.names=F,col.names=F)
}
