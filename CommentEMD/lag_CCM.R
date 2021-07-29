#####################################################################################################3
# R code for lagged CCM analysis
# 2021-07-29, Chun-Wei Chang
# Empirical time series are open-access and downloaded from the websites offered in Yang et al. (2019) Nature Communications
rm(list = ls())
if(F){ #Execute the code if you have not installed the packages
   install.packages("devtools")
   devtools::install_github("SugiharaLab/rEDM")
   install.packages("deSolve")
   install.packages("dplyr") 
}

library(rEDM)
library(deSolve)
library(dplyr)
setwd('D:\\data\\Comment\\causal_EMD\\CommentEMD\\Data_files') # Path to data files

cri='rho' # Criteria for parameter selection
lag.s=seq(-6,0) #
# Function for making lagged time series
make_lag=function(x,y,lag.seq,name=NULL){
   if(length(x)!=length(y)){print('unequal length of x and y')}else{
      n=length(x)
      lagda=list()
      for(i in 1:length(lag.seq)){
         lag.t=lag.seq[i];x.lag=x;y.lag=y
         if(lag.t<0){
            x.lag=x.lag[-c(1:abs(lag.t))];
            y.lag=y.lag[-c((n-abs(lag.t)+1):n)]
         }else if(lag.t>0){
            x.lag=x.lag[-c((n-lag.t+1):n)];
            y.lag=y.lag[-c(1:lag.t)]
         }
         da.t=cbind(x.lag,y.lag)
         if(!is.null(name)){colnames(da.t)=name}
         lagda[[i]]= da.t
      }# end of i
      names(lagda)=lag.seq
      return(lagda)
   }
}

##################################################
# LV predator-prey model
##################################################
parameters <- c(alpha = 1, bet = 0.05, delta = 0.02, gamm=0.5 )
state <- c(x=0.5, y=0.5)
LV<-function(t, state, parameters) {
   with(as.list(c(state, parameters)),{
     # rate of change
       dx <- alpha*x - bet*x*y
       dy <- delta*x*y - gamm*y
       # return the rate of change
        list(c(dx, dy))
      }) 
   }
times <- seq(0, 1140, by = 0.01)

# ODE solution  
sub.time <- seq(1001,1140,1)
out <- ode(y = state, times = times, func = LV, parms = parameters)
out <- as.data.frame(out)
out.2 <- filter(out,time%in%sub.time)
out.s <- out.2;out.s[,-1]=apply(out.s[,-1],2,scale)


#######################################################################
#### Lag CCM for LV model
#######################################################################
x.s <- out.2;x.s[,-1]=apply(x.s[,-1],2,scale)

# Select embedding dimensions based on hindcast ability
cm.x=cm.y=NULL
for(Ed in 2:8){
   cm.x=rbind(cm.x,ccm(x.s, E = Ed, random_libs = F,
                       lib_column = "x", target_column = "y", 
                       lib_sizes = nrow(x.s), tp=-1))
   cm.y=rbind(cm.y,ccm(x.s, E = Ed, random_libs = F,
                       lib_column = "y", target_column = "x", 
                       lib_sizes = nrow(x.s), tp=-1))
}
(Eb1=cm.x[which.max(cm.x[,cri]),'E'])
(Eb2=cm.y[which.max(cm.y[,cri]),'E'])

lag.s=seq(lag.max,0)
lagdata_x=make_lag(x=x.s[,'x'],y=x.s[,'y'],lag.seq=lag.s,name=c('x','y'))
lagdata_y=make_lag(x=x.s[,'y'],y=x.s[,'x'],lag.seq=lag.s,name=c('y','x'))

# Select the optimal lag 
rho.terminal=NULL
for(i in 1:length(lag.s)){
   x.lag=lagdata_x[[i]]
   y.lag=lagdata_y[[i]]
   lsize.x=nrow(x.lag);lsize.y=nrow(y.lag)
   x_xmap_y <- ccm(x.lag, E = Eb1, 
                   lib_column = "x", target_column = "y", 
                   lib_sizes = lsize.x, num_samples = 100,RNGseed = 17568)
   y_xmap_x <- ccm(y.lag, E = Eb2, 
                   lib_column = "y", target_column = "x", 
                   lib_sizes = lsize.y, num_samples = 100,RNGseed = 17568)
   
   xxmapy <- ccm_means(x_xmap_y)
   yxmapx <- ccm_means(y_xmap_x)
   rho.terminal=rbind(rho.terminal,
                      c(lag=lag.s[i],
                        rho_x_term=xxmapy[nrow(xxmapy),cri],
                        rho_y_term=yxmapx[nrow(yxmapx),cri]
                      ))
}

(lag_xxmapy=rho.terminal[which.max(rho.terminal[,'rho_x_term']),'lag'])
(lag_yxmapx=rho.terminal[which.max(rho.terminal[,'rho_y_term']),'lag'])

# time series with optimal lags
nolagda=lagdata_x[[which(lag.s==0)]]
lagda_x_xmap_y=lagdata_x[[which(lag.s==lag_xxmapy)]]
lagda_y_xmap_x=lagdata_y[[which(lag.s==lag_yxmapx)]]


win.graph(600,200);par(mfcol=c(1,3))
# P1 No lag
lsize=unique(c(seq(5, nrow(nolagda), by = 2),nrow(nolagda)))
x_xmap_y <- ccm(nolagda, E = Eb1, 
                lib_column = "x", target_column = "y", 
                lib_sizes = lsize, num_samples = 100,RNGseed = 17568)
y_xmap_x <- ccm(nolagda, E = Eb2, 
                lib_column = "y", target_column = "x", 
                lib_sizes = lsize, num_samples = 100,RNGseed = 17568)

xxmapy <- ccm_means(x_xmap_y)
yxmapx <- ccm_means(y_xmap_x)

plot(rho~lib_size,xxmapy,
     xlab='Library size (L)',ylab='Correlation coefficient (rho)',
     main='LV model system (no lag)',col=4,
     type='l',lwd=2,ylim=c(min(c(xxmapy[,cri],yxmapx[,cri])),max(rho.terminal[,-1]))
     )
lines(rho~lib_size,yxmapx,col=2,lwd=2)
legend('bottomright',
       c(paste('x_xmap_y (top-down;lag=',0,')',sep=''),
         paste('y_xmap_x (bottom-up;lag=',0,')',sep='')),
       col=c(4,2),lty=c(1,1),bty='n')

# P2 Lags vs skill
plot(rho_x_term~lag,rho.terminal,type='l',ylim=range(rho.terminal[,-1]),lwd=2,ylab='Correlation coefficient (rho)',col=4)
lines(rho_y_term~lag,rho.terminal,type='l',col=2,lwd=2)
abline(v=c(lag_xxmapy),lty=2,col=4);abline(v=c(lag_yxmapx),lty=2,col=2)
legend('bottomright',c('x_xmap_y','y_xmap_x'),
       col=c(4,2),lty=c(1,1),bty='n')

# P3 Lags
lsize=unique(c(seq(5, nrow(lagda_x_xmap_y), by = 2),nrow(lagda_x_xmap_y)))
x_xmap_y <- ccm(lagda_x_xmap_y, E = Eb1, 
                lib_column = "x", target_column = "y", 
                lib_sizes = lsize, num_samples = 100,RNGseed = 17568)
lsize=unique(c(seq(5, nrow(lagda_y_xmap_x), by = 2),nrow(lagda_y_xmap_x)))
y_xmap_x <- ccm(lagda_y_xmap_x, E = Eb2, 
                lib_column = "y", target_column = "x", 
                lib_sizes = lsize, num_samples = 100,RNGseed = 17568)

xxmapy <- ccm_means(x_xmap_y)
yxmapx <- ccm_means(y_xmap_x)

plot(rho~lib_size,xxmapy,
     xlab='Library size (L)',ylab='Correlation coefficient (rho)',
     main='LV model system (lagged)',col=4,
     type='l',lwd=2,ylim=c(min(c(xxmapy[,cri],yxmapx[,cri])),max(rho.terminal[,-1]))
)
lines(rho~lib_size,yxmapx,col=2,lwd=2)
legend('bottomright',
       c(paste('x_xmap_y (top-down;lag=',lag_xxmapy,')',sep=''),
         paste('y_xmap_x (bottom-up;lag=',lag_yxmapx,')',sep='')),
       col=c(4,2),lty=c(1,1),bty='n')


#######################################################################
#### Paramecium-Didinium interaction
Paramecium=read.csv('para_didi_0500_1.csv',header=T)
Paramecium.s <- Paramecium;Paramecium.s[,-1]=apply(Paramecium.s[,-1],2,scale)

# Select embedding dimensions based on hindcast ability
cm.x=cm.y=NULL
for(Ed in 2:8){
   cm.x=rbind(cm.x,ccm(Paramecium.s, E = Ed, random_libs = F,
             lib_column = "Paramecium", target_column = "Didinium", 
             lib_sizes = nrow(Paramecium.s), tp=-1))
   cm.y=rbind(cm.y,ccm(Paramecium.s, E = Ed, random_libs = F,
             lib_column = "Didinium", target_column = "Paramecium", 
             lib_sizes = nrow(Paramecium.s), tp=-1))
}
(Eb1=cm.x[which.max(cm.x[,cri]),'E'])
(Eb2=cm.y[which.max(cm.y[,cri]),'E'])

# Select the optimal lag 
lagdata_x=make_lag(x=Paramecium.s[,'Paramecium'],y=Paramecium.s[,'Didinium'],lag.seq=lag.s,name=c('Paramecium','Didinium'))
lagdata_y=make_lag(x=Paramecium.s[,'Didinium'],y=Paramecium.s[,'Paramecium'],lag.seq=lag.s,name=c('Didinium','Paramecium'))
rho.terminal=NULL
for(i in 1:length(lag.s)){
   x.lag=lagdata_x[[i]]
   y.lag=lagdata_y[[i]]
   lsize.x=nrow(x.lag);lsize.y=nrow(y.lag)
   x_xmap_y <- ccm(x.lag, E = Eb1, 
                   lib_column = "Paramecium", target_column = "Didinium", 
                   lib_sizes = lsize.x, num_samples = 100,RNGseed = 17568)
   y_xmap_x <- ccm(y.lag, E = Eb2,                   
                   lib_column = "Didinium", target_column = "Paramecium", 
                   lib_sizes = lsize.y, num_samples = 100,RNGseed = 17568)
   
   xxmapy <- ccm_means(x_xmap_y)
   yxmapx <- ccm_means(y_xmap_x)
   rho.terminal=rbind(rho.terminal,
                      c(lag=lag.s[i],
                        rho_x_term=xxmapy[nrow(xxmapy),cri],
                        rho_y_term=yxmapx[nrow(yxmapx),cri]
                      ))
}

(lag_xxmapy=rho.terminal[which.max(rho.terminal[,'rho_x_term']),'lag'])
(lag_yxmapx=rho.terminal[which.max(rho.terminal[,'rho_y_term']),'lag'])

# time series with optimal lags
nolagda=lagdata_x[[which(lag.s==0)]]
lagda_x_xmap_y=lagdata_x[[which(lag.s==lag_xxmapy)]]
lagda_y_xmap_x=lagdata_y[[which(lag.s==lag_yxmapx)]]


win.graph(600,200);par(mfcol=c(1,3))
# P1 No lag
lsize=unique(c(seq(5, nrow(nolagda), by = 2),nrow(nolagda)))
x_xmap_y <- ccm(nolagda, E = Eb1, 
                lib_column = "Paramecium", target_column = "Didinium", 
                lib_sizes = lsize, num_samples = 100,RNGseed = 17568)
y_xmap_x <- ccm(nolagda, E = Eb2, 
                lib_column = "Didinium", target_column = "Paramecium", 
                lib_sizes = lsize, num_samples = 100,RNGseed = 17568)

xxmapy <- ccm_means(x_xmap_y)
yxmapx <- ccm_means(y_xmap_x)

plot(rho~lib_size,xxmapy,
     xlab='Library size (L)',ylab='Correlation coefficient (rho)',
     main='Didinium-Paramecium system (no lag)',col=4,
     type='l',lwd=2,ylim=c(min(c(xxmapy[,cri],yxmapx[,cri])),0.9)
)
lines(rho~lib_size,yxmapx,col=2,lwd=2)
legend('bottomright',
       c(paste('Paramecium_xmap_Didinium (top-down;lag=',0,')',sep=''),
         paste('Didinium_xmap_Paramecium (bottom-up;lag=',0,')',sep='')),
       col=1:2,lty=c(1,1),bty='n')


# P2 Lags vs skill
plot(rho_x_term~lag,rho.terminal,type='l',ylim=range(rho.terminal[,-1]),lwd=2,ylab='Correlation coefficient (rho)',col=4)
lines(rho_y_term~lag,rho.terminal,type='l',col=2,lwd=2)
abline(v=c(lag_xxmapy),lty=2,col=4);abline(v=c(lag_yxmapx),lty=2,col=2)
legend('bottomright',c('Paramecium_xmap_Didinium','Didinium_xmap_Paramecium'),
       col=c(4,2),lty=c(1,1),bty='n')


# P3 Lags
lsize=unique(c(seq(5, nrow(lagda_x_xmap_y), by = 2),nrow(lagda_x_xmap_y)))
x_xmap_y <- ccm(lagda_x_xmap_y, E = Eb1, 
                lib_column = "Paramecium", target_column = "Didinium", 
                lib_sizes = lsize, num_samples = 100,RNGseed = 17568)
lsize=unique(c(seq(5, nrow(lagda_y_xmap_x), by = 2),nrow(lagda_y_xmap_x)))
y_xmap_x <- ccm(lagda_y_xmap_x, E = Eb2, 
                lib_column = "Didinium", target_column = "Paramecium", 
                lib_sizes = lsize, num_samples = 100,RNGseed = 17568)

xxmapy <- ccm_means(x_xmap_y)
yxmapx <- ccm_means(y_xmap_x)

plot(rho~lib_size,xxmapy,
     xlab='Library size (L)',ylab='Correlation coefficient (rho)',
     main='Didinium-Paramecium system (lagged)',col=4,
     type='l',lwd=2,ylim=c(min(c(xxmapy[,cri],yxmapx[,cri])),0.9)
)
lines(rho~lib_size,yxmapx,col=2,lwd=2)
legend('bottomright',
       c(paste('Paramecium_xmap_Didinium (top-down;lag=',lag_xxmapy,')',sep=''),
         paste('Didinium_xmap_Paramecium (bottom-up;lag=',lag_yxmapx,')',sep='')),
       col=c(4,2),lty=c(1,1),bty='n')




#######################################################################3
#### Moose-Wolf interaction
moose=read.csv('Moose.csv',header=T)
moose.s <- moose;moose.s[,-1]=apply(moose.s[,-1],2,scale)

# Select embedding dimensions based on hindcast ability
cm.x=cm.y=NULL
for(Ed in 2:8){
   cm.x=rbind(cm.x,ccm(moose.s, E = Ed, random_libs = F,
                       lib_column = "Moose", target_column = "Wolves",
                       lib_sizes = nrow(moose.s), tp=-1))
   cm.y=rbind(cm.y,ccm(moose.s, E = Ed, random_libs = F,
                       lib_column = "Wolves", target_column = "Moose", 
                       lib_sizes = nrow(moose.s), tp=-1))
}
(Eb1=cm.x[which.max(cm.x[,cri]),'E'])
(Eb2=cm.y[which.max(cm.y[,cri]),'E'])

# Select the optimal lag 
lagdata_x=make_lag(x=moose.s[,'Moose'],y=moose.s[,'Wolves'],lag.seq=lag.s,name=c('Moose','Wolves'))
lagdata_y=make_lag(x=moose.s[,'Wolves'],y=moose.s[,'Moose'],lag.seq=lag.s,name=c('Wolves','Moose'))

rho.terminal=NULL
for(i in 1:length(lag.s)){
   x.lag=lagdata_x[[i]]
   y.lag=lagdata_y[[i]]
   lsize.x=nrow(x.lag);lsize.y=nrow(y.lag)
   x_xmap_y <- ccm(x.lag, E = Eb1, 
                   lib_column = "Moose", target_column = "Wolves", 
                   lib_sizes = lsize.x, num_samples = 100,RNGseed = 17568)
   y_xmap_x <- ccm(y.lag, E = Eb2, 
                   lib_column = "Wolves", target_column = "Moose", 
                   lib_sizes = lsize.y, num_samples = 100,RNGseed = 17568)
   
   xxmapy <- ccm_means(x_xmap_y)
   yxmapx <- ccm_means(y_xmap_x)
   rho.terminal=rbind(rho.terminal,
                      c(lag=lag.s[i],
                        rho_x_term=xxmapy[nrow(xxmapy),cri],
                        rho_y_term=yxmapx[nrow(yxmapx),cri]
                      ))
}

(lag_xxmapy=rho.terminal[which.max(rho.terminal[,'rho_x_term']),'lag'])
(lag_yxmapx=rho.terminal[which.max(rho.terminal[,'rho_y_term']),'lag'])

# time series with optimal lags
nolagda=lagdata_x[[which(lag.s==0)]]
lagda_x_xmap_y=lagdata_x[[which(lag.s==lag_xxmapy)]]
lagda_y_xmap_x=lagdata_y[[which(lag.s==lag_yxmapx)]]




win.graph(600,200);par(mfcol=c(1,3))
# P1 No lag
lsize=unique(c(seq(5, nrow(nolagda), by = 2),nrow(nolagda)))
x_xmap_y <- ccm(nolagda, E = Eb1, 
                lib_column = "Moose", target_column = "Wolves", 
                lib_sizes = lsize, num_samples = 100,RNGseed = 17568)
y_xmap_x <- ccm(nolagda, E = Eb2, 
                lib_column = "Wolves", target_column = "Moose", 
                lib_sizes = lsize, num_samples = 100,RNGseed = 17568)

xxmapy <- ccm_means(x_xmap_y)
yxmapx <- ccm_means(y_xmap_x)

plot(rho~lib_size,xxmapy,
     xlab='Library size (L)',ylab='Correlation coefficient (rho)',
     main='Wolf-moose system (no lag)',col=4,
     type='l',lwd=2,ylim=c(min(c(xxmapy[,cri],yxmapx[,cri])),max(rho.terminal[,-1]))
)
lines(rho~lib_size,yxmapx,col=2,lwd=2)
legend('bottomright',
       c(paste('moose_xmap_wolf (top-down;lag=',0,')',sep=''),
         paste('wolf_xmap_moose (bottom-up;lag=',0,')',sep='')),
       col=c(4,2),lty=c(1,1),bty='n')

# P2 Lags vs skill
plot(rho_x_term~lag,rho.terminal,type='l',ylim=range(rho.terminal[,-1]),lwd=2,ylab='Correlation coefficient (rho)',col=4)
lines(rho_y_term~lag,rho.terminal,type='l',col=2,lwd=2)
abline(v=c(lag_xxmapy),lty=2,col=4);abline(v=c(lag_yxmapx),lty=2,col=2)
legend('bottomright',c('moose_xmap_wolf','wolf_xmap_moose'),
       col=c(4,2),lty=c(1,1),bty='n')

# P3 Lags
lsize=unique(c(seq(5, nrow(lagda_x_xmap_y), by = 2),nrow(lagda_x_xmap_y)))
x_xmap_y <- ccm(lagda_x_xmap_y, E = Eb1, 
                lib_column = "Moose", target_column = "Wolves", 
                lib_sizes = lsize, num_samples = 100,RNGseed = 17568)
lsize=unique(c(seq(5, nrow(lagda_y_xmap_x), by = 2),nrow(lagda_y_xmap_x)))
y_xmap_x <- ccm(lagda_y_xmap_x, E = Eb2, 
                lib_column = "Wolves", target_column = "Moose", 
                lib_sizes = lsize, num_samples = 100,RNGseed = 17568)

xxmapy <- ccm_means(x_xmap_y)
yxmapx <- ccm_means(y_xmap_x)

plot(rho~lib_size,xxmapy,
     xlab='Library size (L)',ylab='Correlation coefficient (rho)',
     main='Wolf-moose system (lagged)',col=4,
     type='l',lwd=2,ylim=c(min(c(xxmapy[,cri],yxmapx[,cri])),max(rho.terminal[,-1]))
)
lines(rho~lib_size,yxmapx,col=2,lwd=2)
legend('bottomright',
       c(paste('moose_xmap_wolf (top-down;lag=',lag_xxmapy,')',sep=''),
         paste('wolf_xmap_moose (bottom-up;lag=',lag_yxmapx,')',sep='')),
       col=c(4,2),lty=c(1,1),bty='n')



#######################################################################3
#### Hares-Lynx interaction
Hares=read.csv('Hares_Lynx.csv',header=T)
Hares.s <- Hares;Hares.s[,-1]=apply(Hares.s[,-1],2,scale)

# Select embedding dimensions based on hindcast ability
cm.x=cm.y=NULL
for(Ed in 2:8){
   cm.x=rbind(cm.x,ccm(Hares.s, E = Ed, random_libs = F,
                       lib_column = "Hares", target_column = "Lynx", 
                       lib_sizes = nrow(Hares.s), tp=-1))
   cm.y=rbind(cm.y,ccm(Hares.s, E = Ed, random_libs = F,
                       lib_column = "Lynx", target_column = "Hares", 
                       lib_sizes = nrow(Hares.s), tp=-1))
}
(Eb1=cm.x[which.max(cm.x[,cri]),'E'])
(Eb2=cm.y[which.max(cm.y[,cri]),'E'])


# Select the optimal lag 
lagdata_x=make_lag(x=Hares.s[,'Hares'],y=Hares.s[,'Lynx'],lag.seq=lag.s,name=c('Hares','Lynx'))
lagdata_y=make_lag(x=Hares.s[,'Lynx'],y=Hares.s[,'Hares'],lag.seq=lag.s,name=c('Lynx','Hares'))

rho.terminal=NULL
for(i in 1:length(lag.s)){
   x.lag=lagdata_x[[i]]
   y.lag=lagdata_y[[i]]
   lsize.x=nrow(x.lag);lsize.y=nrow(y.lag)
   x_xmap_y <- ccm(x.lag, E = Eb1, 
                   lib_column = "Hares", target_column = "Lynx", 
                   lib_sizes = lsize.x, num_samples = 100,RNGseed = 17568)
   y_xmap_x <- ccm(y.lag, E = Eb2, 
                   lib_column = "Lynx", target_column = "Hares", 
                   lib_sizes = lsize.y, num_samples = 100,RNGseed = 17568)
   
   xxmapy <- ccm_means(x_xmap_y)
   yxmapx <- ccm_means(y_xmap_x)
   rho.terminal=rbind(rho.terminal,
                      c(lag=lag.s[i],
                        rho_x_term=xxmapy[nrow(xxmapy),cri],
                        rho_y_term=yxmapx[nrow(yxmapx),cri]
                      ))
}

(lag_xxmapy=rho.terminal[which.max(rho.terminal[,'rho_x_term']),'lag'])
(lag_yxmapx=rho.terminal[which.max(rho.terminal[,'rho_y_term']),'lag'])

# time series with optimal lags
nolagda=lagdata_x[[which(lag.s==0)]]
lagda_x_xmap_y=lagdata_x[[which(lag.s==lag_xxmapy)]]
lagda_y_xmap_x=lagdata_y[[which(lag.s==lag_yxmapx)]]


win.graph(600,200);par(mfcol=c(1,3))
# P1 No lag
lsize=unique(c(seq(5, nrow(nolagda), by = 2),nrow(nolagda)))
x_xmap_y <- ccm(nolagda, E = Eb1, 
                lib_column = "Hares", target_column = "Lynx", 
                lib_sizes = lsize, num_samples = 100,RNGseed = 17568)
y_xmap_x <- ccm(nolagda, E = Eb2, 
                lib_column = "Lynx", target_column = "Hares", 
                lib_sizes = lsize, num_samples = 100,RNGseed = 17568)

xxmapy <- ccm_means(x_xmap_y)
yxmapx <- ccm_means(y_xmap_x)

plot(rho~lib_size,xxmapy,
     xlab='Library size (L)',ylab='Correlation coefficient (rho)',
     main='Lynx-Hares system (no lag)',col=4,
     type='l',lwd=2,ylim=c(min(c(xxmapy[,cri],yxmapx[,cri])),max(rho.terminal[,-1]))
     )
lines(rho~lib_size,yxmapx,col=2,lwd=2)
legend('bottomright',
       c(paste('Hares_xmap_Lynx (top-down;lag=',0,')',sep=''),
         paste('Lynx_xmap_Hares (bottom-up;lag=',0,')',sep='')),
       col=c(4,2),lty=c(1,1),bty='n')

# P3 Lags vs skill
plot(rho_x_term~lag,rho.terminal,type='l',ylim=range(rho.terminal[,-1]),lwd=2,ylab='Correlation coefficient (rho)',col=4)
lines(rho_y_term~lag,rho.terminal,type='l',col=2,lwd=2)
abline(v=c(lag_xxmapy),lty=2,col=4);abline(v=c(lag_yxmapx),lty=2,col=2)
legend('bottomright',c('Hares_xmap_Lynx','Lynx_xmap_Hares'),
       col=c(4,2),lty=c(1,1),bty='n')


# P2 Lags
lsize=unique(c(seq(5, nrow(lagda_x_xmap_y), by = 2),nrow(lagda_x_xmap_y)))
x_xmap_y <- ccm(lagda_x_xmap_y, E = Eb1, 
                lib_column = "Hares", target_column = "Lynx", 
                lib_sizes = lsize, num_samples = 100,RNGseed = 17568)
lsize=unique(c(seq(5, nrow(lagda_y_xmap_x), by = 2),nrow(lagda_y_xmap_x)))
y_xmap_x <- ccm(lagda_y_xmap_x, E = Eb2, 
                lib_column = "Lynx", target_column = "Hares", 
                lib_sizes = lsize, num_samples = 100,RNGseed = 17568)

xxmapy <- ccm_means(x_xmap_y)
yxmapx <- ccm_means(y_xmap_x)

plot(rho~lib_size,xxmapy,
     xlab='Library size (L)',ylab='Correlation coefficient (rho)',
     main='Lynx-Hares system (lagged)',col=4,
     type='l',lwd=2,ylim=c(min(c(xxmapy[,cri],yxmapx[,cri])),max(rho.terminal[,-1]))
     )
lines(rho~lib_size,yxmapx,col=2,lwd=2)
legend('bottomright',
       c(paste('Hares_xmap_Lynx (top-down;lag=',lag_xxmapy,')',sep=''),
         paste('Lynx_xmap_Hares (bottom-up;lag=',lag_yxmapx,')',sep='')),
       col=c(4,2),lty=c(1,1),bty='n')




