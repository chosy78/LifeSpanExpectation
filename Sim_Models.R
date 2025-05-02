setwd('')
library(mvtnorm)

n=100 ## Change to 10000 for real run


alpha=1

nsim=5

MSE=matrix(NA,nr=60*27*3,nc=12)
Bias=matrix(NA,nr=60*27*3,nc=12)
MSE_outlier=matrix(NA,nr=60*27*3,nc=12)
Bias_outlier=matrix(NA,nr=60*27*3,nc=12)

param=matrix(NA,nr=60*27*3,nc=7)
kk=1
sigma2=0.5
delta=0.3

set.seed(1)

for (sigma in c(1:5)*.4){
  for (rho in c(0,0.2,0.5,0.8)){
    for (p in c(10,20,50)){
      for (p1 in c(0,1,3)){
        for (p2 in c(0,1,3)){
          for (p3 in c(0,1,3)){
            for (sigma2 in c(0.5,1,2)){
             
        
      
      #sigma=1
      #rho=0.5
      
      
      Sigma=matrix(rho,nr=p,nc=p)
      diag(Sigma)=1
      beta=rep(alpha,p)
      
      res=matrix(NA,nr=nsim,nc=12)
      out21=rep(0,nsim)
      out22=rep(0,nsim)
      out23=rep(0,nsim)
      out31=rep(0,nsim)
      out32=rep(0,nsim)
      out33=rep(0,nsim)
      for (sim in 1:nsim){
        X=rmvnorm(n,sigma=Sigma)
        e=rnorm(n)*sigma
        Y=X%*%beta+e
        X0=rmvnorm(1,sigma=Sigma)
        e0=rnorm(1)*sigma
        Y0=X0%*%beta+e0
        a1=lm(Y~X)
        if (p1+p2>0){
          a3=lm(Y~X[,-(1:(p1+p2))])
        }else{
          a3=a1
        }
        
        X01=X0
        X02=X0
        X03=X0
        if (p2+p3>0){
          X02[(p1+1):(p1+p2+p3)]=X02[(p1+1):(p1+p2+p3)]+rnorm(p2+p3)*sigma2
          X03[(p1+1):(p1+p2+p3)]=X03[(p1+1):(p1+p2+p3)]+rnorm(p2+p3)*sigma2+delta
        } 
        
        Y11=X01%*%coef(a1)[-1]+coef(a1)[1]
        Y12=X02%*%coef(a1)[-1]+coef(a1)[1]
        Y13=X03%*%coef(a1)[-1]+coef(a1)[1]
        ## Method 1 X<0 vs X>0##
        ### Accurate X01
        if (p1+p2>0){
        ind=which(X01[1:(p1+p2)]<0)
        if (length(ind)==0){
          Y21=Y11
        }else{a2=lm(Y~X[,-ind])
        out21[sim]=1
        Y21=X01[-ind]%*%coef(a2)[-1]+coef(a2)[1]}
        ### Noise X02
        ind=which(X02[1:(p1+p2)]<0)
        if (length(ind)==0){
          Y22=Y12
        }else{a2=lm(Y~X[,-ind])
        out22[sim]=1
        Y22=X02[-ind]%*%coef(a2)[-1]+coef(a2)[1]}
        ### Noise X03 with Bias
        ind=which(X03[1:(p1+p2)]<0)
        if (length(ind)==0){
          Y23=Y13
        }else{a2=lm(Y~X[,-ind])
        out23[sim]=1
        Y23=X03[-ind]%*%coef(a2)[-1]+coef(a2)[1]}
        ## Method 2 |X|>1 vs |X|<1##
        ### Accurate X01
        ind=which(abs(X01[1:(p1+p2)])>1)
        if (length(ind)==0){
          Y31=Y11
        }else{a2=lm(Y~X[,-ind])
        out31[sim]=1
        Y31=X01[-ind]%*%coef(a2)[-1]+coef(a2)[1]}
        ### Noise X02
        ind=which(abs(X02[1:(p1+p2)])>1)
        if (length(ind)==0){
          Y32=Y12
        }else{a2=lm(Y~X[,-ind])
        out32[sim]=1
        Y32=X02[-ind]%*%coef(a2)[-1]+coef(a2)[1]}
        ### Noise X03 with Bias
        ind=which(abs(X03[1:(p1+p2)])>1)
        if (length(ind)==0){
          Y33=Y13
        }else{a2=lm(Y~X[,-ind])
        out33[sim]=1
        Y33=X03[-ind]%*%coef(a2)[-1]+coef(a2)[1]}
        ## Method 3 X##
        #Y11=X01%*%coef(a1)[-1]+coef(a1)[1]
        #Y12=X02%*%coef(a1)[-1]+coef(a1)[1]
        #Y13=X03%*%coef(a1)[-1]+coef(a1)[1]
        ## Method 4 1##
        Y41=X01[,-(1:(p1+p2))]%*%coef(a3)[-1]+coef(a3)[1]
        Y42=X02[,-(1:(p1+p2))]%*%coef(a3)[-1]+coef(a3)[1]
        Y43=X03[,-(1:(p1+p2))]%*%coef(a3)[-1]+coef(a3)[1]
        }else{
          Y21=Y11
          Y31=Y11
          Y41=Y11
          Y22=Y12
          Y32=Y12
          Y42=Y12
          Y23=Y13
          Y33=Y13
          Y43=Y13
        }
        res[sim,]=c(Y21,Y22,Y23,Y31,Y32,Y33,Y11,Y12,Y13,Y41,Y42,Y43)-Y0
      }
      
      for (j in 1:12){
        MSE[kk,j]=mean(res[,j]^2)
        Bias[kk,j]=mean(res[,j])
        MSE_outlier[kk,j]=mean(res[,j]^2)
        Bias_outlier[kk,j]=mean(res[,j])
      }
      MSE_outlier[kk,1]=mean(res[which(out21==1),1]^2)
      Bias_outlier[kk,1]=mean(res[which(out21==1),1])
      MSE_outlier[kk,2]=mean(res[which(out22==1),2]^2)
      Bias_outlier[kk,2]=mean(res[which(out22==1),2])
      MSE_outlier[kk,3]=mean(res[which(out23==1),3]^2)
      Bias_outlier[kk,3]=mean(res[which(out23==1),3])
      MSE_outlier[kk,4]=mean(res[which(out31==1),4]^2)
      Bias_outlier[kk,4]=mean(res[which(out31==1),4])
      MSE_outlier[kk,5]=mean(res[which(out32==1),5]^2)
      Bias_outlier[kk,5]=mean(res[which(out32==1),5])
      MSE_outlier[kk,6]=mean(res[which(out33==1),6]^2)
      Bias_outlier[kk,6]=mean(res[which(out33==1),6])
      
      param[kk,]=c(sigma,rho,p,p1,p2,p3,sigma2)
      kk=kk+1
      print(kk)
          }}}}
    }
  }
}



write.csv(MSE,'MSE_linear.csv',row.names=F,quote=F)
write.csv(Bias,'Bias_linear.csv',row.names=F,quote=F)
write.csv(MSE_outlier,'MSE_o_linear.csv',row.names=F,quote=F)
write.csv(Bias_outlier,'Bias_o_linear.csv',row.names=F,quote=F)

MSE1=MSE_outlier
Bias1=Bias_outlier






library(mvtnorm)

n=1000


alpha=1

nsim=1000

MSE=matrix(NA,nr=160,nc=4)
Bias=matrix(NA,nr=160,nc=4)
kk=1

for (sigma in c(1:10)*.2){
  for (rho in c(0,0.2,0.5,0.8)){
    for (p in c(1,2,5,10)){
      
      #sigma=1
      #rho=0.5
      
      
      Sigma=matrix(rho,nr=p,nc=p)
      diag(Sigma)=1
      beta=rep(alpha,p)
      
      res=matrix(NA,nr=nsim,nc=4)
      for (sim in 1:nsim){
        X=rmvnorm(n,sigma=Sigma)
        X1=X
        e=rnorm(n)*sigma
        Y=X%*%beta+(X[,1]>1)*(X[,1]-1)*alpha+(X[,1]<=-1)*(X[,1]+1)*alpha+e
        X0=rmvnorm(1,sigma=Sigma)
        e0=rnorm(1)*sigma
        Y0=X0%*%beta+e0
        a1=lm(Y~X1)
        if (p==1){
          a2=lm(Y~1)
        }else{
          a2=lm(Y~X1[,-1])
        }
        
        Y1=X0%*%coef(a1)[-1]+coef(a1)[1]
        Y2=X0[-1]%*%coef(a2)[-1]+coef(a2)[1]
        ## Method 1 X<0 vs X>0##
        if (X0[1]<0){
          res[sim,1]=Y2-Y0
        } else {
          res[sim,1]=Y1-Y0
        }
        ## Method 2 |X|>1 vs |X|<1##
        if (abs(X0[1])<1){
          res[sim,2]=Y2-Y0
        } else {
          res[sim,2]=Y1-Y0
        }
        ## Method 3 X##
        res[sim,3]=Y1-Y0
        ## Method 4 1##
        res[sim,4]=Y2-Y0
      }
      
      MSE[kk,1]=mean(res[,1]^2)
      MSE[kk,2]=mean(res[,2]^2)
      MSE[kk,3]=mean(res[,3]^2)
      MSE[kk,4]=mean(res[,4]^2)
      Bias[kk,1]=mean(res[,1])
      Bias[kk,2]=mean(res[,2])
      Bias[kk,3]=mean(res[,3])
      Bias[kk,4]=mean(res[,4])
      kk=kk+1
      print(kk)
    }
  }
}

write.csv(MSE,'MSE_nonlinear.csv',row.names=F,quote=F)
write.csv(Bias,'Bias_nonlinear.csv',row.names=F,quote=F)

MSE2=MSE
Bias2=Bias

SSS=matrix(NA,nr=160,nc=3)
kk=1

for (sigma in c(1:10)*.2){
  for (rho in c(0,0.2,0.5,0.8)){
    for (p in c(1,2,5,10)){
      SSS[kk,1]=sigma
      SSS[kk,2]=rho
      SSS[kk,3]=p
      kk=kk+1
    }}}

sig=c(1:10)*.2
pdf('MSE_linear.pdf',width=20,height=20)
par(mfrow=c(4,4),mar=c(6,6,6,3))
for (rho in c(0,0.2,0.5,0.8)){
  for (p in c(1,2,5,10)){
    if (p==1 & rho!=0){
      plot(1, type = "n", xlab = "",
           ylab = "", xlim = c(0, 2),
           ylim = c(0, 3))
    } else{
      ind=which(SSS[,2]==rho & SSS[,3]==p)
      dat=MSE1[ind,]
      sig=SSS[ind,1]
      plot(dat[,1]~sig,type='l',col='red',lwd=2,ylim=c(0,5),xlab='Noise',ylab='MSE',main=paste0('rho = ',rho,', p = ',p,sep=''),cex.lab=2,cex.main=3)
      par(new=T)
      plot(dat[,2]~sig,type='l',col='blue',ylim=c(0,5),lwd=2,xlab='',ylab='',main='',cex.lab=2,cex.main=3)
      par(new=T)
      plot(dat[,3]~sig,type='l',col='green',ylim=c(0,5),lwd=2,xlab='',ylab='',main='',cex.lab=2,cex.main=3)
      par(new=T)
      plot(dat[,4]~sig,type='l',col='black',ylim=c(0,5),lwd=2,xlab='',ylab='',main='',cex.lab=2,cex.main=3)
      legend('topleft',legend = c('Method 1','Method 2','Method 3','Method 4'),lty=1,lwd=2,col=c('red','blue','green','black'),cex=1.5)
    }
  }
}
dev.off()



pdf('MSE_nonlinear.pdf',width=20,height=20)
par(mfrow=c(4,4),mar=c(6,6,6,3))
for (rho in c(0,0.2,0.5,0.8)){
  for (p in c(1,2,5,10)){
    if (p==1 & rho!=0){
      plot(1, type = "n", xlab = "",
           ylab = "", xlim = c(0, 2),
           ylim = c(0, 3))
    } else{
      ind=which(SSS[,2]==rho & SSS[,3]==p)
      dat=MSE2[ind,]
      sig=SSS[ind,1]
      plot(dat[,1]~sig,type='l',col='red',lwd=2,ylim=c(0,5),xlab='Noise',ylab='MSE',main=paste0('rho = ',rho,', p = ',p,sep=''),cex.lab=2,cex.main=3)
      par(new=T)
      plot(dat[,2]~sig,type='l',col='blue',ylim=c(0,5),lwd=2,xlab='',ylab='',main='',cex.lab=2,cex.main=3)
      par(new=T)
      plot(dat[,3]~sig,type='l',col='green',ylim=c(0,5),lwd=2,xlab='',ylab='',main='',cex.lab=2,cex.main=3)
      par(new=T)
      plot(dat[,4]~sig,type='l',col='black',ylim=c(0,5),lwd=2,xlab='',ylab='',main='',cex.lab=2,cex.main=3)
      legend('topleft',legend = c('Method 1','Method 2','Method 3','Method 4'),lty=1,lwd=2,col=c('red','blue','green','black'),cex=1.5)
    }
  }
}
dev.off()



sig=c(1:10)*.2
pdf('Bias_linear.pdf',width=20,height=20)
par(mfrow=c(4,4),mar=c(6,6,6,3))
for (rho in c(0,0.2,0.5,0.8)){
  for (p in c(1,2,5,10)){
    if (p==1 & rho!=0){
      plot(1, type = "n", xlab = "",
           ylab = "", xlim = c(0, 2),
           ylim = c(0, 3))
    } else{
      ind=which(SSS[,2]==rho & SSS[,3]==p)
      dat=Bias1[ind,]
      sig=SSS[ind,1]
      plot(dat[,1]~sig,type='l',col='red',lwd=2,ylim=c(-1,1),xlab='Noise',ylab='Bias',main=paste0('rho = ',rho,', p = ',p,sep=''),cex.lab=2,cex.main=3)
      par(new=T)
      plot(dat[,2]~sig,type='l',col='blue',ylim=c(-1,1),lwd=2,xlab='',ylab='',main='',cex.lab=2,cex.main=3)
      par(new=T)
      plot(dat[,3]~sig,type='l',col='green',ylim=c(-1,1),lwd=2,xlab='',ylab='',main='',cex.lab=2,cex.main=3)
      par(new=T)
      plot(dat[,4]~sig,type='l',col='black',ylim=c(-1,1),lwd=2,xlab='',ylab='',main='',cex.lab=2,cex.main=3)
      abline(h=0,lty=2)
      legend('bottomleft',legend = c('Method 1','Method 2','Method 3','Method 4'),lty=1,lwd=2,col=c('red','blue','green','black'),cex=1.5)
    }
  }
}
dev.off()



sig=c(1:10)*.2
pdf('Bias_nonlinear.pdf',width=20,height=20)
par(mfrow=c(4,4),mar=c(6,6,6,3))
for (rho in c(0,0.2,0.5,0.8)){
  for (p in c(1,2,5,10)){
    if (p==1 & rho!=0){
      plot(1, type = "n", xlab = "",
           ylab = "", xlim = c(0, 2),
           ylim = c(0, 3))
    } else{
      ind=which(SSS[,2]==rho & SSS[,3]==p)
      dat=Bias2[ind,]
      sig=SSS[ind,1]
      plot(dat[,1]~sig,type='l',col='red',lwd=2,ylim=c(-1,1),xlab='Noise',ylab='Bias',main=paste0('rho = ',rho,', p = ',p,sep=''),cex.lab=2,cex.main=3)
      par(new=T)
      plot(dat[,2]~sig,type='l',col='blue',ylim=c(-1,1),lwd=2,xlab='',ylab='',main='',cex.lab=2,cex.main=3)
      par(new=T)
      plot(dat[,3]~sig,type='l',col='green',ylim=c(-1,1),lwd=2,xlab='',ylab='',main='',cex.lab=2,cex.main=3)
      par(new=T)
      plot(dat[,4]~sig,type='l',col='black',ylim=c(-1,1),lwd=2,xlab='',ylab='',main='',cex.lab=2,cex.main=3)
      abline(h=0,lty=2)
      legend('bottomleft',legend = c('Method 1','Method 2','Method 3','Method 4'),lty=1,lwd=2,col=c('red','blue','green','black'),cex=1.5)
    }
  }
}
dev.off()