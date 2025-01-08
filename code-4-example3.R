#for figures in Example 3 in Super-efficiency-paper, on 11/28/24 by Y. Zuo
# named as code-4-example3.R
library(mvtnorm)
library(robustbase)

n=100; p=2#
epsilon=0.3; n1=floor(epsilon*n)
m1=rmvnorm(n, mean=(rep(0, p)),sigma=matrix(c(1, 0.88, 0.88, 1), ncol=2, byrow=T))
m111=m1
m1=m111

par(mfrow=c(1,2))  
plot(m1[,1], m1[,2], xlim=c(-3,8), ylim=c(-2, 7), xlab="x-axis", ylab="y-axis", pch=20)
fit0<-ltsReg(m1[,1:(p-1)], m1[,p])
abline(fit0$coefficients[[1]], fit0$coefficients[[2]], col=1, lty=1, lwd=2)
#abline(AA1_main(m1,1, 0, 500, 0.001), col=2, lty=2, lwd=1)
fit2<-lm(m1[,2]~m1[,1]) #p is assume to be 2
abline(fit2, col=2, lty=2, lwd=2)
fit3=lmrob(m1[,2]~m1[,1])
abline(fit3, lty=3, col=3, lwd=2)
lst_beta=lstReg(m1, 1, 1, 20)
abline(lst_beta, lty=4, col=4, lwd=2)
legend(x = "topleft", legend=c("LTS", "LS", "MM","LST"),    
       col=1:4, lty=1:4, cex=1,
       title="Line types", text.font=2, bg='lightblue')

if(n1>=1)
{
  m2=rmvnorm(n1, mean=c(7,-2),sigma=diag(rep(0.1, p)))
  m1[sample(1:n,n1),]<-m2 
}

m222=m1
m1=m222

plot(m1[,1], m1[,2], xlim=c(-3,8), ylim=c(-2, 7), xlab="x-axis", ylab="y-axis", pch=20)
fit1<-ltsReg(m1[,1:(p-1)], m1[,p])
abline(fit1$coefficients[[1]], fit1$coefficients[[2]], col=1, lty=1, lwd=1)
#abline(AA1_main(m1,1, 0, 500, 0.001), col=2, lty=2, lwd=1)
fit2<-lm(m1[,2]~m1[,1]) #p is assume to be 2
abline(fit2, col=2, lty=2, lwd=1)
fit3=lmrob(m1[,2]~m1[,1])
abline(fit3, lty=3, col=3, lwd=1)
lst_beta=lstReg(m1, 1, 1, 20)
abline(lst_beta, lty=4, col=4, lwd=2)
#text(4., -1.6, expression ('LTS'))
#text(6, 5.7, expression('LST'))
#text(4, -1, expression('LS'))

legend(x = "topleft", legend=c("LTS", "LS", "MM","LST"),
       col=1:4, lty=1:4, cex=1,
       title="Line types", text.font=2, bg='lightblue')
