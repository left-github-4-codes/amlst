# revision for TAS by Y. Zuo on 04/18/24
x=c(5, 5.5, 4, 3.5, 3, 2.5, -2)
y=c(-.5, -.5, 6, 4, 2.4, 2, .5)
m1=Zn=cbind(x,y)
x=m1[,1]; y=m1[,2]
par(mfrow=c(1,2))
plot(x, y, pch=16, xlim=c(-3,8), ylim=c(-2, 7))
abline(0,0, lty=1, col=1, lwd=1)
abline(0, 1, lty=2, col=2, lwd=1)
text(x, y-0.3, as.character(seq(1:7)))
legend(x = "topleft", legend=c("y=0", "y=x"), col=1:2, lty=1:2, cex=1, text.font=3, bg='lightblue',title="reference lines")

cs=c(1,2,7); m1=Zn
m1<-m1[-cs,]

alpha=3
#lst_beta= lstRegv1(m1,3,1,20)
lst_beta=AA1_main_new_lst_v2(m1,3,1,20)
#fit3<-lm(m1[,2]~m1[,1]) #p is assume to be 2
abline(lst_beta,lty=2, col=2, lwd=1) #abline(fit3, col=3, lty=3, lwd=1)

#for RSS of LTS based on ls calculation
fit3<-lm(m1[,2]~m1[,1]) #p is assume to be 2
abline(fit3, col=3, lty=3, lwd=1)
ls_bet=fit3$coefficients
RS=resid(fit3) ^ 2
index=sort(RS)[1:5]
RSS=sum(index)
print(RSS)

##################### for RSS of LTS based on the LST line 
lst_beta= lstRegv1(m1,3,1,20)
bt = lst_beta
#bt=ls_bet
z=m1;p=2
rn = z[, p] -data.matrix(cbind(rep(1, n), z[, 1:(p - 1)])) %*% as.matrix(bt)
RS=rn^2
index=sort(RS)[1:4]
RSS=sum(index)
print(RSS)

RSS=deviance(fit3)
print(c(RSS, cs))

#for RSS of LST
#alpha=3
z=Zn; n=7; p=2
rn = rep(0, n)
I_bet = rep(0, n)
sd = rep(1, n)
#bt=c(0,0)
bt = lst_beta
#bt=c(1,1)

rn = z[, p] -data.matrix(cbind(rep(1, n), z[, 1:(p - 1)])) %*% as.matrix(bt)
med = median(rn)
madd = mad(rn)
abdv = abs(rn - med)
sd = abdv / madd
indx_bet = which(sd <= alpha)
K = length(indx_bet)
print("K"); print(K)  #total number of squared residuals employed
RSS=sum(sd^2) #RSS=deviance(fit3)
print(c(RSS, cs))

######################### 04/16/24 calculate the SSR for y=x or y=0 foe the 7 points data set in Figure 1 of ZZ23
x=c(5, 5.5, 4, 3.5, 3, 2.5, -2)
y=c(-.5, -.5, 6, 4, 2.4, 2, .5)
m1=Zn=cbind(x,y)

LTS
Q1=sum(z[c(1,2,6,7),][,2]^2)
print(Q1)

vv=z[c(3,4,5,6),]
dxy=vv[,1]-vv[,2]
Q2=sum(dxy^2)
print(Q2)

LST
Q1=sum(z[c(4,5,6,7),][,2]^2)
print(Q1)

vv=z[c(3,2,5,4),]
dxy=vv[,1]-vv[,2]
Q2=sum(dxy^2)
print(Q2)