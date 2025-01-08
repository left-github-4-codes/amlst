# code for example 2 amstat by Y. Zuo on 11/23/2023
# changed from code-4-example2&3-v1.R to code-4-table1-v1.R for 
# robust penalized lST for high-dimensional data paper # on 12/18/2023
# now this been used for Super-efficiency article and named as code-4-table-1.R on 11/21/24 by Y.zuo

library(mvtnorm)
library(robustbase)

R=RepN=1000
n=200 #tuning this n for your need
p=30 #tuning
alpha=3
N_ls=1
delta=1
# Cc=6
# k=6
# cutoff=0.0001
#beta0=rep(1, p)
beta0=c(rep(1, 2), rep(0, p-2)) #other choices (1, rep(0, p-1)), (1,1, rep(0, p-2))
#beta0=c(rep(0,5), rep(0,p-5))

beta_aa1=beta_aa2=beta_aa3=beta_aa4=matrix(0, nrow=R, ncol=p)  #aa1-mm, aa2-lst, aa3-LTS, aa4-ls
t_aa1=t_aa2=t_aa3=t_aa4=rep(0, R)

epsilon=0.02  #tuning
m=n1=floor(epsilon*n)

for (i in 1:R)
{ ss=matrix(0.0, nrow=p, ncol=p); diag(ss)<-1
m1=rmvnorm(n, mean=(rep(0, p)),ss)
#m1=rmvnorm(n, mean=(rep(0, p)),sigma=matrix(c(1, 0.9, 0.9, 1), ncol=2, byrow=T))
#m1=rmvnorm(n, mean=(rep(0, p)), sigma=diag(rep(1,p)))

# #comtamination scheme I  
# if(n1>=1)
# {
#   mu=rep(7,p); mu[p]=-2
#   m2=rmvnorm(n1, mu,sigma=diag(rep(0.1, p)))
#   m1[sample(1:n,n1),]<-m2 
# }
z=m1
# 

#contamination scheme II 
#  x <- matrix(rnorm(n*(p-1), sigma),nrow=n, ncol=p-1)
#  x<-cbind(rep(1,n), x) #n by p matrix, this should be the design of the paper
# e <- rnorm(n,0,1)  # error terms
#eout[1:m] <- eout[1:m] + 10 # vertical outliers

#Contamination scheme one:
# eout <- e;  #eout[1:m] <- eout[1:m] + 20 # vertical outliers
# yout <- c(x %*% matrix(beta0) + sigma * eout) # response
# y <- yout; #xout[1:m,] <- xout[1:m,] + 20 # bad leverage points

#Contamination scheme two:
# eout <- e; eout[1:m] <- eout[1:m] + 20 # vertical outliers
# ind<-sample(n, m)
# yout <- c(x %*% beta0 + sigma * eout) # response
# xout <- x; xout[ind,] <- xout[ind,] + 20 # bad leverage points

#contamination scheme three (based on the proof of Theorem 4.2)
# eout <- e; eout[1:m] <- eout[1:m] + 20 # vertical outliers
# e <- rnorm(n,0,1) 
# ind<-sample(n, m); Del=3.5 # 10^1;
# Kap= 3 #10^1 # yout <- c(x %*% beta0 + sigma * eout) # response
# x<- m1[,1:(p-1)] #;  y<-m1[,p] 
# xx<-cbind(rep(1,n), x)
# y<-xx%*%matrix(bet0)+matrix(e)
# x[ind,]<- c(rep(Del, p-1)) #c(Del,rep(0,(p-2)));
# y[ind]<- 1*Del
#z=cbind(x,y)  
#m1=z

# ind<-sample(n, m); 
#  Del=7 # tuning from 3, 4, 4.5, to 7
# # # yout <- c(x %*% beta0 + sigma * eout) # response
#  x<- m1[,1:(p-1)];  y<-m1[,p]
#  x[ind,]<-c(Del, rep(Del,(p-2))); y[ind]<- -1*Del^1
#  z=cbind(x,y)  

e <- rnorm(n,0,1) 
ind<-sample(n, m); Del=4.5;
# # #Kap= 3 #10^1 # yout <- c(x %*% beta0 + sigma * eout) # response
x<- m1[,1:(p-1)] #;  y<-m1[,p] 
 xx<-cbind(rep(1,n), x)
y<-xx%*%matrix(beta0)+matrix(e)
x[ind,]<- c(rep(Del, p-1)) #c(Del,rep(0,(p-2)));
y[ind]<- -1*Del
z=cbind(x,y)  


t1=Sys.time()  
# beta_aa1[i,]= AA1_main_new_lst_v2(z,alpha, delta, N_ls)
fit0=lmrob(z[,p]~z[,1:(p-1)])        #MM
beta_aa1[i,]=as.numeric(fit0$coefficients)
t2=Sys.time()-t1
t_aa1[i]=t2
# 
t1=Sys.time()  
beta_aa2[i,]=lstReg(z, alpha, delta, N_ls) #LST
t2=Sys.time()-t1
t_aa2[i]=t2

t1=Sys.time()  
fit1<-ltsReg(z[,1:(p-1)], z[,p])
beta_aa3[i,]=as.numeric(fit1$coefficients) #LTS
t2=Sys.time()-t1
t_aa3[i]=t2

t1=Sys.time()  
fit2<-lm(z[,p]~z[,1:(p-1)])
beta_aa4[i,]=as.numeric(fit2$coefficients) #LS
t2=Sys.time()-t1
t_aa4[i]=t2

print(i)
}

beta_aa1_mean=colMeans(beta_aa1)
deviat_beta_aa1=beta_aa1-matrix(beta_aa1_mean, byrow=T, nrow=RepN, ncol=p)
beta_aa2_mean=colMeans(beta_aa2)
deviat_beta_aa2=beta_aa2-matrix(beta_aa2_mean, byrow=T,nrow=RepN, ncol=p)
beta_aa3_mean=colMeans(beta_aa3)
deviat_beta_aa3=beta_aa3-matrix(beta_aa3_mean, byrow=T,nrow=RepN, ncol=p)
beta_aa4_mean=colMeans(beta_aa4)
deviat_beta_aa4=beta_aa4-matrix(beta_aa4_mean, byrow=T,nrow=RepN, ncol=p)


SVAR_aa1=(deviat_beta_aa1*deviat_beta_aa1) #/(RepN-1)
SVAR_aa2=(deviat_beta_aa2*deviat_beta_aa2) #/(RepN-1)
SVAR_aa3=(deviat_beta_aa3*deviat_beta_aa3) #/(RepN-1)
SVAR_aa4=(deviat_beta_aa4*deviat_beta_aa4) #/(RepN-1)

# SVAR_aa1=sum(deviat_beta_aa1*deviat_beta_aa1)/(RepN-1)
# SVAR_aa2=sum(deviat_beta_aa2*deviat_beta_aa2)/(RepN-1)
# SVAR_aa3=sum(deviat_beta_aa3*deviat_beta_aa3)/(RepN-1)
# SVAR_aa4=sum(deviat_beta_aa4*deviat_beta_aa4)/(RepN-1)
SVAR_aa11=sum(SVAR_aa1)/(R-1)
SVAR_aa21=sum(SVAR_aa2)/(R-1)
SVAR_aa31=sum(SVAR_aa3)/(R-1)
SVAR_aa41=sum(SVAR_aa4)/(R-1)

if (sum(beta0)==0) 
{EMSE_aa1=(beta_aa1*beta_aa1) #/R  #treat beta0 as 0 vector
EMSE_aa2=(beta_aa2*beta_aa2) #/R
EMSE_aa3=(beta_aa3*beta_aa3) #/R
EMSE_aa4=(beta_aa4*beta_aa4) #/R

# EMSE_aa1=sum(beta_aa1*beta_aa1)/R  #treat beta0 as 0 vector
# EMSE_aa2=sum(beta_aa2*beta_aa2)/R
# EMSE_aa3=sum(beta_aa3*beta_aa3)/R
# EMSE_aa4=sum(beta_aa4*beta_aa4)/R

EMSE_aa11=sum(EMSE_aa1)/R
EMSE_aa21=sum(EMSE_aa2)/R
EMSE_aa31=sum(EMSE_aa3)/R
EMSE_aa41=sum(EMSE_aa4)/R
} 
if (sum(beta0)!=0)
{ print("beta0")
  beta_matrix=matrix(beta0, byrow=T, nrow=R, ncol=p) #this is situation when beta0 is given and non-zero
  beta=beta_matrix
  beta_aa1=beta_aa1-beta
  beta_aa2=beta_aa2-beta
  beta_aa3=beta_aa3-beta
  beta_aa4=beta_aa4-beta
  # 
  EMSE_aa1=(beta_aa1*beta_aa1) #/R  #treat beta0 as 0 vector
  EMSE_aa2=(beta_aa2*beta_aa2) #/R
  EMSE_aa3=(beta_aa3*beta_aa3) #/R
  EMSE_aa4=(beta_aa4*beta_aa4) #/R
  
  EMSE_aa11=sum(EMSE_aa1)/R
  EMSE_aa21=sum(EMSE_aa2)/R
  EMSE_aa31=sum(EMSE_aa3)/R
  EMSE_aa41=sum(EMSE_aa4)/R
}  

t_aa11=sum(t_aa1)
t_aa21=sum(t_aa2)
t_aa31=sum(t_aa3)
t_aa41=sum(t_aa4)

print(c(R,n,p,N_ls,alpha,epsilon, Del) )
print( c(EMSE_aa11, EMSE_aa21, EMSE_aa31, EMSE_aa41))
print( c(SVAR_aa11, SVAR_aa21, SVAR_aa31, SVAR_aa41))
print(c(t_aa11,t_aa21,t_aa31,t_aa41))
print(EMSE_aa41/c(EMSE_aa11, EMSE_aa21, EMSE_aa31,  EMSE_aa41))
print(SVAR_aa41/c(SVAR_aa11, SVAR_aa21, SVAR_aa31,  SVAR_aa41))
 
# EMSE_aa1=sum(beta_aa1*beta_aa1)/R
# EMSE_aa2=sum(beta_aa2*beta_aa2)/R
# EMSE_aa3=sum(beta_aa3*beta_aa3)/R
# #EMSE_aa4=sum(beta_aa3*beta_aa3)/R
# EMSE_aa4=sum(beta_aa4*beta_aa4)/R #modified on 04/20/24 by Y.Zuo during revision of TAS
# 
# print(c(R,n,p,N_ls,alpha,epsilon, cutoff, Del, Kap) )
# print( c(#EMSE_aa1, 
#   EMSE_aa2, EMSE_aa3, EMSE_aa4))
# print(c(#t_aa1,
#   t_aa2,t_aa2,t_aa4))
# print(EMSE_aa4/c(#EMSE_aa1, 
#   EMSE_aa2, EMSE_aa3,  EMSE_aa4))


#the following code is for the boxplots

# #SS_aa1=beta_aa1*beta_aa1
# SS_aa2=beta_aa2*beta_aa2
# SS_aa3=beta_aa3*beta_aa3
# #SS_aa4=beta_aa3*beta_aa3
# SS_aa4=beta_aa4*beta_aa4  #modified on 04/20/24 by Y.Zuo during revision of TAS

TT_aa1=t_aa1
TT_aa2=t_aa2
TT_aa3=t_aa3
TT_aa4=t_aa4

RE_aa1=EMSE_aa41/EMSE_aa11
RE_aa2=EMSE_aa41/EMSE_aa21
RE_aa3=EMSE_aa41/EMSE_aa31
RE_aa4=EMSE_aa41/EMSE_aa41

par(mfrow = c(2,2))
boxplot(cbind(EMSE_aa1[1:1000], 
  EMSE_aa2[1:1000], EMSE_aa3[1:1000], EMSE_aa4[1:1000]), 
  names=c("MM","LST","LTS","LS"),
  xlab="procedure",
  ylab="EMSE",
  color=c("magenta","blue","green", "purple"),
  border="red"
)

boxplot(cbind(SVAR_aa1[1:1000], 
  SVAR_aa2[1:1000], SVAR_aa3[1:1000], SVAR_aa4[1:1000]), 
  names=c("MM", 
    "LST","LTS", "LS"),
  xlab="procedure",
  ylab="SVAR",
  color=c("magenta","blue", "green", "purple"),
  border="blue"
)

boxplot(cbind(TT_aa1, TT_aa2, TT_aa3, TT_aa4), 
        names=c("MM","LST","LTS", "LS"),
        xlab="procedure",
        ylab="TT(seconds)",
        color=c("magenta","blue", "green", "purple"),
        border="green"
)
boxplot(cbind(RE_aa1, 
  RE_aa2, RE_aa3, RE_aa4), 
  names=c("MM", 
    "LST","LTS", "LS"),
  xlab="procedure",
  ylab="RE",
  color=c("magenta","blue", "green", "purple"),
  border="magenta"
)