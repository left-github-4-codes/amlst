library(mvtnorm)
library(robustbase)
library(compiler)
enableJIT(1)

R=RepN=1000
n=300 #tuning this n for your need
p=30
alpha=3
N_ls=1
delta=1
Cc=6
k=6
cutoff=0.0001
#bet0=rep(0, p)
bet0=c(rep(1,15), rep(-1,p-15))

beta_aa1=beta_aa2=beta_aa3=beta_aa4=matrix(0, nrow=R, ncol=p)  #aa1-lst, aa2-wls, aa3-LTS, aa3-ls
t_aa1=t_aa2=t_aa3=t_aa4=0

epsilon=0.00
m=n1=floor(epsilon*n)

for (i in 1:R)
{ ss=matrix(0.9, nrow=p, ncol=p); diag(ss)<-1
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
# z=m1
# 

#contamination scheme II 
#  x <- matrix(rnorm(n*(p-1), sigma),nrow=n, ncol=p-1)
#  x<-cbind(rep(1,n), x) #n by p matrix, this should be the design of the paper
#e <- rnorm(n,0,1)  # error terms
#eout[1:m] <- eout[1:m] + 10 # vertical outliers

#Contamination scheme one:
# eout <- e; eout[1:m] <- eout[1:m] + 20 # vertical outliers
# yout <- c(x %*% beta0 + sigma * eout) # response
# xout <- x; xout[1:m,] <- xout[1:m,] + 20 # bad leverage points

#Contamination scheme two:
# eout <- e; eout[1:m] <- eout[1:m] + 20 # vertical outliers
# ind<-sample(n, m)
# yout <- c(x %*% beta0 + sigma * eout) # response
# xout <- x; xout[ind,] <- xout[ind,] + 20 # bad leverage points

#contamination scheme three (based on the proof of Theorem 4.2)
# eout <- e; eout[1:m] <- eout[1:m] + 20 # vertical outliers
e <- rnorm(n,0,1) 
ind<-sample(n, m); Del=4.5 # 10^1;
Kap= 4 #10^1 # yout <- c(x %*% beta0 + sigma * eout) # response
x<- m1[,1:(p-1)] #;  y<-m1[,p] 
xx<-cbind(rep(1,n), x)
y<-xx%*%matrix(bet0)+matrix(e)
x[ind,]<- c(rep(Del, p-1)) #c(Del,rep(0,(p-2)));
y[ind]<- -1*Del
z=cbind(x,y)  


# t1=Sys.time()  
# beta_aa1[i,]= AA1_main_new_lst_v2(z,alpha, delta, N_ls)
# t2=Sys.time()-t1
# t_aa1=t_aa1+t2
# 
# t1=Sys.time()  
# beta_aa2[i,]=wls(z, Cc, k, cutoff)
# t2=Sys.time()-t1
# t_aa2=t_aa2+t2

t1=Sys.time()  
beta_aa2[i,]=lstReg(z, alpha, delta, N_ls)
t2=Sys.time()-t1
t_aa2=t_aa2+t2


t1=Sys.time()  
fit1<-ltsReg(z[,1:(p-1)], z[,p])
beta_aa3[i,]=as.numeric(fit1$coefficients) 
t2=Sys.time()-t1
t_aa3=t_aa3+t2

t1=Sys.time()  
fit2<-lm(z[,p]~z[,1:(p-1)])
beta_aa4[i,]=as.numeric(fit2$coefficients) 
t2=Sys.time()-t1
t_aa4=t_aa4+t2

print(i)
}

# EMSE_aa1=sum(beta_aa1*beta_aa1)/R
# EMSE_aa2=sum(beta_aa2*beta_aa2)/R
# EMSE_aa3=sum(beta_aa3*beta_aa3)/R
# EMSE_aa4=sum(beta_aa4*beta_aa4)/R

#beta_aa1_mean=colMeans(beta_aa1)
#deviat_beta_aa1=beta_aa1-matrix(beta_aa1_mean, byrow=T, nrow=RepN, ncol=p)
beta_aa2_mean=colMeans(beta_aa2)
deviat_beta_aa2=beta_aa2-matrix(beta_aa2_mean, byrow=T,nrow=RepN, ncol=p)
beta_aa3_mean=colMeans(beta_aa3)
deviat_beta_aa3=beta_aa3-matrix(beta_aa3_mean, byrow=T,nrow=RepN, ncol=p)
beta_aa4_mean=colMeans(beta_aa4)
deviat_beta_aa4=beta_aa4-matrix(beta_aa4_mean, byrow=T,nrow=RepN, ncol=p)


#SVAR_aa1=sum(deviat_beta_aa1*deviat_beta_aa1)/(RepN-1)
SVAR_aa2=sum(deviat_beta_aa2*deviat_beta_aa2)/(RepN-1)
SVAR_aa3=sum(deviat_beta_aa3*deviat_beta_aa3)/(RepN-1)
SVAR_aa4=sum(deviat_beta_aa4*deviat_beta_aa4)/(RepN-1)

beta_matrix=matrix(bet0, byrow=T, nrow=R, ncol=p)
beta=beta_matrix
#beta_aa1=beta_aa1-beta
beta_aa2=beta_aa2-beta
beta_aa3=beta_aa3-beta
beta_aa4=beta_aa4-beta

#EMSE_aa1=sum(beta_aa1*beta_aa1)/R
EMSE_aa2=sum(beta_aa2*beta_aa2)/R
EMSE_aa3=sum(beta_aa3*beta_aa3)/R
EMSE_aa4=sum(beta_aa4*beta_aa4)/R  #corrected on 04/20/24 by Y. Zuo 3 on RHS changed to 4 

print(c(R,n,p,N_ls,alpha,epsilon, Del, Kap) )
print( c(EMSE_aa2, EMSE_aa3, EMSE_aa4))
print(c(t_aa2,t_aa3,t_aa4))
print(EMSE_aa4/c(EMSE_aa2, EMSE_aa3,  EMSE_aa4))
print(SVAR_aa4/c(SVAR_aa2, SVAR_aa3,  SVAR_aa4))
print(c(SVAR_aa2, SVAR_aa3,  SVAR_aa4))

#SS_aa1=beta_aa1*beta_aa1
SS_aa2=beta_aa2*beta_aa2
SS_aa3=beta_aa3*beta_aa3
SS_aa4=beta_aa4*beta_aa4 #corrected on 04/20/24 by Y. Zuo 3 on RHS changed to 4 

#TT_aa1=t_aa1
TT_aa2=t_aa2
TT_aa3=t_aa3
TT_aa4=t_aa4

#RE_aa1=EMSE_aa4/EMSE_aa1
RE_aa2=EMSE_aa4/EMSE_aa2
RE_aa3=EMSE_aa4/EMSE_aa3
RE_aa4=EMSE_aa4/EMSE_aa4

par(mfrow = c(2,2))
boxplot(cbind(#EMSE_aa1[1:1000], 
  EMSE_aa2[1:1000], EMSE_aa3[1:1000], EMSE_aa4[1:1000]), 
  names=c("LST","LTS","LS"),
  xlab="procedure",
  ylab="EMSE",
  color=c("blue","green", "purple"),
  border="red"
)

boxplot(cbind(#SS_aa1[1:1000], 
  SS_aa2[1:1000], SS_aa3[1:1000], SS_aa4[1:1000]), 
  names=c(#"LST", 
    "LST","LTS", "LS"),
  xlab="procedure",
  ylab="SD",
  color=c("blue", "green", "purple"),
  border="blue"
)

boxplot(cbind(TT_aa2, TT_aa3, TT_aa4), 
    names=c("LST","LTS", "LS"),
  xlab="procedure",
  ylab="TT(seconds)",
  color=c("blue", "green", "purple"),
  border="green"
)
boxplot(cbind(#RE_aa1, 
  RE_aa2, RE_aa3, RE_aa4), 
  names=c(#"LST", 
    "LST","LTS", "LS"),
  xlab="procedure",
  ylab="RE",
  color=c("blue", "green", "purple"),
  border="purple"
)