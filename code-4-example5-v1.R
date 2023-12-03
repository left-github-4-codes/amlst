# code 4 amstat by Y. Zuo on 11/26/23
###########################################part(a)###################
#rm(list=ls())
library(mvtnorm)
library(robustbase)
library(compiler)
enableJIT(1)

# library(AER)
# #library(MASS)
# #data(CPS1985)
# data(CPS1985)
# #attach(CPS1985)
# z=CPS1985
# c1=colnames(z)[2:11]
# ctemp=rep("wage",11)
# ctemp[1:10]<-c1
# z=z[ctemp] # switch the 1st colum of z to the last colum since wage is the reponse variable (y) 

# library(gamlss.data)
# data(plasma)
# z=plasma
# colnames(z)<-NULL #remove the column names
# #data(Boston, package="MASS")
# options(warn=-1)

###########################################################################
##### just entire data set and calulate the total time consumed and sample variance
###########################################################################
#z=Boston
PM10<-scan("PM10.txt")
pm10=matrix(PM10, byrow=T, nrow=500)
z=pm10

temp=z[,2:8]; col1=z[,1]
z=cbind(temp, col1)
colnames(z)<-NULL #remove the column names
n=dim(z)[1]; p=dim(z)[2]
# z=matrix(unlist(z), nrow=n, ncol=p)
# zz=z
# plasma=z

NO2<-scan("NO2.txt")
NO2=matrix(NO2, byrow=T, ncol=8)
z=NO2
tt=z[, 2:8]; t1=z[,1]
z=cbind(tt,t1)
n=dim(z)[1]; p=dim(z)[2]

RepN=1000; R=RepN
alpha=1 #h will automaticaly be set to be the default value: h=floor((n+p+1)/2) when p>2
c=0
cut_off=10^{-3}
epsilon=cut_off
N=200
gamma=10^2
beta_aa1=beta_aa2=beta_aa3=matrix(0, nrow=R, ncol=p) #aa3 is actually replace by LTS, aa1-LS, aa2-LST
t_aa1=t_aa2=t_aa3=0
#i=1
for (i in 1:RepN)
{
  t1=Sys.time()  
  fit2<-lm(z[,p]~z[,1:(p-1)])
  beta_aa1[i,]=as.numeric(fit2$coefficients) 
  #t1=Sys.time()  
  #beta_aa1[i,]=AA1_main_lts(z, alpha, c,N, cut_off )
  t2=Sys.time()-t1
  t_aa1=t_aa1+t2
  
  t1=Sys.time()  
  beta_aa2[i,]=lstReg(z, alpha, delta, N_ls)
  #beta_aa2[i,]=lstRegv1(z, alpha, delta, N_ls)
  t2=Sys.time()-t1
  t_aa2=t_aa2+t2
  
  t1=Sys.time()  
  fit1<-ltsReg(z[,1:(p-1)], z[,p])
  beta_aa3[i,]=as.numeric(fit1$coefficients) 
  t2=Sys.time()-t1
  t_aa3=t_aa3+t2
  
  print(i)
}
beta_aa1_mean=colMeans(beta_aa1)
deviat_beta_aa1=beta_aa1-matrix(beta_aa1_mean, byrow=T, nrow=RepN, ncol=p)
beta_aa2_mean=colMeans(beta_aa2)
deviat_beta_aa2=beta_aa2-matrix(beta_aa2_mean, byrow=T, nrow=RepN, ncol=p)
beta_aa3_mean=colMeans(beta_aa3)
deviat_beta_aa3=beta_aa3-matrix(beta_aa3_mean, byrow=T, nrow=RepN, ncol=p)

EMSE_aa1=sum(deviat_beta_aa1*deviat_beta_aa1)/RepN
EMSE_aa2=sum(deviat_beta_aa2*deviat_beta_aa2)/RepN
EMSE_aa3=sum(deviat_beta_aa3*deviat_beta_aa3)/RepN

print(c(R, n,p,c, N,alpha, gamma,epsilon,cut_off) )
print(c(t_aa1, t_aa2, t_aa3))
print(c(EMSE_aa1, EMSE_aa2, EMSE_aa3))
print(c(EMSE_aa1/c(EMSE_aa1, EMSE_aa2, EMSE_aa3)))



####################################part(b)#####################################
#rm(list=ls())
library(mvtnorm)
library(robustbase)
library(compiler)
enableJIT(1)
library(MASS)
library(gamlss.data)
# data(plasma)
# z=plasma
# colnames(z)<-NULL #remove the column names
# #data(Boston, package="MASS")
# options(warn=-1)

PM10<-scan("PM10.txt")
pm10=matrix(PM10, byrow=T, nrow=500)
z=pm10

temp=z[,2:8]; col1=z[,1]
z=cbind(temp, col1)
colnames(z)<-NULL #remove the column names
n=dim(z)[1]; p=dim(z)[2]
#z=pm10
zzz=z #zz used below and changes every time so reserve zz in zzz
z=zzz 
#z=Boston
#Bos_minus=zz[,c(-3,-7)] # delete indus and age two predictor since they are not significant with p-value <0.001
#z=Bos_minus

NO2<-scan("NO2.txt")
NO2=matrix(NO2, byrow=T, ncol=8)
z=NO2
tt=z[, 2:8]; t1=z[,1]
z=cbind(tt,t1)

n=dim(z)[1]; p=dim(z)[2]
RepN=1000; R=RepN
m=100  # tune this 200, 300, 400
alpha=5
N_ls=1
delta=1
# m=150  # tune this 200, 300, 400
# alpha=1/2 # one or three?
# c=0
# cut_off=10^{-3}
# N=100
# gamma=10^2
beta_aa1=beta_aa2=beta_aa3=matrix(0, nrow=R, ncol=p) #aa3 is actually replace by ltsReg, aa1--LS, aa2-LST
rss_aa1=rss_aa2=rss_aa3=rep(0, R)
res_aa1=res_aa2=res_aa3=matrix(0, nrow=n-m, ncol=1)
t_aa1=t_aa2=t_aa3=t1=t2=0
i=1
for (i in 1:RepN)
{
  index=sample(1:n, m)
  z=zzz[index,]; z=data.matrix(z)
  if(m==n) {zz=z} else{zz=zzz[-index,]} # could drop this and use all data points
  nn=dim(zz)[1]
  
  t1=Sys.time()  
  fit1<-lmrob(z[,p]~z[,1:(p-1)])
  #fit1<-rlm(z[,p]~z[,1:(p-1)])
  #fit0<-rreg(z[,1:(p-1)], z[,p])
  #  fit0<-lmsreg(z[,1:(p-1)], z[,p])
  # fit1<-ltsReg(z[,1:(p-1)], z[,p])
  beta_aa3[i,]=as.numeric(fit1$coefficients) 
  res_aa3=zz[,p]-data.matrix(cbind(matrix(1, nrow=nn), zz[, 1:(p-1)]))%*%matrix(beta_aa3[i,])
  
 # if(sum(beta_aa3[i,]!=0)) #ltsReg often output the initialized zero vector for beta
 #  {  
    rss_aa3[i]=sum(res_aa3*res_aa3)
    #rss_aa3=rss_aa3+rss_temp
    t2=Sys.time()-t1
    t_aa3=t_aa3+t2
    
    t1=Sys.time()  
    fit2<-lm(z[,p]~z[,1:(p-1)])
    beta_aa1[i,]=as.numeric(fit2$coefficients) 
    #beta_aa1[i,]=AA1_main_lts(z, alpha, c, N, cut_off )
    res_aa1=zz[,p]-data.matrix(cbind(matrix(1, nrow=nn), zz[, 1:(p-1)]))%*%matrix(beta_aa1[i,], nrow=p, ncol=1)
    
    rss_aa1[i]=sum(res_aa1*res_aa1)
    #rss_aa1=rss_aa1+rss_temp
    t2=Sys.time()-t1
    t_aa1=t_aa1+t2
    
    t1=Sys.time()  
    beta_aa2[i,]=lstReg(z, alpha, delta, N_ls)
    res_aa2=zz[,p]-data.matrix(cbind(matrix(1, nrow=nn), zz[, 1:(p-1)]))%*%matrix(beta_aa2[i,], nrow=p, ncol=1)
    
    rss_aa2[i]=sum(res_aa2*res_aa2)
    #rss_aa2=rss_aa2+rss_temp
    t2=Sys.time()-t1
    t_aa2=t_aa2+t2
#  } # end of if 
  # t1=Sys.time()  
  # fit1<-ltsReg(z[,1:(p-1)], z[,p])
  # beta_aa3[i,]=as.numeric(fit1$coefficients) 
  # res_aa3=zz[,p]-as.matrix(cbind(matrix(1, nrow=n-m), zz[, 1:(p-1)]))%*%matrix(beta_aa3[i,], nrow=p, ncol=1)
  # #res_aa3=Bos[,p]-as.matrix(cbind(matrix(1, nrow=n), Bos[, 1:(p-1)]))%*%matrix(beta_aa3[i,], nrow=p, ncol=1)
  
  print(i)
} # end of for loop

print(c(R, n,p, m,c, N,alpha, delta, N_ls) )
print(c(t_aa1, t_aa2, t_aa3))
print(c(sum(rss_aa1), sum(rss_aa2), sum(rss_aa3))/RepN)