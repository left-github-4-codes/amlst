#lst for am-stat by Y. Zuo on 11/22/2023

install.packages("lm.beta")         # Install lm.beta package
library("lm.beta")  
library(robustbase)
library(compiler)
enableJIT(1)

lstReg1 = function(z, alpha, delta, N_ls)
{
  #z is the input data set, with n rows and p columns the last column is the yi
  #alpha, a tuning parameter, default is one, the larger, the faster but less efficient
  #N_ls , a tuning parameter, the total number of local LS estimation decided to perform
  #usually is set to be one to ten, default=1
  #delta is the amount of small perturbation applied to beta on the boundary, default=.5
  
  #Initialization step
  z=data.matrix(z)
  n = dim(z)[1]  
  p = dim(z)[2]
  N1 = n * (n - 1) / 2
  N2 = choose(n, floor((n + 1) / 2))
  N = min(N_ls, N2)
  M_combi = matrix(0, nrow = 2, ncol = N1)
  M_bet = matrix(0, nrow = p, ncol = (4 * p + 2))
  temp = rep(0, (p - 1))
  bet_0 = bet_1 = P_zero = bet_final =rep(0, p)    #matrix(0, nrow = p, ncol = 1)
  CR = 1 #current row number of M_index
  SS_min = 10 ^ 8 #Q initial value
  
  
  M_combi = combn(1:n, 2)   # all combinations of two district indices
  
  for (k in 1:N1) #big for loop
  {
    i = M_combi[1, k]       
    j = M_combi[2, k]
    xi = z[i, 1:(p - 1)]    
    xj = z[j, 1:(p - 1)]
    yi = z[i, p]           
    yj = z[j, p]
    
    temp = xi - xj
    temp1 = which(temp != 0)
    l=m = temp1[1]

    P_zero[(l + 1)] <- (yi - yj) / (xi[l] - xj[l])
    bet_0 =bet_1 <- P_zero
    bet_0[1] <- 1
    bet_0=as.numeric(bet_0); bet_1=as.numeric(bet_1)
    #
    # create the 4pbeta matrix
    #
    for (j in 1:p)
    {
      v0 = v00 = bet_0
      v1 = v11 = bet_1
      
      
      v0[j] <- (v0[j] + delta)
      v1[j] <- (v1[j] + delta)
      M_bet[, j] <- v0          
      M_bet[, j + p] <- v1
      v00[j] <- v00[j] - delta
      v11[j] <- v11[j] - delta
      M_bet[, j + 2 * p] <- v00    
      M_bet[, j + 3 * p] <- v11
      
    }
    M_bet[, (4 * p + 1)] = bet_0
    M_bet[, (4 * p + 2)] = bet_1 #adding the two initial betas
    #
    # calculate local LS and update the beta_lst and the SS_min
    #
    for (l in 1:(4 * p + 2))
    {
      rn = I_bet= rep(0, n)
      sd = rep(1, n); bt = M_bet[, l]
      
      rn = z[, p] -data.matrix(cbind(rep(1, n), z[, 1:(p - 1)])) %*% as.matrix(bt)
      med = median(rn)
      madd = mad(rn)
      abdv = abs(rn - med)
      sd = abdv / madd
      
      indx_bet = which(sd <= alpha)
      K = length(indx_bet)
      D_seqn <- sd[indx_bet]
      subdata = z[indx_bet, ]
      
      if (length(unique(D_seqn)) != K) {
        break
      }
      else
      {
        CR = CR + 1
        fit <- lm(subdata[, p] ~ subdata[, 1:(p - 1)])
        SS_ls = sum(resid(fit) ^ 2)
        cof <- fit$coefficients
        names(cof) <- NULL
        bet_ls = cof
        
        if (SS_ls < SS_min)
        {
          SS_min = SS_ls
          bet_final = bet_ls
        }
      }#end of else
      
    }# end of (4p+2) loop
    if (CR > N) { break }
  }#big for N1 loop end
  
  return(bet_final)
} #function end
lstReg=cmpfun(lstReg1)
############################################################
n=80; p=2#
epsilon=0.3; n1=floor(epsilon*n)
m1=rmvnorm(n, mean=(rep(0, p)),sigma=matrix(c(1, 0.88, 0.88, 1), ncol=2, byrow=T))
m111=m1
m1=m111

par(mfrow=c(1,2))  
plot(m1[,1], m1[,2], xlim=c(-3,8), ylim=c(-2, 7), xlab="x-axis", ylab="y-axis", pch=20)
fit0<-ltsReg(m1[,1:(p-1)], m1[,p])
abline(fit0$coefficients[[1]], fit0$coefficients[[2]], col=1, lty=1, lwd=1)
#abline(AA1_main(m1,1, 0, 500, 0.001), col=2, lty=2, lwd=1)
fit2<-lm(m1[,2]~m1[,1]) #p is assume to be 2
abline(fit2, col=2, lty=2, lwd=1)
#wls_beta=wls(m1, 5, 6, 0.00001)
#abline(wls_beta, lty=3, col=3, lwd=1)
lst_beta=lstReg(m1, 1, 1, 1)
abline(lst_beta, lty=3, col=3, lwd=2)
legend(x = "topleft", legend=c("LTS", "LS", "LST"),    
       col=1:3, lty=1:3, cex=0.8,
       title="Line types", text.font=1, bg='lightblue')

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
#wls_beta=wls(m1, 5, 6, 0.00001)
#abline(wls_beta, lty=3, col=3, lwd=1)
lst_beta=lstReg(m1, 1, 1, 1)
abline(lst_beta, lty=3, col=3, lwd=2)
#text(4., -1.6, expression ('LTS'))
#text(6, 5.7, expression('LST'))
#text(4, -1, expression('LS'))

legend(x = "topleft", legend=c("LTS", "LS", "LST"),    
       col=1:3, lty=1:3, cex=0.8,
       title="Line types", text.font=1, bg='lightblue')
