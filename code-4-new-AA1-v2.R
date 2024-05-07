#code for new algorithm for LST based on Theorem 2.1 by Y. Zuo 10/27/22
install.packages("lm.beta")         # Install lm.beta package
library("lm.beta")  
library(robustbase)

AA1_main_new_lst_v2 = function(z, alpha, delta, N_ls)
{
  #z is the input data set, with n rows and p columns the last column is the yi
  #alpha, a tuning parameter, default is one, the larger, the faster but less robust
  #N_ls , a tuning parameter, the total number of local LS estimation decoded to perform
  #usually is set to be 500 or 1000
  #delta is the amount of small perturbation applied to beta on the boundary, default=.5
  
  #Initialization step
  z=data.matrix(z)
  n = dim(z)[1]
  p = dim(z)[2]
  N1 = n * (n - 1) / 2
  N2 = choose(n, floor((n + 1) / 2))
  N = min(N_ls, N2)
  M_combi = matrix(0, nrow = 2, ncol = N1)
  #M_index = matrix(1, nrow = N, ncol = n)
  M_bet = matrix(0, nrow = p, ncol = (4 * p + 2))
  temp = rep(0, (p - 1))
  bet_0 = bet_1 = P_zero =rep(0, p)    #matrix(0, nrow = p, ncol = 1)
  CR = 1 #current row number of M_index
  SS_min = 10 ^ 8 #Q initial value
  T_ls = 0
  bet_final = rep(0, p)
  
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
    #  M=length(temp1)
    m = temp1[1]
    #  for(m in 1:M)  #2nd for loop
    #   {
    l <- which(temp != 0)[m] # implicitly assume that xi is not identical to xj
    P_zero[(l + 1)] <- (yi - yj) / (xi[l] - xj[l])
    bet_0 <- P_zero
    bet_0[1] <- 1
    bet_1 <- bet_0
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
      rn = rep(0, n)
      I_bet = rep(0, n)
      sd = rep(1, n)
      bt = M_bet[, l]
      
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
        # if(is.element(indx_bet, M_index))
        # {break}
        # else
        # {
        #   M_index[CR,]<-indx_bet;
        CR = CR + 1
        fit <- lm(subdata[, p] ~ subdata[, 1:(p - 1)])
        T_ls = T_ls + 1
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
    #    } #2nd for end
    if (T_ls > N) { break }
  }#big for end
  
  return(bet_final)
} #function end