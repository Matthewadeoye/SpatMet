#Gradients (Fisher identity), Hessian (Louis identity), loglikelihood
RemoveLatergradmultstrainLoglikelihood2<- function(y, e_it, nstrain, r, s, u, Gamma, B, Bits, a_k, Model, Q_r, Q_s, Q_u){
  ndept <- nrow(y[,,1])
  time  <- ncol(y[,,1])
  nstate <- 2^nstrain
  if(Model == 0){
    month_indexes <- ((1:time - 1) %% 12) + 1

    r_mat <- matrix(rep(r, each = ndept), nrow = ndept)
    s_mat <- matrix(rep(s[month_indexes], each = ndept), nrow = ndept)
    u_mat <- matrix(rep(u, times = time), nrow = ndept)

    log_risk <- r_mat + s_mat + u_mat

    poisMean<- matrix(0, nrow = ndept, ncol = time)
    allPoisMean<- array(NA, dim = c(ndept,time,nstrain))
    delta<- matrix(0, nrow = ndept, ncol = time)
    for(k in 1:nstrain){
      delta<- delta + (y[,,k] - (e_it * exp(log_risk + a_k[k])))
      poisMean<- poisMean + e_it * exp(log_risk + a_k[k])
      allPoisMean[,,k]<- e_it * exp(log_risk + a_k[k])
    }
    loglike<- sum(dpois(y, lambda = allPoisMean, log = T))

    # Temporal trend r
    grad_r <- colSums(delta) - as.numeric(Q_r %*% r)
    cov_r<- solve(diag(colSums(poisMean)) + Q_r + diag(1e-8, time))

    # Seasonal s
    fishervec_s<- numeric(12)
    grad_s <- numeric(12)
    for (month_index in 1:12) {
      t_idx <- which(((1:time - 1) %% 12 + 1) == month_index)
      grad_s[month_index] <- sum(delta[, t_idx])
      fishervec_s[month_index]<- sum(poisMean[, t_idx])
    }
    grad_s <- grad_s - as.numeric(Q_s %*% s)
    cov_s<- solve(diag(fishervec_s) + Q_s)

    # Spatial u
    grad_u <- rowSums(delta) - as.numeric(Q_u %*% u)

    poisMean4GibbsUpdate<- sum(e_it * exp(log_risk))

    return(list(loglike = loglike, grad_r = grad_r, grad_s = grad_s, grad_u = grad_u, cov_r=cov_r, cov_s=cov_s, poisMean4GibbsUpdate=poisMean4GibbsUpdate))
  }else{

    loglike_total <- 0

    JointTPM <- JointTransitionMatrix(gamma = Gamma, K = nstrain)
    logJointTPM <- log(JointTPM)

    E_lambda_tk <- array(0, dim = c(ndept,time,nstrain))   #Expected Poisson mean
    #Louis_integral_mat<- matrix(0, nrow = time, ncol = time) #smoothed_margin integral in Louis identity
    E_lambda_tk2 <- array(0, dim = c(ndept,time,nstrain))   #for Gibbs update of a_k's

    for(i in 1:ndept){

      logEmissions <- matrix(NA, nrow = time, ncol = nstate)
      lambda_array  <- array(0, dim = c(time, nstate, nstrain))
      #Louis_Grad_array  <- array(0, dim = c(time, nstate, nstrain)) #for Louise identity
      lambda_array2  <- array(0, dim = c(time, nstate, nstrain)) #for Gibbs update of a_k's

      for(t in 1:time){
        month_index <- (t-1) %% 12 + 1
        for(n in 1:nstate){
          for(k in 1:nstrain){
            newB<- rep(0, nstrain)
            newB[k]<- B[k]
            lambda_array[t,n,k] <- e_it[i,t] * exp(a_k[k] + r[t] + s[month_index] + u[i] + as.numeric(newB %*% Bits[n, ]))
           # Louis_Grad_array[t,n,k] <- y[i,t,k] - e_it[i,t] * exp(a_k[k] + r[t] + s[month_index] + u[i] + as.numeric(newB %*% Bits[n, ])) #for Louise identity
            lambda_array2[t,n,k] <- e_it[i,t] * exp(r[t] + s[month_index] + u[i] + as.numeric(newB %*% Bits[n, ]))  #for Gibbs update of a_k's
          }
          logEmissions[t,n] <- sum(dpois(y[i,t,], lambda = lambda_array[t,n,], log = TRUE))
        }
      }

      # forward pass
      loginit <- log(stationarydist(JointTPM))
      logalpha <- matrix(-Inf, nrow = time, ncol = nstate)
      logalpha[1, ] <- loginit + logEmissions[1, ]
      for(t in 2:time){
        logalpha[t, ] <- logspace_vecmatmult(logalpha[t-1, ], logJointTPM) + logEmissions[t, ]
      }

      loglik_i <- logSumExp_cpp(logalpha[time, ])
      loglike_total <- loglike_total + loglik_i

      # backward pass
      logbeta <- matrix(-Inf, nrow = time, ncol = nstate)
      logbeta[time, ] <- 0
      for(t in seq(time-1, 1, by = -1)){
        logbeta[t, ] <- logspace_vecmatmult(logEmissions[t+1, ] + logbeta[t+1, ], t(logJointTPM))
      }

      #Marginal posterior probabilities P_s
      logP_s<- (logalpha + logbeta) - loglik_i
      P_s<- exp(logP_s)

      #Louis grad_r
      #GradR_allstates<- matrix(0, nrow=time, ncol=nstate)
      #HessianR_allstates<- matrix(0, nrow=time, ncol=nstate)
      #for(k in 1:nstrain){
      #GradR_allstates<- GradR_allstates + Louis_Grad_array[,,k]
      #HessianR_allstates<- HessianR_allstates + lambda_array[,,k]
      #}
      #HessianRmat_allstates<- array(NA, dim=c(time,time,nstate))

      #for(n in 1:nstate){
      #  HessianRmat_allstates[,,n] <- P_s[,n] * (diag(HessianR_allstates[, n]) + GradR_allstates[,n] %*% t(GradR_allstates[,n]))
      #  Louis_integral_mat<- Louis_integral_mat + HessianRmat_allstates[,,n]
      #}

      #Expected Poisson mean
      for(t in 1:time){
        for(k in 1:nstrain){
          E_lambda_tk[i,t,k] <- sum(P_s[t, ] * lambda_array[t, , k])
          E_lambda_tk2[i,t,k] <- sum(P_s[t, ] * lambda_array2[t, , k])  #for Gibbs update of a_k's
        }
      }
    }

    poisMean4GibbsUpdate<- numeric(nstrain)
    poisMean<- matrix(0, nrow = ndept, ncol = time)
    delta<- matrix(0, nrow = ndept, ncol = time)
    for(k in 1:nstrain){
      delta<- delta + (y[,,k] - E_lambda_tk[,,k])
      poisMean<- poisMean + E_lambda_tk[,,k]
      poisMean4GibbsUpdate[k]<- sum(E_lambda_tk2[,,k])  #for Gibbs update of a_k's
    }

    # Temporal trend r
    grad_r <- colSums(delta) - as.numeric(Q_r %*% r)
    cov_r<- solve(diag(colSums(poisMean)) + Q_r + diag(1e-8, time))
    #newcov_r<- solve(Louis_integral_mat - grad_r %*% t(grad_r) + Q_r + diag(1e-8, time))

    # Seasonal s
    fishervec_s<- numeric(12)
    grad_s <- numeric(12)
    newcov_s<- matrix(NA, 12, 12)
    for (month_index in 1:12) {
      t_idx <- which(((1:time - 1) %% 12 + 1) == month_index)
      grad_s[month_index] <- sum(delta[, t_idx])
      fishervec_s[month_index]<- sum(poisMean[, t_idx])
    }
    grad_s <- grad_s - as.numeric(Q_s %*% s)
    cov_s<- solve(diag(fishervec_s) + Q_s)

    #month_indexes <- ((1:time - 1) %% 12) + 1
    #newcov_s<- solve(fishermat_s - grad_s%*%t(grad_s) + Q_s)

    # Spatial u
    grad_u <- rowSums(delta) - as.numeric(Q_u %*% u)

    return(list(loglike = loglike_total, grad_r = grad_r, grad_s = grad_s, grad_u = grad_u, cov_r=cov_r, cov_s=cov_s, poisMean4GibbsUpdate=poisMean4GibbsUpdate, newcov_r=cov_r))
  }
}

#Gradients, Hessian, loglikelihood --- FFBS
gradmultstrainLoglikelihood2<- function(y, e_it, nstrain, r, s, u, Gamma, B, Bits, a_k, Model, Q_r, Q_s, Q_u){
  ndept <- nrow(y[,,1])
  time  <- ncol(y[,,1])
  nstate <- 2^nstrain
  if(Model == 0){
    month_indexes <- ((1:time - 1) %% 12) + 1

    r_mat <- matrix(rep(r, each = ndept), nrow = ndept)
    s_mat <- matrix(rep(s[month_indexes], each = ndept), nrow = ndept)
    u_mat <- matrix(rep(u, times = time), nrow = ndept)

    log_risk <- r_mat + s_mat + u_mat

    poisMean<- matrix(0, nrow = ndept, ncol = time)
    allPoisMean<- array(NA, dim = c(ndept,time,nstrain))
    delta<- matrix(0, nrow = ndept, ncol = time)
    for(k in 1:nstrain){
      delta<- delta + (y[,,k] - (e_it * exp(log_risk + a_k[k])))
      poisMean<- poisMean + e_it * exp(log_risk + a_k[k])
      allPoisMean[,,k]<- e_it * exp(log_risk + a_k[k])
    }
    loglike<- sum(dpois(y, lambda = allPoisMean, log = T))

    # Temporal trend r
    grad_r <- colSums(delta) - as.numeric(Q_r %*% r)
    cov_r<- solve(diag(colSums(poisMean)) + Q_r + diag(1e-8, time))

    # Seasonal s
    fishervec_s<- numeric(12)
    grad_s <- numeric(12)
    for (month_index in 1:12) {
      t_idx <- which(((1:time - 1) %% 12 + 1) == month_index)
      grad_s[month_index] <- sum(delta[, t_idx])
      fishervec_s[month_index]<- sum(poisMean[, t_idx])
    }
    grad_s <- grad_s - as.numeric(Q_s %*% s)
    cov_s<- solve(diag(fishervec_s) + Q_s)

    # Spatial u
    grad_u <- rowSums(delta) - as.numeric(Q_u %*% u)

    poisMean4GibbsUpdate<- sum(e_it * exp(log_risk))

    return(list(loglike = loglike, grad_r = grad_r, grad_s = grad_s, grad_u = grad_u, cov_r=cov_r, cov_s=cov_s, poisMean4GibbsUpdate=poisMean4GibbsUpdate))
  }else{

    loglike_total <- 0
    Actual_likelihood<- 0

    JointTPM <- JointTransitionMatrix(gamma = Gamma, K = nstrain)
    logJointTPM <- log(JointTPM)

    Actual_lambda_tk <- array(0, dim = c(ndept,time,nstrain))   #Actual Poisson mean
    Actual_lambda_tk2 <- array(0, dim = c(ndept,time,nstrain))   #for Gibbs update of a_k's

    for(i in 1:ndept){

      logEmissions <- matrix(NA, nrow = time, ncol = nstate)
      lambda_array  <- array(0, dim = c(time, nstate, nstrain))
      lambda_array2  <- array(0, dim = c(time, nstate, nstrain)) #for Gibbs update of a_k's

      for(t in 1:time){
        month_index <- (t-1) %% 12 + 1
        for(n in 1:nstate){
          for(k in 1:nstrain){
            newB<- rep(0, nstrain)
            newB[k]<- B[k]
            lambda_array[t,n,k] <- e_it[i,t] * exp(a_k[k] + r[t] + s[month_index] + u[i] + as.numeric(newB %*% Bits[n, ]))
            lambda_array2[t,n,k] <- e_it[i,t] * exp(r[t] + s[month_index] + u[i] + as.numeric(newB %*% Bits[n, ]))  #for Gibbs update of a_k's
          }
          logEmissions[t,n] <- sum(dpois(y[i,t,], lambda = lambda_array[t,n,], log = TRUE))
        }
      }

      # forward filtering
      loginit <- log(stationarydist(JointTPM))
      logalpha <- matrix(-Inf, nrow = time, ncol = nstate)
      logalpha[1, ] <- loginit + logEmissions[1, ]
      for(t in 2:time){
        logalpha[t, ] <- logspace_vecmatmult(logalpha[t-1, ], logJointTPM) + logEmissions[t, ]
      }
      loglik_i <- logSumExp_cpp(logalpha[time, ])
      loglike_total <- loglike_total + loglik_i

      #Backward sampling for t =T,T-1,...,1
      states <- numeric(time)
      #final state P(s_T | observations)
      logP_T <- logalpha[time, ] - loglik_i
      P_T <- exp(logP_T - logSumExp_cpp(logP_T))
      states[time] <- sample(1:nstate, 1, prob = P_T)

      for (t in (time-1):1) {
        logP_t <- logalpha[t, ] + logJointTPM[, states[t+1]]
        P_t <- exp(logP_t - logSumExp_cpp(logP_t))
        states[t] <- sample(1:nstate, 1, prob = P_t)
      }

      #Actual Poisson mean
      for(t in 1:time){
        month_index <- (t-1) %% 12 + 1
        for(k in 1:nstrain){
          newB<- rep(0, nstrain)
          newB[k]<- B[k]
          Actual_lambda_tk[i,t,k] <- e_it[i,t] * exp(a_k[k] + r[t] + s[month_index] + u[i] + as.numeric(newB %*% Bits[states[t], ]))
          Actual_lambda_tk2[i,t,k] <- e_it[i,t] * exp(r[t] + s[month_index] + u[i] + as.numeric(newB %*% Bits[states[t], ]))  #for Gibbs update of a_k's
        }
      }
     # Actual_likelihood<- Actual_likelihood + sum(dpois(y[i,,], lambda = Actual_lambda_tk[i,,], log = TRUE))
    }

    poisMean4GibbsUpdate<- numeric(nstrain)
    poisMean<- matrix(0, nrow = ndept, ncol = time)
    delta<- matrix(0, nrow = ndept, ncol = time)
    for(k in 1:nstrain){
      delta<- delta + (y[,,k] - Actual_lambda_tk[,,k])
      poisMean<- poisMean + Actual_lambda_tk[,,k]
      poisMean4GibbsUpdate[k]<- sum(Actual_lambda_tk2[,,k])  #for Gibbs update of a_k's
    }

    # Temporal trend r
    grad_r <- colSums(delta) - as.numeric(Q_r %*% r)
    cov_r<- solve(diag(colSums(poisMean)) + Q_r + diag(1e-8, time))

    # Seasonal s
    fishervec_s<- numeric(12)
    grad_s <- numeric(12)
    for (month_index in 1:12) {
      t_idx <- which(((1:time - 1) %% 12 + 1) == month_index)
      grad_s[month_index] <- sum(delta[, t_idx])
      fishervec_s[month_index]<- sum(poisMean[, t_idx])
    }
    grad_s <- grad_s - as.numeric(Q_s %*% s)
    cov_s<- solve(diag(fishervec_s) + Q_s)

    # Spatial u
    grad_u <- rowSums(delta) - as.numeric(Q_u %*% u)

    return(list(loglike = loglike_total, grad_r = grad_r, grad_s = grad_s, grad_u = grad_u, cov_r=cov_r, cov_s=cov_s, poisMean4GibbsUpdate=poisMean4GibbsUpdate, newcov_r=cov_r))
  }
}

#Riemann Manifold Langevin updates
multMMALAInference<- function(y, e_it, Model, adjmat, step_sizes, num_iteration = 15000, sdGs=0.1, sdBs=0.03, sdAs=0.03) {
  start_time <- Sys.time()
  ndept <- nrow(e_it)
  time <- ncol(e_it)
  nstrain<- dim(y)[3]
  Bits<- encodeBits(nstrain)

  R<- -1 * adjmat
  diag(R)<- -rowSums(R, na.rm = T)
  rankdef<- nrow(R)-qr(R)$rank

  RW1PrecMat<- matrix(0, nrow=12, ncol=12)
  RW1PrecMat[1, ]<- c(2,-1, rep(0, 12-3), -1)
  RW1PrecMat[2, ]<- c(-1,2,-1, rep(0, 12-3))
  RW1PrecMat[3, ]<- c(0, -1,2,-1, rep(0, 12-4))
  RW1PrecMat[(12-1), ]<- c(rep(0, 12-3), -1,2,-1)
  RW1PrecMat[12, ]<- c(-1, rep(0, 12-3), -1, 2)
  for(i in 3:(12-3)){
    RW1PrecMat[i+1, ((i):(i+2))]<- c(-1,2,-1)
  }

  RW2PrecMat<- matrix(0, nrow=time, ncol=time)
  RW2PrecMat[1,(1:3)]<- c(1,-2,1)
  RW2PrecMat[2,(1:4)]<- c(-2,5,-4,1)
  RW2PrecMat[3,(1:5)]<- c(1,-4,6,-4,1)
  RW2PrecMat[(time-1),((time-3):time)]<- c(1,-4,5,-2)
  RW2PrecMat[time,((time-2):time)]<- c(1,-2,1)
  for(i in 3:(time-3)){
    RW2PrecMat[i+1, ((i-1):(i+3))]<- c(1,-4,6,-4,1)
  }

  SumYk_vec<- numeric(nstrain)
  SumYk_vec[1]<- sum(y[,,1])
  sumY<- y[,,1]
  for(k in 2:nstrain){
    sumY<- sumY + y[,,k]
    SumYk_vec[k]<- sum(y[,,k])
  }

  crudeResults<- DetectOutbreaks:::crudeEst(sumY, e_it)
  crudeR<- crudeResults[[1]] - mean(crudeResults[[1]])
  crudeS<- crudeResults[[2]]
  crudeU<- crudeResults[[3]]
  crudeU<- ifelse(is.nan(crudeU), mean(crudeU[is.finite(crudeU)]), crudeU)
  crudeblock<- floor(time/12)
  crudeblock<- ((crudeblock*12)-11):(crudeblock*12)

  MC_chain<- matrix(NA, nrow=num_iteration, ncol=5+time+12+ndept+nstrain+nstrain+1)
  initG12<- runif(1)
  initG21<- runif(1)
  initstateD<- stationarydist(G(initG12, initG21))[2]
  MC_chain[1,]<- c(initG12, initG21, 1/var(crudeR), 1/var(crudeS), 1/var(crudeU), crudeR, crudeS[crudeblock-12], crudeU, rep(0, nstrain), rep(mean(crudeResults[[1]]), nstrain), initstateD)

  #flag missing data for inference
  y<- ifelse(is.na(y), -1, y)
  yflat<- as.numeric(aperm(y, c(2,1,3)))

  Q_r<- MC_chain[1,3] * RW2PrecMat
  Q_s<- MC_chain[1,4] * RW1PrecMat
  Q_u<- MC_chain[1,5] * R

  #Compute gradients
  Allquantities<- gradmultstrainLoglikelihood2(y=y, e_it=e_it, nstrain=nstrain,  r=MC_chain[1, 5+(1:time)], s=MC_chain[1, 5+time+(1:12)], u=MC_chain[1, 5+time+12+(1:ndept)], Gamma=G(MC_chain[1,1],MC_chain[1,2]), B=MC_chain[1, 5+time+12+ndept+(1:nstrain)], Bits=Bits, a_k=MC_chain[1, 5+time+12+ndept+nstrain+(1:nstrain)], Model=Model,Q_r=Q_r,Q_s = Q_s,Q_u=Q_u)
  likelihoodcurrent<- Allquantities$loglike
  priorcurrentRcomps<- randomwalk2(MC_chain[1, 5+(1:time)], MC_chain[1, 3])
  priorcurrentScomps<- seasonalComp2(MC_chain[1, 5+time+(1:12)], MC_chain[1, 4], RW1PrecMat)
  priorcurrentUcomps<- logIGMRF1(MC_chain[1, 5+time+12+(1:ndept)], MC_chain[1, 5], R, rankdef)

  grad_current <- list(grad_r=Allquantities$grad_r, grad_s=Allquantities$grad_s, grad_u=Allquantities$grad_u, cov_r=Allquantities$cov_r, cov_s=Allquantities$cov_s)

  for (i in 2:num_iteration) {

    MC_chain[i,3]<- rgamma(1, shape = 1 + (time-2)/2, rate = 0.0001 + (t(MC_chain[i-1, 5+(1:time)]) %*% RW2PrecMat %*% MC_chain[i-1, 5+(1:time)])/2)
    MC_chain[i,4]<- rgamma(1, shape = 1 + 11/2, rate = 0.001 + (t(MC_chain[i-1, 5+time+(1:12)]) %*% RW1PrecMat %*% MC_chain[i-1, 5+time+(1:12)])/2)
    MC_chain[i,5]<- rgamma(1, shape = 1 + (ndept-1)/2, rate = 0.01 + (t(MC_chain[i-1, 5+time+12+(1:ndept)]) %*% R %*% MC_chain[i-1, 5+time+12+(1:ndept)])/2)

    Q_r<- MC_chain[i,3] * RW2PrecMat
    Q_s<- MC_chain[i,4] * RW1PrecMat
    Q_u<- MC_chain[i,5] * R

    current_r <- MC_chain[i-1, 5+(1:time)]
    current_s <- MC_chain[i-1, 5+time+(1:12)]
    current_u <- MC_chain[i-1, 5+time+12+(1:ndept)]

    #Update s
    Mmatcs<- as.numeric(grad_current$cov_s %*% grad_current$grad_s)

    eps_s <- rnorm(12)
    proposedScomps <- as.numeric(MC_chain[i-1, 5+time+(1:12)] + 0.5 * step_sizes$s^2 * Mmatcs + step_sizes$s * chol(grad_current$cov_s) %*% eps_s)
    proposedScomps<- proposedScomps - mean(proposedScomps)

    Allquantities<- gradmultstrainLoglikelihood2(y=y, e_it=e_it, nstrain=nstrain,  r=current_r, s=proposedScomps, u=current_u, Gamma=G(MC_chain[i-1,1],MC_chain[i-1,2]), B=MC_chain[i-1, 5+time+12+ndept+(1:nstrain)], Bits=Bits, a_k=MC_chain[i-1, 5+time+12+ndept+nstrain+(1:nstrain)], Model=Model,Q_r=Q_r,Q_s = Q_s,Q_u=Q_u)
    grad_proposed <- list(grad_r=Allquantities$grad_r, grad_s=Allquantities$grad_s, grad_u=Allquantities$grad_u, cov_r=Allquantities$cov_r, cov_s=Allquantities$cov_s)

    Mmatps<- as.numeric(grad_proposed$cov_s %*% grad_proposed$grad_s)

    likelihoodproposed<- Allquantities$loglike
    q_prop <- mvnfast::dmvn(proposedScomps, mu = MC_chain[i-1, 5+time+(1:12)] + 0.5 * step_sizes$s^2 * Mmatcs, sigma = grad_current$cov_s * step_sizes$s^2, log = TRUE)
    q_curr <- mvnfast::dmvn(MC_chain[i-1, 5+time+(1:12)], mu = proposedScomps + 0.5 * step_sizes$s^2 * Mmatps, sigma = grad_proposed$cov_s * step_sizes$s^2, log = TRUE)
    priorproposedScomps<- seasonalComp2(proposedScomps, MC_chain[i, 4], RW1PrecMat)

    log_alpha_s <- likelihoodproposed + priorproposedScomps + q_curr - likelihoodcurrent - priorcurrentScomps - q_prop
    if (is.finite(log_alpha_s) && log(runif(1)) < log_alpha_s){
      MC_chain[i, 5+time+(1:12)]<- proposedScomps
      likelihoodcurrent<- likelihoodproposed
      priorcurrentScomps<- priorproposedScomps
      grad_current<- grad_proposed
    }else{
      MC_chain[i, 5+time+(1:12)]<- MC_chain[i-1, 5+time+(1:12)]
    }

    #Update r
    eps_r <- rnorm(time)
    Mmatrc<- as.numeric(grad_current$cov_r %*% grad_current$grad_r)

    proposedRcomps <- as.numeric(current_r + 0.5 * step_sizes$r^2 * Mmatrc + step_sizes$r * chol(grad_current$cov_r) %*% eps_r)
    proposedRcomps<- proposedRcomps - mean(proposedRcomps)

    Allquantities<- gradmultstrainLoglikelihood2(y=y, e_it=e_it, nstrain=nstrain,  r=proposedRcomps, s=MC_chain[i, 5+time+(1:12)], u=current_u, Gamma=G(MC_chain[i-1,1],MC_chain[i-1,2]), B=MC_chain[i-1, 5+time+12+ndept+(1:nstrain)], Bits=Bits, a_k=MC_chain[i-1, 5+time+12+ndept+nstrain+(1:nstrain)], Model=Model,Q_r=Q_r,Q_s = Q_s,Q_u=Q_u)
    grad_proposed <- list(grad_r=Allquantities$grad_r, grad_s=Allquantities$grad_s, grad_u=Allquantities$grad_u, cov_r=Allquantities$cov_r, cov_s=Allquantities$cov_s)

    Mmatrp<- as.numeric(grad_proposed$cov_r %*% grad_proposed$grad_r)

    q_prop <- mvnfast::dmvn(proposedRcomps, mu = current_r + 0.5 * step_sizes$r^2 * Mmatrc, sigma = grad_current$cov_r * step_sizes$r^2, log = TRUE)
    q_curr <- mvnfast::dmvn(current_r, mu = proposedRcomps + 0.5 * step_sizes$r^2 * Mmatrp, sigma = grad_proposed$cov_r * step_sizes$r^2, log = TRUE)

    likelihoodproposed<- Allquantities$loglike
    priorproposedRcomps <- randomwalk2(proposedRcomps, MC_chain[i, 3])

    log_alpha_r <- likelihoodproposed + priorproposedRcomps + q_curr - likelihoodcurrent - priorcurrentRcomps - q_prop

    if (is.finite(log_alpha_r) && log(runif(1)) < log_alpha_r){
      MC_chain[i, 5+(1:time)] <- proposedRcomps
      likelihoodcurrent<- likelihoodproposed
      priorcurrentRcomps<- priorproposedRcomps
      grad_current<- grad_proposed
    }else{
      MC_chain[i, 5+(1:time)]<- MC_chain[i-1, 5+(1:time)]
    }

    #Update u
    eps_u <- rnorm(ndept)

    proposedUcomps <- as.numeric(MC_chain[i-1, 5+time+12+(1:ndept)] + 0.5 * step_sizes$u^2 * grad_current$grad_u + step_sizes$u * eps_u)
    proposedUcomps<- proposedUcomps - mean(proposedUcomps)

    Allquantities<- gradmultstrainLoglikelihood2(y=y, e_it=e_it, nstrain=nstrain,  r=MC_chain[i, 5+(1:time)], s=MC_chain[i, 5+time+(1:12)], u=proposedUcomps, Gamma=G(MC_chain[i-1,1],MC_chain[i-1,2]), B=MC_chain[i-1, 5+time+12+ndept+(1:nstrain)], Bits=Bits, a_k=MC_chain[i-1, 5+time+12+ndept+nstrain+(1:nstrain)], Model=Model,Q_r=Q_r,Q_s = Q_s,Q_u=Q_u)
    grad_proposed <- list(grad_r=Allquantities$grad_r, grad_s=Allquantities$grad_s, grad_u=Allquantities$grad_u, cov_r=Allquantities$cov_r, cov_s=Allquantities$cov_s)

    likelihoodproposed<- Allquantities$loglike

    q_prop <- sum(dnorm(proposedUcomps, mean = MC_chain[i-1, 5+time+12+(1:ndept)] + 0.5 * step_sizes$u^2 * grad_current$grad_u, sd = step_sizes$u, log = TRUE))
    q_curr <- sum(dnorm(MC_chain[i-1, 5+time+12+(1:ndept)], mean = proposedUcomps + 0.5 * step_sizes$u^2 * grad_proposed$grad_u, sd = step_sizes$u, log = TRUE))

    priorproposedUcomps<- logIGMRF1(proposedUcomps, MC_chain[i, 5], R, rankdef)
    priorcurrentUcomps<- logIGMRF1(MC_chain[i-1, 5+time+12+(1:ndept)], MC_chain[i, 5], R, rankdef)

    log_alpha_u <- likelihoodproposed + priorproposedUcomps + q_curr - likelihoodcurrent - priorcurrentUcomps - q_prop
    if (is.finite(log_alpha_u) && log(runif(1)) < log_alpha_u){
      MC_chain[i, 5+time+12+(1:ndept)]<- proposedUcomps
      likelihoodcurrent<- likelihoodproposed
      #priorcurrentUcomps<- priorproposedUcomps
      grad_current<- grad_proposed
    }else{
      MC_chain[i, 5+time+12+(1:ndept)]<- MC_chain[i-1, 5+time+12+(1:ndept)]
    }

    if(Model == 0){
      MC_chain[i, 5+time+12+ndept+(1:nstrain)]<-  MC_chain[i-1, 5+time+12+ndept+(1:nstrain)]
      MC_chain[i, 1:2]<- MC_chain[i-1, 1:2]
    }else{
      proposedB <- abs(rnorm(nstrain, mean = MC_chain[i-1, 5+time+12+ndept+(1:nstrain)], sd = rep(sdBs, nstrain)))
      priorcurrentB<- sum(dgamma(MC_chain[i-1, 5+time+12+ndept+(1:nstrain)], shape = rep(2, nstrain), rate = rep(2,nstrain), log=TRUE))
      priorproposedB<- sum(dgamma(proposedB, shape = rep(2, nstrain), rate = rep(2, nstrain), log=TRUE))

    Allquantities<- gradmultstrainLoglikelihood2(y=y, e_it=e_it, nstrain=nstrain,  r=MC_chain[i, 5+(1:time)], s=MC_chain[i, 5+time+(1:12)], u=MC_chain[i, 5+time+12+(1:ndept)], Gamma=G(MC_chain[i-1,1],MC_chain[i-1,2]), B=proposedB, Bits=Bits, a_k=MC_chain[i-1, 5+time+12+ndept+nstrain+(1:nstrain)], Model=Model,Q_r=Q_r,Q_s = Q_s,Q_u=Q_u)
    grad_proposed <- list(grad_r=Allquantities$grad_r, grad_s=Allquantities$grad_s, grad_u=Allquantities$grad_u, cov_r=Allquantities$cov_r, cov_s=Allquantities$cov_s)

    likelihoodproposed<- Allquantities$loglike

    mh.ratio<- exp(likelihoodproposed + priorproposedB
                   - likelihoodcurrent - priorcurrentB)

    #print(paste("mh.ratioB = ", mh.ratio))

    if(!is.na(mh.ratio) && runif(1) < mh.ratio){
      MC_chain[i, 5+time+12+ndept+(1:nstrain)]<- proposedB
      likelihoodcurrent<- likelihoodproposed
      grad_current<- grad_proposed
    }
    else{
      MC_chain[i, 5+time+12+ndept+(1:nstrain)]<- MC_chain[i-1, 5+time+12+ndept+(1:nstrain)]
    }

    proposedGs<- abs(rnorm(2,mean=c(MC_chain[i-1,1], MC_chain[i-1, 2]), sd=c(sdGs, sdGs)))
    if(proposedGs[1]>1) proposedGs[1]=2-proposedGs[1]
    if(proposedGs[2]>1) proposedGs[2]=2-proposedGs[2]

    priorcurrentGs<- sum(dbeta(MC_chain[i-1,1:2], shape1 = c(2,2), shape2 = c(2,2), log=TRUE))
    priorproposedGs<- sum(dbeta(proposedGs, shape1 = c(2,2), shape2 = c(2,2), log=TRUE))

    Allquantities<- gradmultstrainLoglikelihood2(y=y, e_it=e_it, nstrain=nstrain,  r=MC_chain[i, 5+(1:time)], s=MC_chain[i, 5+time+(1:12)], u=MC_chain[i, 5+time+12+(1:ndept)], Gamma=G(proposedGs[1],proposedGs[2]), B=MC_chain[i, 5+time+12+ndept+(1:nstrain)], Bits=Bits, a_k=MC_chain[i-1, 5+time+12+ndept+nstrain+(1:nstrain)], Model=Model,Q_r=Q_r,Q_s = Q_s,Q_u=Q_u)
    grad_proposed <- list(grad_r=Allquantities$grad_r, grad_s=Allquantities$grad_s, grad_u=Allquantities$grad_u, cov_r=Allquantities$cov_r, cov_s=Allquantities$cov_s)

    likelihoodproposed<- Allquantities$loglike

    mh.ratio<- exp(likelihoodproposed + priorproposedGs
                   - likelihoodcurrent - priorcurrentGs)

    #print(mh.ratio)

    if(!is.na(mh.ratio) && runif(1) < mh.ratio){
      MC_chain[i, 1:2]<- proposedGs
      likelihoodcurrent<- likelihoodproposed
      grad_current<- grad_proposed
    }
    else{
      MC_chain[i, 1:2]<- MC_chain[i-1,1:2]
    }
  }
    MC_chain[i, 5+time+12+ndept+nstrain+nstrain+1]<- stationarydist(G(MC_chain[i, 1], MC_chain[i, 2]))[2]

    #Gibbs A_k's update
    MC_chain[i, 5+time+12+ndept+nstrain+(1:nstrain)]<- log(rgamma(nstrain, shape = 0.01+SumYk_vec, rate = Allquantities$poisMean4GibbsUpdate + 0.01/exp(-15)))

    if(i %% 1000 == 0) cat("Iteration:", i, "\n")
  }
  colnames(MC_chain) <- paste(c("G12", "G21", "kappa_r", "kappa_s", "kappa_u", paste("r", 1:time, sep=""), paste("s", 1:12, sep=""), paste("u", 1:ndept, sep=""), paste("B", 1:nstrain, sep=""), paste("a_k", 1:nstrain, sep=""), "StationaryDistribution"))
  MC_chain<- as.data.frame(MC_chain)
  end_time <- Sys.time()
  time_taken<- end_time - start_time
  print(time_taken)
  return(MC_chain)
}


#set.seed(212);multmod0nstrain5<- Multstrain.simulate(Model = 0, time=60, adj.matrix = sim_adjmat, nstrain=5, B=c(1.65,0.95,1.4,1.1,1.7))
#set.seed(212);multmod1nstrain5<- Multstrain.simulate(Model = 1, time=60, adj.matrix = sim_adjmat, Modeltype = 1, nstrain=5, B=c(1.65,0.95,1.4,1.1,1.7))
#set.seed(212);perstrainmultmod1nstrain5<- Multstrain.simulate(Model = 1, time=60, adj.matrix = sim_adjmat, Modeltype = 2, nstrain=5, B=c(1.65,0.95,1.4,1.1,1.7))
#set.seed(212);dependentmultmod1nstrain5<- Multstrain.simulate(Model = 1, time=60, adj.matrix = sim_adjmat, Modeltype = 3, nstrain=5, B=c(1.65,0.95,1.4,1.1,1.7))

#multMmalaRes2<- CPPmultMMALAInference(y=multmod0nstrain5[[1]], e_it = multmod0nstrain5[[2]], Model = 0, adjmat = sim_adjmat, step_sizes = list("r"=0.3,"s"=0.3,"u"=0.025), num_iteration = 20000)
#gradmultstrainLoglikelihood2(y=multmod0nstrain5[["y"]], e_it=multmod0nstrain5[["e_it"]], nstrain=5, r=multmod0nstrain5[["r"]], s=multmod0nstrain5[["s"]], u=multmod0nstrain5[["u"]], Gamma=G(0.1,0.2), B=rep(0,5), Bits=Bits, a_k=multmod0nstrain5[["a_k"]], Model=0, Q_r=RW2PrecMat, Q_s=RW1PrecMat, Q_u=R)$loglike


#Riemann Manifold Langevin updates
CPPmultMMALAInference<- function(y, e_it, Model, adjmat, step_sizes, num_iteration = 15000, sdGs=0.1, sdBs=0.03, sdAs=0.03) {
  start_time <- Sys.time()
  ndept <- nrow(e_it)
  time <- ncol(e_it)
  nstrain<- dim(y)[3]
  Bits<- encodeBits(nstrain)

  R<- -1 * adjmat
  diag(R)<- -rowSums(R, na.rm = T)
  rankdef<- nrow(R)-qr(R)$rank

  RW1PrecMat<- matrix(0, nrow=12, ncol=12)
  RW1PrecMat[1, ]<- c(2,-1, rep(0, 12-3), -1)
  RW1PrecMat[2, ]<- c(-1,2,-1, rep(0, 12-3))
  RW1PrecMat[3, ]<- c(0, -1,2,-1, rep(0, 12-4))
  RW1PrecMat[(12-1), ]<- c(rep(0, 12-3), -1,2,-1)
  RW1PrecMat[12, ]<- c(-1, rep(0, 12-3), -1, 2)
  for(i in 3:(12-3)){
    RW1PrecMat[i+1, ((i):(i+2))]<- c(-1,2,-1)
  }

  RW2PrecMat<- matrix(0, nrow=time, ncol=time)
  RW2PrecMat[1,(1:3)]<- c(1,-2,1)
  RW2PrecMat[2,(1:4)]<- c(-2,5,-4,1)
  RW2PrecMat[3,(1:5)]<- c(1,-4,6,-4,1)
  RW2PrecMat[(time-1),((time-3):time)]<- c(1,-4,5,-2)
  RW2PrecMat[time,((time-2):time)]<- c(1,-2,1)
  for(i in 3:(time-3)){
    RW2PrecMat[i+1, ((i-1):(i+3))]<- c(1,-4,6,-4,1)
  }

  SumYk_vec<- numeric(nstrain)
  SumYk_vec[1]<- sum(y[,,1])
  sumY<- y[,,1]
  for(k in 2:nstrain){
    sumY<- sumY + y[,,k]
    SumYk_vec[k]<- sum(y[,,k])
  }

  crudeResults<- DetectOutbreaks:::crudeEst(sumY, e_it)
  crudeR<- crudeResults[[1]] - mean(crudeResults[[1]])
  crudeS<- crudeResults[[2]]
  crudeU<- crudeResults[[3]]
  crudeU<- ifelse(is.nan(crudeU), mean(crudeU[is.finite(crudeU)]), crudeU)
  crudeblock<- floor(time/12)
  crudeblock<- ((crudeblock*12)-11):(crudeblock*12)

  MC_chain<- matrix(NA, nrow=num_iteration, ncol=5+time+12+ndept+nstrain+nstrain+1)
  initG12<- runif(1)
  initG21<- runif(1)
  initstateD<- stationarydist(G(initG12, initG21))[2]
  MC_chain[1,]<- c(initG12, initG21, 1/var(crudeR), 1/var(crudeS), 1/var(crudeU), crudeR, crudeS[crudeblock-12], crudeU, rep(0, nstrain), rep(mean(crudeResults[[1]]), nstrain), initstateD)

  #flag missing data for inference
  y<- ifelse(is.na(y), -1, y)
  yflat<- as.numeric(aperm(y, c(2,1,3)))

  Q_r<- MC_chain[1,3] * RW2PrecMat
  Q_s<- MC_chain[1,4] * RW1PrecMat
  Q_u<- MC_chain[1,5] * R

  #Compute gradients
  Allquantities<- gradmultstrainLoglikelihood2_cpp(y=y, e_it=e_it, nstrain=nstrain,  r=MC_chain[1, 5+(1:time)], s=MC_chain[1, 5+time+(1:12)], u=MC_chain[1, 5+time+12+(1:ndept)], Gamma=G(MC_chain[1,1],MC_chain[1,2]), B=MC_chain[1, 5+time+12+ndept+(1:nstrain)], Bits=Bits, a_k=MC_chain[1, 5+time+12+ndept+nstrain+(1:nstrain)], Model=Model,Q_r=Q_r,Q_s = Q_s,Q_u=Q_u)
  likelihoodcurrent<- Allquantities$loglike
  priorcurrentRcomps<- randomwalk2(MC_chain[1, 5+(1:time)], MC_chain[1, 3])
  priorcurrentScomps<- seasonalComp2(MC_chain[1, 5+time+(1:12)], MC_chain[1, 4], RW1PrecMat)
  priorcurrentUcomps<- logIGMRF1(MC_chain[1, 5+time+12+(1:ndept)], MC_chain[1, 5], R, rankdef)

  grad_current <- list(grad_r=as.numeric(Allquantities$grad_r), grad_s=as.numeric(Allquantities$grad_s), grad_u=as.numeric(Allquantities$grad_u), cov_r=Allquantities$cov_r, cov_s=Allquantities$cov_s)

  for (i in 2:num_iteration) {

    MC_chain[i,3]<- rgamma(1, shape = 1 + (time-2)/2, rate = 0.0001 + (t(MC_chain[i-1, 5+(1:time)]) %*% RW2PrecMat %*% MC_chain[i-1, 5+(1:time)])/2)
    MC_chain[i,4]<- rgamma(1, shape = 1 + 11/2, rate = 0.001 + (t(MC_chain[i-1, 5+time+(1:12)]) %*% RW1PrecMat %*% MC_chain[i-1, 5+time+(1:12)])/2)
    MC_chain[i,5]<- rgamma(1, shape = 1 + (ndept-1)/2, rate = 0.01 + (t(MC_chain[i-1, 5+time+12+(1:ndept)]) %*% R %*% MC_chain[i-1, 5+time+12+(1:ndept)])/2)

    Q_r<- MC_chain[i,3] * RW2PrecMat
    Q_s<- MC_chain[i,4] * RW1PrecMat
    Q_u<- MC_chain[i,5] * R

    current_r <- MC_chain[i-1, 5+(1:time)]
    current_s <- MC_chain[i-1, 5+time+(1:12)]
    current_u <- MC_chain[i-1, 5+time+12+(1:ndept)]

    #Update s
    Mmatcs<- as.numeric(grad_current$cov_s %*% grad_current$grad_s)

    eps_s <- rnorm(12)
    proposedScomps <- as.numeric(MC_chain[i-1, 5+time+(1:12)] + 0.5 * step_sizes$s^2 * Mmatcs + step_sizes$s * chol(grad_current$cov_s) %*% eps_s)
    proposedScomps<- proposedScomps - mean(proposedScomps)

    Allquantities<- gradmultstrainLoglikelihood2_cpp(y=y, e_it=e_it, nstrain=nstrain,  r=current_r, s=proposedScomps, u=current_u, Gamma=G(MC_chain[i-1,1],MC_chain[i-1,2]), B=MC_chain[i-1, 5+time+12+ndept+(1:nstrain)], Bits=Bits, a_k=MC_chain[i-1, 5+time+12+ndept+nstrain+(1:nstrain)], Model=Model,Q_r=Q_r,Q_s = Q_s,Q_u=Q_u)
    grad_proposed <- list(grad_r=as.numeric(Allquantities$grad_r), grad_s=as.numeric(Allquantities$grad_s), grad_u=as.numeric(Allquantities$grad_u), cov_r=Allquantities$cov_r, cov_s=Allquantities$cov_s)

    Mmatps<- as.numeric(grad_proposed$cov_s %*% grad_proposed$grad_s)

    likelihoodproposed<- Allquantities$loglike
    q_prop <- mvnfast::dmvn(proposedScomps, mu = MC_chain[i-1, 5+time+(1:12)] + 0.5 * step_sizes$s^2 * Mmatcs, sigma = grad_current$cov_s * step_sizes$s^2, log = TRUE)
    q_curr <- mvnfast::dmvn(MC_chain[i-1, 5+time+(1:12)], mu = proposedScomps + 0.5 * step_sizes$s^2 * Mmatps, sigma = grad_proposed$cov_s * step_sizes$s^2, log = TRUE)
    priorproposedScomps<- seasonalComp2(proposedScomps, MC_chain[i, 4], RW1PrecMat)

    log_alpha_s <- likelihoodproposed + priorproposedScomps + q_curr - likelihoodcurrent - priorcurrentScomps - q_prop
    if (is.finite(log_alpha_s) && log(runif(1)) < log_alpha_s){
      MC_chain[i, 5+time+(1:12)]<- proposedScomps
      likelihoodcurrent<- likelihoodproposed
      priorcurrentScomps<- priorproposedScomps
      grad_current<- grad_proposed
    }else{
      MC_chain[i, 5+time+(1:12)]<- MC_chain[i-1, 5+time+(1:12)]
    }

    #Update r
    eps_r <- rnorm(time)
    Mmatrc<- as.numeric(grad_current$cov_r %*% grad_current$grad_r)

    proposedRcomps <- as.numeric(current_r + 0.5 * step_sizes$r^2 * Mmatrc + step_sizes$r * chol(grad_current$cov_r) %*% eps_r)
    proposedRcomps<- proposedRcomps - mean(proposedRcomps)

    Allquantities<- gradmultstrainLoglikelihood2_cpp(y=y, e_it=e_it, nstrain=nstrain,  r=proposedRcomps, s=MC_chain[i, 5+time+(1:12)], u=current_u, Gamma=G(MC_chain[i-1,1],MC_chain[i-1,2]), B=MC_chain[i-1, 5+time+12+ndept+(1:nstrain)], Bits=Bits, a_k=MC_chain[i-1, 5+time+12+ndept+nstrain+(1:nstrain)], Model=Model,Q_r=Q_r,Q_s = Q_s,Q_u=Q_u)
    grad_proposed <- list(grad_r=as.numeric(Allquantities$grad_r), grad_s=as.numeric(Allquantities$grad_s), grad_u=as.numeric(Allquantities$grad_u), cov_r=Allquantities$cov_r, cov_s=Allquantities$cov_s)

    Mmatrp<- as.numeric(grad_proposed$cov_r %*% grad_proposed$grad_r)

    q_prop <- mvnfast::dmvn(proposedRcomps, mu = current_r + 0.5 * step_sizes$r^2 * Mmatrc, sigma = grad_current$cov_r * step_sizes$r^2, log = TRUE)
    q_curr <- mvnfast::dmvn(current_r, mu = proposedRcomps + 0.5 * step_sizes$r^2 * Mmatrp, sigma = grad_proposed$cov_r * step_sizes$r^2, log = TRUE)

    likelihoodproposed<- Allquantities$loglike
    priorproposedRcomps <- randomwalk2(proposedRcomps, MC_chain[i, 3])

    log_alpha_r <- likelihoodproposed + priorproposedRcomps + q_curr - likelihoodcurrent - priorcurrentRcomps - q_prop

    if (is.finite(log_alpha_r) && log(runif(1)) < log_alpha_r){
      MC_chain[i, 5+(1:time)] <- proposedRcomps
      likelihoodcurrent<- likelihoodproposed
      priorcurrentRcomps<- priorproposedRcomps
      grad_current<- grad_proposed
    }else{
      MC_chain[i, 5+(1:time)]<- MC_chain[i-1, 5+(1:time)]
    }

    #Update u
    eps_u <- rnorm(ndept)

    proposedUcomps <- as.numeric(MC_chain[i-1, 5+time+12+(1:ndept)] + 0.5 * step_sizes$u^2 * grad_current$grad_u + step_sizes$u * eps_u)
    proposedUcomps<- proposedUcomps - mean(proposedUcomps)

    Allquantities<- gradmultstrainLoglikelihood2_cpp(y=y, e_it=e_it, nstrain=nstrain,  r=MC_chain[i, 5+(1:time)], s=MC_chain[i, 5+time+(1:12)], u=proposedUcomps, Gamma=G(MC_chain[i-1,1],MC_chain[i-1,2]), B=MC_chain[i-1, 5+time+12+ndept+(1:nstrain)], Bits=Bits, a_k=MC_chain[i-1, 5+time+12+ndept+nstrain+(1:nstrain)], Model=Model,Q_r=Q_r,Q_s = Q_s,Q_u=Q_u)
    grad_proposed <- list(grad_r=as.numeric(Allquantities$grad_r), grad_s=as.numeric(Allquantities$grad_s), grad_u=as.numeric(Allquantities$grad_u), cov_r=Allquantities$cov_r, cov_s=Allquantities$cov_s)

    likelihoodproposed<- Allquantities$loglike

    q_prop <- sum(dnorm(proposedUcomps, mean = MC_chain[i-1, 5+time+12+(1:ndept)] + 0.5 * step_sizes$u^2 * grad_current$grad_u, sd = step_sizes$u, log = TRUE))
    q_curr <- sum(dnorm(MC_chain[i-1, 5+time+12+(1:ndept)], mean = proposedUcomps + 0.5 * step_sizes$u^2 * grad_proposed$grad_u, sd = step_sizes$u, log = TRUE))

    priorproposedUcomps<- logIGMRF1(proposedUcomps, MC_chain[i, 5], R, rankdef)
    priorcurrentUcomps<- logIGMRF1(MC_chain[i-1, 5+time+12+(1:ndept)], MC_chain[i, 5], R, rankdef)

    log_alpha_u <- likelihoodproposed + priorproposedUcomps + q_curr - likelihoodcurrent - priorcurrentUcomps - q_prop
    if (is.finite(log_alpha_u) && log(runif(1)) < log_alpha_u){
      MC_chain[i, 5+time+12+(1:ndept)]<- proposedUcomps
      likelihoodcurrent<- likelihoodproposed
      #priorcurrentUcomps<- priorproposedUcomps
      grad_current<- grad_proposed
    }else{
      MC_chain[i, 5+time+12+(1:ndept)]<- MC_chain[i-1, 5+time+12+(1:ndept)]
    }

    if(Model == 0){
      MC_chain[i, 5+time+12+ndept+(1:nstrain)]<-  MC_chain[i-1, 5+time+12+ndept+(1:nstrain)]
      MC_chain[i, 1:2]<- MC_chain[i-1, 1:2]
    }else{
      proposedB <- abs(rnorm(nstrain, mean = MC_chain[i-1, 5+time+12+ndept+(1:nstrain)], sd = rep(sdBs, nstrain)))
      priorcurrentB<- sum(dgamma(MC_chain[i-1, 5+time+12+ndept+(1:nstrain)], shape = rep(2, nstrain), rate = rep(2,nstrain), log=TRUE))
      priorproposedB<- sum(dgamma(proposedB, shape = rep(2, nstrain), rate = rep(2, nstrain), log=TRUE))

      Allquantities<- gradmultstrainLoglikelihood2_cpp(y=y, e_it=e_it, nstrain=nstrain,  r=MC_chain[i, 5+(1:time)], s=MC_chain[i, 5+time+(1:12)], u=MC_chain[i, 5+time+12+(1:ndept)], Gamma=G(MC_chain[i-1,1],MC_chain[i-1,2]), B=proposedB, Bits=Bits, a_k=MC_chain[i-1, 5+time+12+ndept+nstrain+(1:nstrain)], Model=Model,Q_r=Q_r,Q_s = Q_s,Q_u=Q_u)
      grad_proposed <- list(grad_r=as.numeric(Allquantities$grad_r), grad_s=as.numeric(Allquantities$grad_s), grad_u=as.numeric(Allquantities$grad_u), cov_r=Allquantities$cov_r, cov_s=Allquantities$cov_s)

      likelihoodproposed<- Allquantities$loglike

      mh.ratio<- exp(likelihoodproposed + priorproposedB
                     - likelihoodcurrent - priorcurrentB)

      #print(paste("mh.ratioB = ", mh.ratio))

      if(!is.na(mh.ratio) && runif(1) < mh.ratio){
        MC_chain[i, 5+time+12+ndept+(1:nstrain)]<- proposedB
        likelihoodcurrent<- likelihoodproposed
        grad_current<- grad_proposed
      }
      else{
        MC_chain[i, 5+time+12+ndept+(1:nstrain)]<- MC_chain[i-1, 5+time+12+ndept+(1:nstrain)]
      }

      proposedGs<- abs(rnorm(2,mean=c(MC_chain[i-1,1], MC_chain[i-1, 2]), sd=c(sdGs, sdGs)))
      if(proposedGs[1]>1) proposedGs[1]=2-proposedGs[1]
      if(proposedGs[2]>1) proposedGs[2]=2-proposedGs[2]

      priorcurrentGs<- sum(dbeta(MC_chain[i-1,1:2], shape1 = c(2,2), shape2 = c(2,2), log=TRUE))
      priorproposedGs<- sum(dbeta(proposedGs, shape1 = c(2,2), shape2 = c(2,2), log=TRUE))

      Allquantities<- gradmultstrainLoglikelihood2_cpp(y=y, e_it=e_it, nstrain=nstrain,  r=MC_chain[i, 5+(1:time)], s=MC_chain[i, 5+time+(1:12)], u=MC_chain[i, 5+time+12+(1:ndept)], Gamma=G(proposedGs[1],proposedGs[2]), B=MC_chain[i, 5+time+12+ndept+(1:nstrain)], Bits=Bits, a_k=MC_chain[i-1, 5+time+12+ndept+nstrain+(1:nstrain)], Model=Model,Q_r=Q_r,Q_s = Q_s,Q_u=Q_u)
      grad_proposed <- list(grad_r=as.numeric(Allquantities$grad_r), grad_s=as.numeric(Allquantities$grad_s), grad_u=as.numeric(Allquantities$grad_u), cov_r=Allquantities$cov_r, cov_s=Allquantities$cov_s)

      likelihoodproposed<- Allquantities$loglike

      mh.ratio<- exp(likelihoodproposed + priorproposedGs
                     - likelihoodcurrent - priorcurrentGs)

      #print(mh.ratio)

      if(!is.na(mh.ratio) && runif(1) < mh.ratio){
        MC_chain[i, 1:2]<- proposedGs
        likelihoodcurrent<- likelihoodproposed
        grad_current<- grad_proposed
      }
      else{
        MC_chain[i, 1:2]<- MC_chain[i-1,1:2]
      }
    }
    MC_chain[i, 5+time+12+ndept+nstrain+nstrain+1]<- stationarydist(G(MC_chain[i, 1], MC_chain[i, 2]))[2]

    #Gibbs A_k's update
    MC_chain[i, 5+time+12+ndept+nstrain+(1:nstrain)]<- log(rgamma(nstrain, shape = 0.01+SumYk_vec, rate = as.numeric(Allquantities$poisMean4GibbsUpdate) + 0.01/exp(-15)))

    if(i %% 1000 == 0) cat("Iteration:", i, "\n")
  }
  colnames(MC_chain) <- paste(c("G12", "G21", "kappa_r", "kappa_s", "kappa_u", paste("r", 1:time, sep=""), paste("s", 1:12, sep=""), paste("u", 1:ndept, sep=""), paste("B", 1:nstrain, sep=""), paste("a_k", 1:nstrain, sep=""), "StationaryDistribution"))
  MC_chain<- as.data.frame(MC_chain)
  end_time <- Sys.time()
  time_taken<- end_time - start_time
  print(time_taken)
  return(MC_chain)
}


#Full Riemann Manifold Langevin updates in C++
fullCPPmultMMALAInference<- function(y, e_it, Model, adjmat, step_sizes, num_iteration = 1000, independentChains=0) {
  start_time <- Sys.time()
  ndept <- nrow(e_it)
  time <- ncol(e_it)
  nstrain<- dim(y)[3]
  Bits<- encodeBits(nstrain)

  R<- -1 * adjmat
  diag(R)<- -rowSums(R, na.rm = T)
  rankdef<- nrow(R)-qr(R)$rank

  RW1PrecMat<- matrix(0, nrow=12, ncol=12)
  RW1PrecMat[1, ]<- c(2,-1, rep(0, 12-3), -1)
  RW1PrecMat[2, ]<- c(-1,2,-1, rep(0, 12-3))
  RW1PrecMat[3, ]<- c(0, -1,2,-1, rep(0, 12-4))
  RW1PrecMat[(12-1), ]<- c(rep(0, 12-3), -1,2,-1)
  RW1PrecMat[12, ]<- c(-1, rep(0, 12-3), -1, 2)
  for(i in 3:(12-3)){
    RW1PrecMat[i+1, ((i):(i+2))]<- c(-1,2,-1)
  }

  RW2PrecMat<- matrix(0, nrow=time, ncol=time)
  RW2PrecMat[1,(1:3)]<- c(1,-2,1)
  RW2PrecMat[2,(1:4)]<- c(-2,5,-4,1)
  RW2PrecMat[3,(1:5)]<- c(1,-4,6,-4,1)
  RW2PrecMat[(time-1),((time-3):time)]<- c(1,-4,5,-2)
  RW2PrecMat[time,((time-2):time)]<- c(1,-2,1)
  for(i in 3:(time-3)){
    RW2PrecMat[i+1, ((i-1):(i+3))]<- c(1,-4,6,-4,1)
  }

  SumYk_vec<- numeric(nstrain)
  SumYk_vec[1]<- sum(y[,,1])
  sumY<- y[,,1]
  for(k in 2:nstrain){
    sumY<- sumY + y[,,k]
    SumYk_vec[k]<- sum(y[,,k])
  }

  crudeResults<- DetectOutbreaks:::crudeEst(sumY, e_it)
  crudeR<- crudeResults[[1]] - mean(crudeResults[[1]])
  crudeS<- crudeResults[[2]]
  crudeU<- crudeResults[[3]]
  crudeU<- ifelse(is.nan(crudeU), mean(crudeU[is.finite(crudeU)]), crudeU)
  crudeblock<- floor(time/12)
  crudeblock<- ((crudeblock*12)-11):(crudeblock*12)
  meanR<- mean(crudeResults[[1]])

  #flag missing data for inference
  y<- ifelse(is.na(y), -1, y)
  yflat<- as.numeric(aperm(y, c(2,1,3)))

  MC_chain<- MMALA_cpp(y, e_it, Model, Bits, crudeR, crudeS[crudeblock-12], crudeU, RW2PrecMat, RW1PrecMat,
                      R, rankdef, independentChains, num_iteration, meanR, step_sizes)

  colnames(MC_chain) <- paste(c("G12", "G21", "kappa_r", "kappa_s", "kappa_u", paste("r", 1:time, sep=""), paste("s", 1:12, sep=""), paste("u", 1:ndept, sep=""), paste("B", 1:nstrain, sep=""), paste("a_k", 1:nstrain, sep=""), "StationaryDistribution"))
  MC_chain<- as.data.frame(MC_chain)
  end_time <- Sys.time()
  time_taken<- end_time - start_time
  print(time_taken)
  return(MC_chain)
}

#simulatedOutbreakMatrix<- list()
#for(i in 1:nrow(Bits)){
#  simulatedOutbreakMatrix[[i]]<- ifelse(multmod1nstrain5[["states"]]==i-1,1,0)
#}
#Outbreakfigures(matrix_list = simulatedOutbreakMatrix, BitsMatrix = Bits, labelLetter = "A")


#decodedOutbreakMatrix<- list()
#for(i in 1:nrow(Bits)){
#  decodedOutbreakMatrix[[i]]<- multstrain.Decoding(y=multmod1nstrain5[["y"]],e_it=multmod1nstrain5[["e_it"]], nstrain = 5, r=multmod1nstrain5[["r"]], s=multmod1nstrain5[["s"]], u=multmod1nstrain5[["u"]], Gamma=G(0.1,0.2), B=multmod1nstrain5[["B"]], Bits = Bits, a_k = multmod1nstrain5[["a_k"]], state=i)
#}


#decodedOutbreakMatrix<- Posteriormultstrain.Decoding(y=multmod1nstrain5[["y"]], e_it=multmod1nstrain5[["e_it"]], inf.object=MMALAResultscorrectmodel, Modeltype=1, thinningL=1000)
#perstraindecodedOutbreakMatrix<- Posteriormultstrain.Decoding(y=perstrainmultmod1nstrain5[["y"]], e_it=perstrainmultmod1nstrain5[["e_it"]], inf.object=perstrainMMALAResults5strainGibbs[-(1:20000),], Modeltype = 2, thinningL=1000)

#Outbreakfigures(matrix_list = decodedOutbreakMatrix, BitsMatrix = Bits, labelLetter = "B")

#perstrainOutbreakfigures(perstOutP, Outbreaktype = "Simulated outbreaks")
#perstrainOutbreakfigures(perstPostOutP, Outbreaktype = "Decoded outbreaks")
#perstOutP<- perstrainOutbreaks(multmod1nstrain5[["states"]], nstrain = 5)
#perstPostOutP<- perstrainPosteriorOutbreaks(decodedOutbreakMatrix, nstrain = 5)
#perstrainOutbreakfigures(Truth_array = perstOutP, matrix_array = perstPostOutP, Outbreaktype = "Independent_")

#GmultMmalaRes2<- GeneralCPPmultMMALAInference(y=perstrainmultmod1nstrain5[[1]], e_it = perstrainmultmod1nstrain5[[2]], Model = 0, adjmat = sim_adjmat, step_sizes = list("r"=0.3,"s"=0.3,"u"=0.025), num_iteration = 2000)

#Riemann Manifold Langevin updates
GeneralCPPmultMMALAInference<- function(y, e_it, Model, adjmat, step_sizes, num_iteration = 15000, sdGs=0.05, sdBs=0.01, sdAs=0.01){
  start_time <- Sys.time()
  ndept <- nrow(e_it)
  time <- ncol(e_it)
  nstrain<- dim(y)[3]
  Bits<- encodeBits(nstrain)

  R<- -1 * adjmat
  diag(R)<- -rowSums(R, na.rm = T)
  rankdef<- nrow(R)-qr(R)$rank

  RW1PrecMat<- matrix(0, nrow=12, ncol=12)
  RW1PrecMat[1, ]<- c(2,-1, rep(0, 12-3), -1)
  RW1PrecMat[2, ]<- c(-1,2,-1, rep(0, 12-3))
  RW1PrecMat[3, ]<- c(0, -1,2,-1, rep(0, 12-4))
  RW1PrecMat[(12-1), ]<- c(rep(0, 12-3), -1,2,-1)
  RW1PrecMat[12, ]<- c(-1, rep(0, 12-3), -1, 2)
  for(i in 3:(12-3)){
    RW1PrecMat[i+1, ((i):(i+2))]<- c(-1,2,-1)
  }

  RW2PrecMat<- matrix(0, nrow=time, ncol=time)
  RW2PrecMat[1,(1:3)]<- c(1,-2,1)
  RW2PrecMat[2,(1:4)]<- c(-2,5,-4,1)
  RW2PrecMat[3,(1:5)]<- c(1,-4,6,-4,1)
  RW2PrecMat[(time-1),((time-3):time)]<- c(1,-4,5,-2)
  RW2PrecMat[time,((time-2):time)]<- c(1,-2,1)
  for(i in 3:(time-3)){
    RW2PrecMat[i+1, ((i-1):(i+3))]<- c(1,-4,6,-4,1)
  }

  SumYk_vec<- numeric(nstrain)
  SumYk_vec[1]<- sum(y[,,1])
  sumY<- y[,,1]
  for(k in 2:nstrain){
    sumY<- sumY + y[,,k]
    SumYk_vec[k]<- sum(y[,,k])
  }

  crudeResults<- DetectOutbreaks:::crudeEst(sumY, e_it)
  crudeR<- crudeResults[[1]] - mean(crudeResults[[1]])
  crudeS<- crudeResults[[2]]
  crudeU<- crudeResults[[3]]
  crudeU<- ifelse(is.nan(crudeU), mean(crudeU[is.finite(crudeU)]), crudeU)
  crudeblock<- floor(time/12)
  crudeblock<- ((crudeblock*12)-11):(crudeblock*12)

  MC_chain<- matrix(NA, nrow=num_iteration, ncol=2*nstrain+3+time+12+ndept+nstrain+nstrain)
  initGs<- runif(2*nstrain)
  MC_chain[1,]<- c(initGs, 1/var(crudeR), 1/var(crudeS), 1/var(crudeU), crudeR, crudeS[crudeblock-12], crudeU, rep(0, nstrain), rep(mean(crudeResults[[1]]), nstrain))

  #flag missing data for inference
  y<- ifelse(is.na(y), -1, y)
  yflat<- as.numeric(aperm(y, c(2,1,3)))

  Q_r<- MC_chain[1,2*nstrain+1] * RW2PrecMat
  Q_s<- MC_chain[1,2*nstrain+2] * RW1PrecMat
  Q_u<- MC_chain[1,2*nstrain+3] * R

  #Compute gradients
  Gamma_lists<- BuildGamma_list(MC_chain[1,1:(2*nstrain)])
  Allquantities<- perstraingradmultstrainLoglikelihood2_cpp(y=y, e_it=e_it, nstrain=nstrain,  r=MC_chain[1, 2*nstrain+3+(1:time)], s=MC_chain[1, 2*nstrain+3+time+(1:12)], u=MC_chain[1, 2*nstrain+3+time+12+(1:ndept)], Gamma=Gamma_lists, B=MC_chain[1, 2*nstrain+3+time+12+ndept+(1:nstrain)], Bits=Bits, a_k=MC_chain[1, 2*nstrain+3+time+12+ndept+nstrain+(1:nstrain)], Model=Model,Q_r=Q_r,Q_s = Q_s,Q_u=Q_u)
  likelihoodcurrent<- Allquantities$loglike
  priorcurrentRcomps<- randomwalk2(MC_chain[1, 2*nstrain+3+(1:time)], MC_chain[1, 2*nstrain+1])
  priorcurrentScomps<- seasonalComp2(MC_chain[1, 2*nstrain+3+time+(1:12)], MC_chain[1, 2*nstrain+2], RW1PrecMat)
  priorcurrentUcomps<- logIGMRF1(MC_chain[1, 2*nstrain+3+time+12+(1:ndept)], MC_chain[1, 2*nstrain+3], R, rankdef)

  grad_current <- list(grad_r=as.numeric(Allquantities$grad_r), grad_s=as.numeric(Allquantities$grad_s), grad_u=as.numeric(Allquantities$grad_u), cov_r=Allquantities$cov_r, cov_s=Allquantities$cov_s)

  for (i in 2:num_iteration) {

    MC_chain[i,2*nstrain+1]<- rgamma(1, shape = 1 + (time-2)/2, rate = 0.0001 + (t(MC_chain[i-1, 2*nstrain+3+(1:time)]) %*% RW2PrecMat %*% MC_chain[i-1, 2*nstrain+3+(1:time)])/2)
    MC_chain[i,2*nstrain+2]<- rgamma(1, shape = 1 + 11/2, rate = 0.001 + (t(MC_chain[i-1, 2*nstrain+3+time+(1:12)]) %*% RW1PrecMat %*% MC_chain[i-1, 2*nstrain+3+time+(1:12)])/2)
    MC_chain[i,2*nstrain+3]<- rgamma(1, shape = 1 + (ndept-1)/2, rate = 0.01 + (t(MC_chain[i-1, 2*nstrain+3+time+12+(1:ndept)]) %*% R %*% MC_chain[i-1, 2*nstrain+3+time+12+(1:ndept)])/2)

    Q_r<- MC_chain[i,2*nstrain+1] * RW2PrecMat
    Q_s<- MC_chain[i,2*nstrain+2] * RW1PrecMat
    Q_u<- MC_chain[i,2*nstrain+3] * R

    current_r <- MC_chain[i-1, 2*nstrain+3+(1:time)]
    current_s <- MC_chain[i-1, 2*nstrain+3+time+(1:12)]
    current_u <- MC_chain[i-1, 2*nstrain+3+time+12+(1:ndept)]

    #Update s
    Mmatcs<- as.numeric(grad_current$cov_s %*% grad_current$grad_s)

    eps_s <- rnorm(12)
    proposedScomps <- as.numeric(MC_chain[i-1, 2*nstrain+3+time+(1:12)] + 0.5 * step_sizes$s^2 * Mmatcs + step_sizes$s * chol(grad_current$cov_s) %*% eps_s)
    proposedScomps<- proposedScomps - mean(proposedScomps)

    Gamma_lists<- BuildGamma_list(MC_chain[i-1,1:(2*nstrain)])

    Allquantities<- perstraingradmultstrainLoglikelihood2_cpp(y=y, e_it=e_it, nstrain=nstrain,  r=current_r, s=proposedScomps, u=current_u, Gamma=Gamma_lists, B=MC_chain[i-1, 2*nstrain+3+time+12+ndept+(1:nstrain)], Bits=Bits, a_k=MC_chain[i-1, 2*nstrain+3+time+12+ndept+nstrain+(1:nstrain)], Model=Model,Q_r=Q_r,Q_s = Q_s,Q_u=Q_u)
    grad_proposed <- list(grad_r=as.numeric(Allquantities$grad_r), grad_s=as.numeric(Allquantities$grad_s), grad_u=as.numeric(Allquantities$grad_u), cov_r=Allquantities$cov_r, cov_s=Allquantities$cov_s)

    Mmatps<- as.numeric(grad_proposed$cov_s %*% grad_proposed$grad_s)

    likelihoodproposed<- Allquantities$loglike
    q_prop <- mvnfast::dmvn(proposedScomps, mu = MC_chain[i-1, 2*nstrain+3+time+(1:12)] + 0.5 * step_sizes$s^2 * Mmatcs, sigma = grad_current$cov_s * step_sizes$s^2, log = TRUE)
    q_curr <- mvnfast::dmvn(MC_chain[i-1, 2*nstrain+3+time+(1:12)], mu = proposedScomps + 0.5 * step_sizes$s^2 * Mmatps, sigma = grad_proposed$cov_s * step_sizes$s^2, log = TRUE)
    priorproposedScomps<- seasonalComp2(proposedScomps, MC_chain[i, 2*nstrain+2], RW1PrecMat)

    log_alpha_s <- likelihoodproposed + priorproposedScomps + q_curr - likelihoodcurrent - priorcurrentScomps - q_prop
    if (is.finite(log_alpha_s) && log(runif(1)) < log_alpha_s){
      MC_chain[i, 2*nstrain+3+time+(1:12)]<- proposedScomps
      likelihoodcurrent<- likelihoodproposed
      priorcurrentScomps<- priorproposedScomps
      grad_current<- grad_proposed
    }else{
      MC_chain[i, 2*nstrain+3+time+(1:12)]<- MC_chain[i-1, 2*nstrain+3+time+(1:12)]
    }

    #Update r
    eps_r <- rnorm(time)
    Mmatrc<- as.numeric(grad_current$cov_r %*% grad_current$grad_r)

    proposedRcomps <- as.numeric(current_r + 0.5 * step_sizes$r^2 * Mmatrc + step_sizes$r * chol(grad_current$cov_r) %*% eps_r)
    proposedRcomps<- proposedRcomps - mean(proposedRcomps)

    Allquantities<- perstraingradmultstrainLoglikelihood2_cpp(y=y, e_it=e_it, nstrain=nstrain,  r=proposedRcomps, s=MC_chain[i, 2*nstrain+3+time+(1:12)], u=current_u, Gamma=Gamma_lists, B=MC_chain[i-1, 2*nstrain+3+time+12+ndept+(1:nstrain)], Bits=Bits, a_k=MC_chain[i-1, 2*nstrain+3+time+12+ndept+nstrain+(1:nstrain)], Model=Model,Q_r=Q_r,Q_s = Q_s,Q_u=Q_u)
    grad_proposed <- list(grad_r=as.numeric(Allquantities$grad_r), grad_s=as.numeric(Allquantities$grad_s), grad_u=as.numeric(Allquantities$grad_u), cov_r=Allquantities$cov_r, cov_s=Allquantities$cov_s)

    Mmatrp<- as.numeric(grad_proposed$cov_r %*% grad_proposed$grad_r)

    q_prop <- mvnfast::dmvn(proposedRcomps, mu = current_r + 0.5 * step_sizes$r^2 * Mmatrc, sigma = grad_current$cov_r * step_sizes$r^2, log = TRUE)
    q_curr <- mvnfast::dmvn(current_r, mu = proposedRcomps + 0.5 * step_sizes$r^2 * Mmatrp, sigma = grad_proposed$cov_r * step_sizes$r^2, log = TRUE)

    likelihoodproposed<- Allquantities$loglike
    priorproposedRcomps <- randomwalk2(proposedRcomps, MC_chain[i, 2*nstrain+1])

    log_alpha_r <- likelihoodproposed + priorproposedRcomps + q_curr - likelihoodcurrent - priorcurrentRcomps - q_prop

    if (is.finite(log_alpha_r) && log(runif(1)) < log_alpha_r){
      MC_chain[i, 2*nstrain+3+(1:time)] <- proposedRcomps
      likelihoodcurrent<- likelihoodproposed
      priorcurrentRcomps<- priorproposedRcomps
      grad_current<- grad_proposed
    }else{
      MC_chain[i, 2*nstrain+3+(1:time)]<- MC_chain[i-1, 2*nstrain+3+(1:time)]
    }

    #Update u
    eps_u <- rnorm(ndept)

    proposedUcomps <- as.numeric(MC_chain[i-1, 2*nstrain+3+time+12+(1:ndept)] + 0.5 * step_sizes$u^2 * grad_current$grad_u + step_sizes$u * eps_u)
    proposedUcomps<- proposedUcomps - mean(proposedUcomps)

    Allquantities<- perstraingradmultstrainLoglikelihood2_cpp(y=y, e_it=e_it, nstrain=nstrain,  r=MC_chain[i, 2*nstrain+3+(1:time)], s=MC_chain[i, 2*nstrain+3+time+(1:12)], u=proposedUcomps, Gamma=Gamma_lists, B=MC_chain[i-1, 2*nstrain+3+time+12+ndept+(1:nstrain)], Bits=Bits, a_k=MC_chain[i-1, 2*nstrain+3+time+12+ndept+nstrain+(1:nstrain)], Model=Model,Q_r=Q_r,Q_s = Q_s,Q_u=Q_u)
    grad_proposed <- list(grad_r=as.numeric(Allquantities$grad_r), grad_s=as.numeric(Allquantities$grad_s), grad_u=as.numeric(Allquantities$grad_u), cov_r=Allquantities$cov_r, cov_s=Allquantities$cov_s)

    likelihoodproposed<- Allquantities$loglike

    q_prop <- sum(dnorm(proposedUcomps, mean = MC_chain[i-1, 2*nstrain+3+time+12+(1:ndept)] + 0.5 * step_sizes$u^2 * grad_current$grad_u, sd = step_sizes$u, log = TRUE))
    q_curr <- sum(dnorm(MC_chain[i-1, 2*nstrain+3+time+12+(1:ndept)], mean = proposedUcomps + 0.5 * step_sizes$u^2 * grad_proposed$grad_u, sd = step_sizes$u, log = TRUE))

    priorproposedUcomps<- logIGMRF1(proposedUcomps, MC_chain[i, 2*nstrain+3], R, rankdef)
    priorcurrentUcomps<- logIGMRF1(MC_chain[i-1, 2*nstrain+3+time+12+(1:ndept)], MC_chain[i, 2*nstrain+3], R, rankdef)

    log_alpha_u <- likelihoodproposed + priorproposedUcomps + q_curr - likelihoodcurrent - priorcurrentUcomps - q_prop
    if (is.finite(log_alpha_u) && log(runif(1)) < log_alpha_u){
      MC_chain[i, 2*nstrain+3+time+12+(1:ndept)]<- proposedUcomps
      likelihoodcurrent<- likelihoodproposed
      #priorcurrentUcomps<- priorproposedUcomps
      grad_current<- grad_proposed
    }else{
      MC_chain[i, 2*nstrain+3+time+12+(1:ndept)]<- MC_chain[i-1, 2*nstrain+3+time+12+(1:ndept)]
    }

    if(Model == 0){
      MC_chain[i, 2*nstrain+3+time+12+ndept+(1:nstrain)]<-  MC_chain[i-1, 2*nstrain+3+time+12+ndept+(1:nstrain)]
      MC_chain[i, 1:(2*nstrain)]<- MC_chain[i-1, 1:(2*nstrain)]
    }else{
      proposedB <- abs(rnorm(nstrain, mean = MC_chain[i-1, 2*nstrain+3+time+12+ndept+(1:nstrain)], sd = rep(sdBs, nstrain)))
      priorcurrentB<- sum(dgamma(MC_chain[i-1, 2*nstrain+3+time+12+ndept+(1:nstrain)], shape = rep(2, nstrain), rate = rep(2,nstrain), log=TRUE))
      priorproposedB<- sum(dgamma(proposedB, shape = rep(2, nstrain), rate = rep(2, nstrain), log=TRUE))

      Allquantities<- perstraingradmultstrainLoglikelihood2_cpp(y=y, e_it=e_it, nstrain=nstrain,  r=MC_chain[i, 2*nstrain+3+(1:time)], s=MC_chain[i, 2*nstrain+3+time+(1:12)], u=MC_chain[i, 2*nstrain+3+time+12+(1:ndept)], Gamma=Gamma_lists, B=proposedB, Bits=Bits, a_k=MC_chain[i-1, 2*nstrain+3+time+12+ndept+nstrain+(1:nstrain)], Model=Model,Q_r=Q_r,Q_s = Q_s,Q_u=Q_u)
      grad_proposed <- list(grad_r=as.numeric(Allquantities$grad_r), grad_s=as.numeric(Allquantities$grad_s), grad_u=as.numeric(Allquantities$grad_u), cov_r=Allquantities$cov_r, cov_s=Allquantities$cov_s)

      likelihoodproposed<- Allquantities$loglike

      mh.ratio<- exp(likelihoodproposed + priorproposedB
                     - likelihoodcurrent - priorcurrentB)

      #print(paste("mh.ratioB = ", mh.ratio))

      if(!is.na(mh.ratio) && runif(1) < mh.ratio){
        MC_chain[i, 2*nstrain+3+time+12+ndept+(1:nstrain)]<- proposedB
        likelihoodcurrent<- likelihoodproposed
        grad_current<- grad_proposed
      }
      else{
        MC_chain[i, 2*nstrain+3+time+12+ndept+(1:nstrain)]<- MC_chain[i-1, 2*nstrain+3+time+12+ndept+(1:nstrain)]
      }

      proposedGs<- abs(rnorm(2*nstrain,mean=MC_chain[i-1,1:(2*nstrain)], sd=rep(sdGs, 2*nstrain)))
      proposedGs<- ifelse(proposedGs<1, proposedGs, 2-proposedGs)

      priorcurrentGs<- sum(dbeta(MC_chain[i-1,1:(2*nstrain)], shape1 = rep(2,2*nstrain), shape2 = rep(2,2*nstrain), log=TRUE))
      priorproposedGs<- sum(dbeta(proposedGs, shape1 = rep(2,2*nstrain), shape2 = rep(2,2*nstrain), log=TRUE))

      Gamma_lists<- BuildGamma_list(proposedGs)

      Allquantities<- perstraingradmultstrainLoglikelihood2_cpp(y=y, e_it=e_it, nstrain=nstrain,  r=MC_chain[i, 2*nstrain+3+(1:time)], s=MC_chain[i, 2*nstrain+3+time+(1:12)], u=MC_chain[i, 2*nstrain+3+time+12+(1:ndept)], Gamma=Gamma_lists, B=MC_chain[i, 2*nstrain+3+time+12+ndept+(1:nstrain)], Bits=Bits, a_k=MC_chain[i-1, 2*nstrain+3+time+12+ndept+nstrain+(1:nstrain)], Model=Model,Q_r=Q_r,Q_s = Q_s,Q_u=Q_u)
      grad_proposed <- list(grad_r=as.numeric(Allquantities$grad_r), grad_s=as.numeric(Allquantities$grad_s), grad_u=as.numeric(Allquantities$grad_u), cov_r=Allquantities$cov_r, cov_s=Allquantities$cov_s)

      likelihoodproposed<- Allquantities$loglike

      mh.ratio<- exp(likelihoodproposed + priorproposedGs
                     - likelihoodcurrent - priorcurrentGs)

      #print(mh.ratio)

      if(!is.na(mh.ratio) && runif(1) < mh.ratio){
        MC_chain[i, 1:(2*nstrain)]<- proposedGs
        likelihoodcurrent<- likelihoodproposed
        grad_current<- grad_proposed
      }
      else{
        MC_chain[i, 1:(2*nstrain)]<- MC_chain[i-1,1:(2*nstrain)]
      }
    }

    #Gibbs A_k's update
    MC_chain[i, 2*nstrain+3+time+12+ndept+nstrain+(1:nstrain)]<- log(rgamma(nstrain, shape = 0.01+SumYk_vec, rate = as.numeric(Allquantities$poisMean4GibbsUpdate) + 0.01/exp(-15)))

    if(i %% 1000 == 0) cat("Iteration:", i, "\n")
  }
  colnames(MC_chain) <- paste(c(paste0(rep(c("G12", "G21"), nstrain), "Strain", rep(1:nstrain,each=2)), "kappa_r", "kappa_s", "kappa_u", paste("r", 1:time, sep=""), paste("s", 1:12, sep=""), paste("u", 1:ndept, sep=""), paste("B", 1:nstrain, sep=""), paste("a_k", 1:nstrain, sep="")))
  MC_chain<- as.data.frame(MC_chain)
  end_time <- Sys.time()
  time_taken<- end_time - start_time
  print(time_taken)
  return(MC_chain)
}

#sourceCpp("cppFun.cpp")
#replicate(300, SpatMet:::perstraingradmultstrainLoglikelihood2_cpp(
#  y=testdata[["y"]], e_it=testdata[["e_it"]], nstrain=5,
#  r=testdata[["r"]], s=testdata[["s"]], u=testdata[["u"]],
#  Gamma=SpatMet:::BuildGamma_list(testdata[["T.probs"]]), B=testdata[["B"]],
#  Bits=SpatMet:::encodeBits(5), a_k=testdata[["a_k"]],
#  Model=1, Q_r=diag(1, 60), Q_s=diag(1, 12), Q_u=diag(1,9)
#)$loglike)
#replicate(300, SpatMet:::gradmultstrainLoglikelihood2_cpp(
#  y=testdata[["y"]], e_it=testdata[["e_it"]], nstrain=5,
#  r=testdata[["r"]], s=testdata[["s"]], u=testdata[["u"]],
#  Gamma=SpatMet:::G(0.1,0.2), B=testdata[["B"]],
#  Bits=SpatMet:::encodeBits(5), a_k=testdata[["a_k"]],
#  Model=1, Q_r=diag(1, 60), Q_s=diag(1, 12), Q_u=diag(1,9)
#)$loglike)


#set.seed(212);dependentmultmod1nstrain2<- Multstrain.simulate(Model = 1, time=60, adj.matrix = sim_adjmat, Modeltype = 3, nstrain=2, B=c(1.65,1.4));
#dependentfit<- dependentCPPmultMMALAInference(y=dependentmultmod1nstrain2[[1]], e_it = dependentmultmod1nstrain2[[2]], Model = 1, adjmat = sim_adjmat, step_sizes = list("r"=0.3,"s"=0.3,"u"=0.025), num_iteration = 250)

#Riemann Manifold Langevin updates
dependentCPPmultMMALAInference<- function(y, e_it, Model, adjmat, step_sizes, num_iteration = 15000, sdBs=0.01){
  start_time <- Sys.time()
  ndept <- nrow(e_it)
  time <- ncol(e_it)
  nstrain<- dim(y)[3]
  nstate<- 2^nstrain
  Bits<- encodeBits(nstrain)

  R<- -1 * adjmat
  diag(R)<- -rowSums(R, na.rm = T)
  rankdef<- nrow(R)-qr(R)$rank

  RW1PrecMat<- matrix(0, nrow=12, ncol=12)
  RW1PrecMat[1, ]<- c(2,-1, rep(0, 12-3), -1)
  RW1PrecMat[2, ]<- c(-1,2,-1, rep(0, 12-3))
  RW1PrecMat[3, ]<- c(0, -1,2,-1, rep(0, 12-4))
  RW1PrecMat[(12-1), ]<- c(rep(0, 12-3), -1,2,-1)
  RW1PrecMat[12, ]<- c(-1, rep(0, 12-3), -1, 2)
  for(i in 3:(12-3)){
    RW1PrecMat[i+1, ((i):(i+2))]<- c(-1,2,-1)
  }

  RW2PrecMat<- matrix(0, nrow=time, ncol=time)
  RW2PrecMat[1,(1:3)]<- c(1,-2,1)
  RW2PrecMat[2,(1:4)]<- c(-2,5,-4,1)
  RW2PrecMat[3,(1:5)]<- c(1,-4,6,-4,1)
  RW2PrecMat[(time-1),((time-3):time)]<- c(1,-4,5,-2)
  RW2PrecMat[time,((time-2):time)]<- c(1,-2,1)
  for(i in 3:(time-3)){
    RW2PrecMat[i+1, ((i-1):(i+3))]<- c(1,-4,6,-4,1)
  }

  SumYk_vec<- numeric(nstrain)
  SumYk_vec[1]<- sum(y[,,1])
  sumY<- y[,,1]
  for(k in 2:nstrain){
    sumY<- sumY + y[,,k]
    SumYk_vec[k]<- sum(y[,,k])
  }

  crudeResults<- DetectOutbreaks:::crudeEst(sumY, e_it)
  crudeR<- crudeResults[[1]] - mean(crudeResults[[1]])
  crudeS<- crudeResults[[2]]
  crudeU<- crudeResults[[3]]
  crudeU<- ifelse(is.nan(crudeU), mean(crudeU[is.finite(crudeU)]), crudeU)
  crudeblock<- floor(time/12)
  crudeblock<- ((crudeblock*12)-11):(crudeblock*12)

  MC_chain<- matrix(NA, nrow=num_iteration, ncol=nstate*nstate+3+time+12+ndept+nstrain+nstrain)
  initGs<- gtools::rdirichlet(nstate, rep(1, nstate))
  MC_chain[1,]<- c(as.numeric(t(initGs)), 1/var(crudeR), 1/var(crudeS), 1/var(crudeU), crudeR, crudeS[crudeblock-12], crudeU, rep(0, nstrain), rep(mean(crudeResults[[1]]), nstrain))

  #flag missing data for inference
  y<- ifelse(is.na(y), -1, y)
  yflat<- as.numeric(aperm(y, c(2,1,3)))

  Q_r<- MC_chain[1,nstate*nstate+1] * RW2PrecMat
  Q_s<- MC_chain[1,nstate*nstate+2] * RW1PrecMat
  Q_u<- MC_chain[1,nstate*nstate+3] * R

  #Compute gradients
  JointTPM<- initGs
  Allquantities<- dependentgradmultstrainLoglikelihood2_cpp(y=y, e_it=e_it, nstrain=nstrain,  r=MC_chain[1, nstate*nstate+3+(1:time)], s=MC_chain[1, nstate*nstate+3+time+(1:12)], u=MC_chain[1, nstate*nstate+3+time+12+(1:ndept)], jointTPM=JointTPM, B=MC_chain[1, nstate*nstate+3+time+12+ndept+(1:nstrain)], Bits=Bits, a_k=MC_chain[1, nstate*nstate+3+time+12+ndept+nstrain+(1:nstrain)], Model=Model,Q_r=Q_r,Q_s = Q_s,Q_u=Q_u, gradients=1)
  likelihoodcurrent<- Allquantities$loglike
  priorcurrentRcomps<- randomwalk2(MC_chain[1, nstate*nstate+3+(1:time)], MC_chain[1, nstate*nstate+1])
  priorcurrentScomps<- seasonalComp2(MC_chain[1, nstate*nstate+3+time+(1:12)], MC_chain[1, nstate*nstate+2], RW1PrecMat)
  priorcurrentUcomps<- logIGMRF1(MC_chain[1, nstate*nstate+3+time+12+(1:ndept)], MC_chain[1, nstate*nstate+3], R, rankdef)

  grad_current <- list(grad_r=as.numeric(Allquantities$grad_r), grad_s=as.numeric(Allquantities$grad_s), grad_u=as.numeric(Allquantities$grad_u), cov_r=Allquantities$cov_r, cov_s=Allquantities$cov_s)

  deltaP<- 1

  for (i in 2:num_iteration) {

    MC_chain[i,nstate*nstate+1]<- rgamma(1, shape = 1 + (time-2)/2, rate = 0.0001 + (t(MC_chain[i-1, nstate*nstate+3+(1:time)]) %*% RW2PrecMat %*% MC_chain[i-1, nstate*nstate+3+(1:time)])/2)
    MC_chain[i,nstate*nstate+2]<- rgamma(1, shape = 1 + 11/2, rate = 0.001 + (t(MC_chain[i-1, nstate*nstate+3+time+(1:12)]) %*% RW1PrecMat %*% MC_chain[i-1, nstate*nstate+3+time+(1:12)])/2)
    MC_chain[i,nstate*nstate+3]<- rgamma(1, shape = 1 + (ndept-1)/2, rate = 0.01 + (t(MC_chain[i-1, nstate*nstate+3+time+12+(1:ndept)]) %*% R %*% MC_chain[i-1, nstate*nstate+3+time+12+(1:ndept)])/2)

    Q_r<- MC_chain[i,nstate*nstate+1] * RW2PrecMat
    Q_s<- MC_chain[i,nstate*nstate+2] * RW1PrecMat
    Q_u<- MC_chain[i,nstate*nstate+3] * R

    current_r <- MC_chain[i-1, nstate*nstate+3+(1:time)]
    current_s <- MC_chain[i-1, nstate*nstate+3+time+(1:12)]
    current_u <- MC_chain[i-1, nstate*nstate+3+time+12+(1:ndept)]

    #Update s
    Mmatcs<- as.numeric(grad_current$cov_s %*% grad_current$grad_s)

    eps_s <- rnorm(12)
    proposedScomps <- as.numeric(MC_chain[i-1, nstate*nstate+3+time+(1:12)] + 0.5 * step_sizes$s^2 * Mmatcs + step_sizes$s * chol(grad_current$cov_s) %*% eps_s)
    proposedScomps<- proposedScomps - mean(proposedScomps)

    JointTPM<- matrix(MC_chain[i-1, 1:(nstate*nstate)], nrow = nstate, byrow = TRUE)

    Allquantities<- dependentgradmultstrainLoglikelihood2_cpp(y=y, e_it=e_it, nstrain=nstrain,  r=current_r, s=proposedScomps, u=current_u, jointTPM=JointTPM, B=MC_chain[i-1, nstate*nstate+3+time+12+ndept+(1:nstrain)], Bits=Bits, a_k=MC_chain[i-1, nstate*nstate+3+time+12+ndept+nstrain+(1:nstrain)], Model=Model,Q_r=Q_r,Q_s = Q_s,Q_u=Q_u, gradients=1)
    grad_proposed <- list(grad_r=as.numeric(Allquantities$grad_r), grad_s=as.numeric(Allquantities$grad_s), grad_u=as.numeric(Allquantities$grad_u), cov_r=Allquantities$cov_r, cov_s=Allquantities$cov_s)

    Mmatps<- as.numeric(grad_proposed$cov_s %*% grad_proposed$grad_s)

    likelihoodproposed<- Allquantities$loglike
    q_prop <- mvnfast::dmvn(proposedScomps, mu = MC_chain[i-1, nstate*nstate+3+time+(1:12)] + 0.5 * step_sizes$s^2 * Mmatcs, sigma = grad_current$cov_s * step_sizes$s^2, log = TRUE)
    q_curr <- mvnfast::dmvn(MC_chain[i-1, nstate*nstate+3+time+(1:12)], mu = proposedScomps + 0.5 * step_sizes$s^2 * Mmatps, sigma = grad_proposed$cov_s * step_sizes$s^2, log = TRUE)
    priorproposedScomps<- seasonalComp2(proposedScomps, MC_chain[i, nstate*nstate+2], RW1PrecMat)

    log_alpha_s <- likelihoodproposed + priorproposedScomps + q_curr - likelihoodcurrent - priorcurrentScomps - q_prop
    if (is.finite(log_alpha_s) && log(runif(1)) < log_alpha_s){
      MC_chain[i, nstate*nstate+3+time+(1:12)]<- proposedScomps
      likelihoodcurrent<- likelihoodproposed
      priorcurrentScomps<- priorproposedScomps
      grad_current<- grad_proposed
    }else{
      MC_chain[i, nstate*nstate+3+time+(1:12)]<- MC_chain[i-1, nstate*nstate+3+time+(1:12)]
    }

    #Update r
    eps_r <- rnorm(time)
    Mmatrc<- as.numeric(grad_current$cov_r %*% grad_current$grad_r)

    proposedRcomps <- as.numeric(current_r + 0.5 * step_sizes$r^2 * Mmatrc + step_sizes$r * chol(grad_current$cov_r) %*% eps_r)
    proposedRcomps<- proposedRcomps - mean(proposedRcomps)

    Allquantities<- dependentgradmultstrainLoglikelihood2_cpp(y=y, e_it=e_it, nstrain=nstrain,  r=proposedRcomps, s=MC_chain[i, nstate*nstate+3+time+(1:12)], u=current_u, jointTPM=JointTPM, B=MC_chain[i-1, nstate*nstate+3+time+12+ndept+(1:nstrain)], Bits=Bits, a_k=MC_chain[i-1, nstate*nstate+3+time+12+ndept+nstrain+(1:nstrain)], Model=Model,Q_r=Q_r,Q_s = Q_s,Q_u=Q_u, gradients=1)
    grad_proposed <- list(grad_r=as.numeric(Allquantities$grad_r), grad_s=as.numeric(Allquantities$grad_s), grad_u=as.numeric(Allquantities$grad_u), cov_r=Allquantities$cov_r, cov_s=Allquantities$cov_s)

    Mmatrp<- as.numeric(grad_proposed$cov_r %*% grad_proposed$grad_r)

    q_prop <- mvnfast::dmvn(proposedRcomps, mu = current_r + 0.5 * step_sizes$r^2 * Mmatrc, sigma = grad_current$cov_r * step_sizes$r^2, log = TRUE)
    q_curr <- mvnfast::dmvn(current_r, mu = proposedRcomps + 0.5 * step_sizes$r^2 * Mmatrp, sigma = grad_proposed$cov_r * step_sizes$r^2, log = TRUE)

    likelihoodproposed<- Allquantities$loglike
    priorproposedRcomps <- randomwalk2(proposedRcomps, MC_chain[i, nstate*nstate+1])

    log_alpha_r <- likelihoodproposed + priorproposedRcomps + q_curr - likelihoodcurrent - priorcurrentRcomps - q_prop

    if (is.finite(log_alpha_r) && log(runif(1)) < log_alpha_r){
      MC_chain[i, nstate*nstate+3+(1:time)] <- proposedRcomps
      likelihoodcurrent<- likelihoodproposed
      priorcurrentRcomps<- priorproposedRcomps
      grad_current<- grad_proposed
    }else{
      MC_chain[i, nstate*nstate+3+(1:time)]<- MC_chain[i-1, nstate*nstate+3+(1:time)]
    }

    #Update u
    eps_u <- rnorm(ndept)

    proposedUcomps <- as.numeric(MC_chain[i-1, nstate*nstate+3+time+12+(1:ndept)] + 0.5 * step_sizes$u^2 * grad_current$grad_u + step_sizes$u * eps_u)
    proposedUcomps<- proposedUcomps - mean(proposedUcomps)

    Allquantities<- dependentgradmultstrainLoglikelihood2_cpp(y=y, e_it=e_it, nstrain=nstrain,  r=MC_chain[i, nstate*nstate+3+(1:time)], s=MC_chain[i, nstate*nstate+3+time+(1:12)], u=proposedUcomps, jointTPM=JointTPM, B=MC_chain[i-1, nstate*nstate+3+time+12+ndept+(1:nstrain)], Bits=Bits, a_k=MC_chain[i-1, nstate*nstate+3+time+12+ndept+nstrain+(1:nstrain)], Model=Model,Q_r=Q_r,Q_s = Q_s,Q_u=Q_u, gradients=1)
    grad_proposed <- list(grad_r=as.numeric(Allquantities$grad_r), grad_s=as.numeric(Allquantities$grad_s), grad_u=as.numeric(Allquantities$grad_u), cov_r=Allquantities$cov_r, cov_s=Allquantities$cov_s)

    likelihoodproposed<- Allquantities$loglike

    q_prop <- sum(dnorm(proposedUcomps, mean = MC_chain[i-1, nstate*nstate+3+time+12+(1:ndept)] + 0.5 * step_sizes$u^2 * grad_current$grad_u, sd = step_sizes$u, log = TRUE))
    q_curr <- sum(dnorm(MC_chain[i-1, nstate*nstate+3+time+12+(1:ndept)], mean = proposedUcomps + 0.5 * step_sizes$u^2 * grad_proposed$grad_u, sd = step_sizes$u, log = TRUE))

    priorproposedUcomps<- logIGMRF1(proposedUcomps, MC_chain[i, nstate*nstate+3], R, rankdef)
    priorcurrentUcomps<- logIGMRF1(MC_chain[i-1, nstate*nstate+3+time+12+(1:ndept)], MC_chain[i, nstate*nstate+3], R, rankdef)

    log_alpha_u <- likelihoodproposed + priorproposedUcomps + q_curr - likelihoodcurrent - priorcurrentUcomps - q_prop
    if (is.finite(log_alpha_u) && log(runif(1)) < log_alpha_u){
      MC_chain[i, nstate*nstate+3+time+12+(1:ndept)]<- proposedUcomps
      likelihoodcurrent<- likelihoodproposed
      #priorcurrentUcomps<- priorproposedUcomps
      grad_current<- grad_proposed
    }else{
      MC_chain[i, nstate*nstate+3+time+12+(1:ndept)]<- MC_chain[i-1, nstate*nstate+3+time+12+(1:ndept)]
    }

    if(Model == 0){
      MC_chain[i, nstate*nstate+3+time+12+ndept+(1:nstrain)]<-  MC_chain[i-1, nstate*nstate+3+time+12+ndept+(1:nstrain)]
      MC_chain[i, 1:(nstate*nstate)]<- MC_chain[i-1, 1:(nstate*nstate)]
    }else{
      proposedB <- abs(rnorm(nstrain, mean = MC_chain[i-1, nstate*nstate+3+time+12+ndept+(1:nstrain)], sd = rep(sdBs, nstrain)))
      priorcurrentB<- sum(dgamma(MC_chain[i-1, nstate*nstate+3+time+12+ndept+(1:nstrain)], shape = rep(2, nstrain), rate = rep(2,nstrain), log=TRUE))
      priorproposedB<- sum(dgamma(proposedB, shape = rep(2, nstrain), rate = rep(2, nstrain), log=TRUE))

      Allquantities<- dependentgradmultstrainLoglikelihood2_cpp(y=y, e_it=e_it, nstrain=nstrain,  r=MC_chain[i, nstate*nstate+3+(1:time)], s=MC_chain[i, nstate*nstate+3+time+(1:12)], u=MC_chain[i, nstate*nstate+3+time+12+(1:ndept)], jointTPM=JointTPM, B=proposedB, Bits=Bits, a_k=MC_chain[i-1, nstate*nstate+3+time+12+ndept+nstrain+(1:nstrain)], Model=Model,Q_r=Q_r,Q_s = Q_s,Q_u=Q_u, gradients=0)
      #grad_proposed <- list(grad_r=as.numeric(Allquantities$grad_r), grad_s=as.numeric(Allquantities$grad_s), grad_u=as.numeric(Allquantities$grad_u), cov_r=Allquantities$cov_r, cov_s=Allquantities$cov_s)

      likelihoodproposed<- Allquantities$loglike

      mh.ratio<- exp(likelihoodproposed + priorproposedB
                     - likelihoodcurrent - priorcurrentB)

      #print(paste("mh.ratioB = ", mh.ratio))

      if(!is.na(mh.ratio) && runif(1) < mh.ratio){
        MC_chain[i, nstate*nstate+3+time+12+ndept+(1:nstrain)]<- proposedB
        likelihoodcurrent<- likelihoodproposed
        #grad_current<- grad_proposed
      }
      else{
        MC_chain[i, nstate*nstate+3+time+12+ndept+(1:nstrain)]<- MC_chain[i-1, nstate*nstate+3+time+12+ndept+(1:nstrain)]
      }

      #Joint transition probability updates
      for(n in 1:nstate){

      index<- nstate * (n-1) + 1

      JointTPM[n, ] <- gtools::rdirichlet(1, rep(1, nstate) + deltaP * MC_chain[i-1, (index:(n*nstate))])

      proposalproposedGs<-  log(gtools::ddirichlet(JointTPM[n, ], MC_chain[i-1, (index:(n*nstate))]))
      proposalcurrentproposedGs<- log(gtools::ddirichlet(MC_chain[i-1, (index:(n*nstate))], JointTPM[n, ]))

      priorcurrentGs<- log(gtools::ddirichlet(MC_chain[i-1, (index:(n*nstate))], rep(1, nstate)))
      priorproposedGs<- log(gtools::ddirichlet(JointTPM[n, ], rep(1, nstate)))

      if(n == nstate){
        Allquantities<- dependentgradmultstrainLoglikelihood2_cpp(y=y, e_it=e_it, nstrain=nstrain,  r=MC_chain[i, nstate*nstate+3+(1:time)], s=MC_chain[i, nstate*nstate+3+time+(1:12)], u=MC_chain[i, nstate*nstate+3+time+12+(1:ndept)], jointTPM=JointTPM, B=MC_chain[i, nstate*nstate+3+time+12+ndept+(1:nstrain)], Bits=Bits, a_k=MC_chain[i-1, nstate*nstate+3+time+12+ndept+nstrain+(1:nstrain)], Model=Model,Q_r=Q_r,Q_s = Q_s,Q_u=Q_u, gradients=1)
        grad_proposed <- list(grad_r=as.numeric(Allquantities$grad_r), grad_s=as.numeric(Allquantities$grad_s), grad_u=as.numeric(Allquantities$grad_u), cov_r=Allquantities$cov_r, cov_s=Allquantities$cov_s)
      }else{
        Allquantities<- dependentgradmultstrainLoglikelihood2_cpp(y=y, e_it=e_it, nstrain=nstrain,  r=MC_chain[i, nstate*nstate+3+(1:time)], s=MC_chain[i, nstate*nstate+3+time+(1:12)], u=MC_chain[i, nstate*nstate+3+time+12+(1:ndept)], jointTPM=JointTPM, B=MC_chain[i, nstate*nstate+3+time+12+ndept+(1:nstrain)], Bits=Bits, a_k=MC_chain[i-1, nstate*nstate+3+time+12+ndept+nstrain+(1:nstrain)], Model=Model,Q_r=Q_r,Q_s = Q_s,Q_u=Q_u, gradients=0)
      }

      likelihoodproposed<- Allquantities$loglike

      mh.ratio<- exp(likelihoodproposed + priorproposedGs + proposalcurrentproposedGs
                     - likelihoodcurrent - priorcurrentGs - proposalproposedGs)

      #print(mh.ratio)

      if(!is.na(mh.ratio) && runif(1) < mh.ratio){
        MC_chain[i, (index:(n*nstate))]<- as.numeric(JointTPM[n, ])
        likelihoodcurrent<- likelihoodproposed
        grad_current<- grad_proposed
        deltaP<- max(0, deltaP-3)
      }
      else{
        MC_chain[i, (index:(n*nstate))]<- MC_chain[i-1, (index:(n*nstate))]
        deltaP<- deltaP + 1
       }
      }
    }

    #Gibbs A_k's update
    MC_chain[i, nstate*nstate+3+time+12+ndept+nstrain+(1:nstrain)]<- log(rgamma(nstrain, shape = 0.01+SumYk_vec, rate = as.numeric(Allquantities$poisMean4GibbsUpdate) + 0.01/exp(-15)))

    if(i %% 1000 == 0) cat("Iteration:", i, "\n")
  }
  colnames(MC_chain) <- paste(c(paste0(rep("G_", nstate*nstate), rep(1:nstate, each=nstate), ",", 1:nstate), "kappa_r", "kappa_s", "kappa_u", paste("r", 1:time, sep=""), paste("s", 1:12, sep=""), paste("u", 1:ndept, sep=""), paste("B", 1:nstrain, sep=""), paste("a_k", 1:nstrain, sep="")))
  MC_chain<- as.data.frame(MC_chain)
  end_time <- Sys.time()
  time_taken<- end_time - start_time
  print(time_taken)
  return(MC_chain)
}


#Riemann Manifold Langevin updates - with copula modelling
CopulaCPPmultMMALAInference<- function(y, e_it, Model, adjmat, step_sizes, num_iteration = 15000, sdGs=0.05, sdBs=0.01, sdAs=0.01, n_copparams=1){
  start_time <- Sys.time()
  ndept <- nrow(e_it)
  time <- ncol(e_it)
  nstrain<- dim(y)[3]
  Bits<- encodeBits(nstrain)

  R<- -1 * adjmat
  diag(R)<- -rowSums(R, na.rm = T)
  rankdef<- nrow(R)-qr(R)$rank

  RW1PrecMat<- matrix(0, nrow=12, ncol=12)
  RW1PrecMat[1, ]<- c(2,-1, rep(0, 12-3), -1)
  RW1PrecMat[2, ]<- c(-1,2,-1, rep(0, 12-3))
  RW1PrecMat[3, ]<- c(0, -1,2,-1, rep(0, 12-4))
  RW1PrecMat[(12-1), ]<- c(rep(0, 12-3), -1,2,-1)
  RW1PrecMat[12, ]<- c(-1, rep(0, 12-3), -1, 2)
  for(i in 3:(12-3)){
    RW1PrecMat[i+1, ((i):(i+2))]<- c(-1,2,-1)
  }

  RW2PrecMat<- matrix(0, nrow=time, ncol=time)
  RW2PrecMat[1,(1:3)]<- c(1,-2,1)
  RW2PrecMat[2,(1:4)]<- c(-2,5,-4,1)
  RW2PrecMat[3,(1:5)]<- c(1,-4,6,-4,1)
  RW2PrecMat[(time-1),((time-3):time)]<- c(1,-4,5,-2)
  RW2PrecMat[time,((time-2):time)]<- c(1,-2,1)
  for(i in 3:(time-3)){
    RW2PrecMat[i+1, ((i-1):(i+3))]<- c(1,-4,6,-4,1)
  }

  SumYk_vec<- numeric(nstrain)
  SumYk_vec[1]<- sum(y[,,1])
  sumY<- y[,,1]
  for(k in 2:nstrain){
    sumY<- sumY + y[,,k]
    SumYk_vec[k]<- sum(y[,,k])
  }

  crudeResults<- DetectOutbreaks:::crudeEst(sumY, e_it)
  crudeR<- crudeResults[[1]] - mean(crudeResults[[1]])
  crudeS<- crudeResults[[2]]
  crudeU<- crudeResults[[3]]
  crudeU<- ifelse(is.nan(crudeU), mean(crudeU[is.finite(crudeU)]), crudeU)
  crudeblock<- floor(time/12)
  crudeblock<- ((crudeblock*12)-11):(crudeblock*12)

  MC_chain<- matrix(NA, nrow=num_iteration, ncol=5+time+12+ndept+nstrain+nstrain+n_copparams)
  initG12<- runif(1)
  initG21<- runif(1)
  initstateD<- stationarydist(G(initG12, initG21))[2]
  MC_chain[1,]<- c(initG12, initG21, 1/var(crudeR), 1/var(crudeS), 1/var(crudeU), crudeR, crudeS[crudeblock-12], crudeU, rep(0, nstrain), rep(mean(crudeResults[[1]]), nstrain), 0.2)

  #flag missing data for inference
  y<- ifelse(is.na(y), -1, y)
  yflat<- as.numeric(aperm(y, c(2,1,3)))

  Q_r<- MC_chain[1,3] * RW2PrecMat
  Q_s<- MC_chain[1,4] * RW1PrecMat
  Q_u<- MC_chain[1,5] * R

  #Compute gradients
  copulaTPM<- JointTransitionMatrix_copula(G(MC_chain[1,1], MC_chain[1,2]), K=nstrain, MC_chain[1, ncol(MC_chain)])
  Allquantities<- copulagradmultstrainLoglikelihood2_cpp(y=y, e_it=e_it, nstrain=nstrain,  r=MC_chain[1, 5+(1:time)], s=MC_chain[1, 5+time+(1:12)], u=MC_chain[1, 5+time+12+(1:ndept)], Gamma=copulaTPM, B=MC_chain[1, 5+time+12+ndept+(1:nstrain)], Bits=Bits, a_k=MC_chain[1, 5+time+12+ndept+nstrain+(1:nstrain)], Model=Model,Q_r=Q_r,Q_s = Q_s,Q_u=Q_u)
  likelihoodcurrent<- Allquantities$loglike
  priorcurrentRcomps<- randomwalk2(MC_chain[1, 5+(1:time)], MC_chain[1, 3])
  priorcurrentScomps<- seasonalComp2(MC_chain[1, 5+time+(1:12)], MC_chain[1, 4], RW1PrecMat)
  priorcurrentUcomps<- logIGMRF1(MC_chain[1, 5+time+12+(1:ndept)], MC_chain[1, 5], R, rankdef)

  grad_current <- list(grad_r=as.numeric(Allquantities$grad_r), grad_s=as.numeric(Allquantities$grad_s), grad_u=as.numeric(Allquantities$grad_u), cov_r=Allquantities$cov_r, cov_s=Allquantities$cov_s)

  #pseudoGibbsAcceptance<- 0

  for (i in 2:num_iteration) {

    MC_chain[i,3]<- rgamma(1, shape = 1 + (time-2)/2, rate = 0.0001 + (t(MC_chain[i-1, 5+(1:time)]) %*% RW2PrecMat %*% MC_chain[i-1, 5+(1:time)])/2)
    MC_chain[i,4]<- rgamma(1, shape = 1 + 11/2, rate = 0.001 + (t(MC_chain[i-1, 5+time+(1:12)]) %*% RW1PrecMat %*% MC_chain[i-1, 5+time+(1:12)])/2)
    MC_chain[i,5]<- rgamma(1, shape = 1 + (ndept-1)/2, rate = 0.01 + (t(MC_chain[i-1, 5+time+12+(1:ndept)]) %*% R %*% MC_chain[i-1, 5+time+12+(1:ndept)])/2)

    Q_r<- MC_chain[i,3] * RW2PrecMat
    Q_s<- MC_chain[i,4] * RW1PrecMat
    Q_u<- MC_chain[i,5] * R

    current_r <- MC_chain[i-1, 5+(1:time)]
    current_s <- MC_chain[i-1, 5+time+(1:12)]
    current_u <- MC_chain[i-1, 5+time+12+(1:ndept)]

    #Update s
    Mmatcs<- as.numeric(grad_current$cov_s %*% grad_current$grad_s)

    eps_s <- rnorm(12)
    proposedScomps <- as.numeric(MC_chain[i-1, 5+time+(1:12)] + 0.5 * step_sizes$s^2 * Mmatcs + step_sizes$s * chol(grad_current$cov_s) %*% eps_s)
    proposedScomps<- proposedScomps - mean(proposedScomps)

    Allquantities<- copulagradmultstrainLoglikelihood2_cpp(y=y, e_it=e_it, nstrain=nstrain,  r=current_r, s=proposedScomps, u=current_u, Gamma=copulaTPM, B=MC_chain[i-1, 5+time+12+ndept+(1:nstrain)], Bits=Bits, a_k=MC_chain[i-1, 5+time+12+ndept+nstrain+(1:nstrain)], Model=Model,Q_r=Q_r,Q_s = Q_s,Q_u=Q_u)
    grad_proposed <- list(grad_r=as.numeric(Allquantities$grad_r), grad_s=as.numeric(Allquantities$grad_s), grad_u=as.numeric(Allquantities$grad_u), cov_r=Allquantities$cov_r, cov_s=Allquantities$cov_s)

    Mmatps<- as.numeric(grad_proposed$cov_s %*% grad_proposed$grad_s)

    likelihoodproposed<- Allquantities$loglike
    q_prop <- mvnfast::dmvn(proposedScomps, mu = MC_chain[i-1, 5+time+(1:12)] + 0.5 * step_sizes$s^2 * Mmatcs, sigma = grad_current$cov_s * step_sizes$s^2, log = TRUE)
    q_curr <- mvnfast::dmvn(MC_chain[i-1, 5+time+(1:12)], mu = proposedScomps + 0.5 * step_sizes$s^2 * Mmatps, sigma = grad_proposed$cov_s * step_sizes$s^2, log = TRUE)
    priorproposedScomps<- seasonalComp2(proposedScomps, MC_chain[i, 4], RW1PrecMat)

    log_alpha_s <- likelihoodproposed + priorproposedScomps + q_curr - likelihoodcurrent - priorcurrentScomps - q_prop
    if (is.finite(log_alpha_s) && log(runif(1)) < log_alpha_s){
      MC_chain[i, 5+time+(1:12)]<- proposedScomps
      likelihoodcurrent<- likelihoodproposed
      priorcurrentScomps<- priorproposedScomps
      grad_current<- grad_proposed
    }else{
      MC_chain[i, 5+time+(1:12)]<- MC_chain[i-1, 5+time+(1:12)]
    }

    #Update r
    eps_r <- rnorm(time)
    Mmatrc<- as.numeric(grad_current$cov_r %*% grad_current$grad_r)

    proposedRcomps <- as.numeric(current_r + 0.5 * step_sizes$r^2 * Mmatrc + step_sizes$r * chol(grad_current$cov_r) %*% eps_r)
    proposedRcomps<- proposedRcomps - mean(proposedRcomps)

    Allquantities<- copulagradmultstrainLoglikelihood2_cpp(y=y, e_it=e_it, nstrain=nstrain,  r=proposedRcomps, s=MC_chain[i, 5+time+(1:12)], u=current_u, Gamma=copulaTPM, B=MC_chain[i-1, 5+time+12+ndept+(1:nstrain)], Bits=Bits, a_k=MC_chain[i-1, 5+time+12+ndept+nstrain+(1:nstrain)], Model=Model,Q_r=Q_r,Q_s = Q_s,Q_u=Q_u)
    grad_proposed <- list(grad_r=as.numeric(Allquantities$grad_r), grad_s=as.numeric(Allquantities$grad_s), grad_u=as.numeric(Allquantities$grad_u), cov_r=Allquantities$cov_r, cov_s=Allquantities$cov_s)

    Mmatrp<- as.numeric(grad_proposed$cov_r %*% grad_proposed$grad_r)

    q_prop <- mvnfast::dmvn(proposedRcomps, mu = current_r + 0.5 * step_sizes$r^2 * Mmatrc, sigma = grad_current$cov_r * step_sizes$r^2, log = TRUE)
    q_curr <- mvnfast::dmvn(current_r, mu = proposedRcomps + 0.5 * step_sizes$r^2 * Mmatrp, sigma = grad_proposed$cov_r * step_sizes$r^2, log = TRUE)

    likelihoodproposed<- Allquantities$loglike
    priorproposedRcomps <- randomwalk2(proposedRcomps, MC_chain[i, 3])

    log_alpha_r <- likelihoodproposed + priorproposedRcomps + q_curr - likelihoodcurrent - priorcurrentRcomps - q_prop

    if (is.finite(log_alpha_r) && log(runif(1)) < log_alpha_r){
      MC_chain[i, 5+(1:time)] <- proposedRcomps
      likelihoodcurrent<- likelihoodproposed
      priorcurrentRcomps<- priorproposedRcomps
      grad_current<- grad_proposed
    }else{
      MC_chain[i, 5+(1:time)]<- MC_chain[i-1, 5+(1:time)]
    }

    #Update u
    eps_u <- rnorm(ndept)

    proposedUcomps <- as.numeric(MC_chain[i-1, 5+time+12+(1:ndept)] + 0.5 * step_sizes$u^2 * grad_current$grad_u + step_sizes$u * eps_u)
    proposedUcomps<- proposedUcomps - mean(proposedUcomps)

    Allquantities<- copulagradmultstrainLoglikelihood2_cpp(y=y, e_it=e_it, nstrain=nstrain,  r=MC_chain[i, 5+(1:time)], s=MC_chain[i, 5+time+(1:12)], u=proposedUcomps, Gamma=copulaTPM, B=MC_chain[i-1, 5+time+12+ndept+(1:nstrain)], Bits=Bits, a_k=MC_chain[i-1, 5+time+12+ndept+nstrain+(1:nstrain)], Model=Model,Q_r=Q_r,Q_s = Q_s,Q_u=Q_u)
    grad_proposed <- list(grad_r=as.numeric(Allquantities$grad_r), grad_s=as.numeric(Allquantities$grad_s), grad_u=as.numeric(Allquantities$grad_u), cov_r=Allquantities$cov_r, cov_s=Allquantities$cov_s)

    likelihoodproposed<- Allquantities$loglike

    q_prop <- sum(dnorm(proposedUcomps, mean = MC_chain[i-1, 5+time+12+(1:ndept)] + 0.5 * step_sizes$u^2 * grad_current$grad_u, sd = step_sizes$u, log = TRUE))
    q_curr <- sum(dnorm(MC_chain[i-1, 5+time+12+(1:ndept)], mean = proposedUcomps + 0.5 * step_sizes$u^2 * grad_proposed$grad_u, sd = step_sizes$u, log = TRUE))

    priorproposedUcomps<- logIGMRF1(proposedUcomps, MC_chain[i, 5], R, rankdef)
    priorcurrentUcomps<- logIGMRF1(MC_chain[i-1, 5+time+12+(1:ndept)], MC_chain[i, 5], R, rankdef)

    log_alpha_u <- likelihoodproposed + priorproposedUcomps + q_curr - likelihoodcurrent - priorcurrentUcomps - q_prop
    if (is.finite(log_alpha_u) && log(runif(1)) < log_alpha_u){
      MC_chain[i, 5+time+12+(1:ndept)]<- proposedUcomps
      likelihoodcurrent<- likelihoodproposed
      #priorcurrentUcomps<- priorproposedUcomps
      grad_current<- grad_proposed
    }else{
      MC_chain[i, 5+time+12+(1:ndept)]<- MC_chain[i-1, 5+time+12+(1:ndept)]
    }

    if(Model == 0){
      MC_chain[i, 5+time+12+ndept+(1:nstrain)]<-  MC_chain[i-1, 5+time+12+ndept+(1:nstrain)]
      MC_chain[i, 1:2]<- MC_chain[i-1, 1:2]
    }else{
      proposedB <- abs(rnorm(nstrain, mean = MC_chain[i-1, 5+time+12+ndept+(1:nstrain)], sd = rep(sdBs, nstrain)))
      priorcurrentB<- sum(dgamma(MC_chain[i-1, 5+time+12+ndept+(1:nstrain)], shape = rep(2, nstrain), rate = rep(2,nstrain), log=TRUE))
      priorproposedB<- sum(dgamma(proposedB, shape = rep(2, nstrain), rate = rep(2, nstrain), log=TRUE))

      Allquantities<- copulagradmultstrainLoglikelihood2_cpp(y=y, e_it=e_it, nstrain=nstrain,  r=MC_chain[i, 5+(1:time)], s=MC_chain[i, 5+time+(1:12)], u=MC_chain[i, 5+time+12+(1:ndept)], Gamma=copulaTPM, B=proposedB, Bits=Bits, a_k=MC_chain[i-1, 5+time+12+ndept+nstrain+(1:nstrain)], Model=Model,Q_r=Q_r,Q_s = Q_s,Q_u=Q_u)
      grad_proposed <- list(grad_r=as.numeric(Allquantities$grad_r), grad_s=as.numeric(Allquantities$grad_s), grad_u=as.numeric(Allquantities$grad_u), cov_r=Allquantities$cov_r, cov_s=Allquantities$cov_s)

      likelihoodproposed<- Allquantities$loglike

      mh.ratio<- exp(likelihoodproposed + priorproposedB
                     - likelihoodcurrent - priorcurrentB)

      #print(paste("mh.ratioB = ", mh.ratio))

      if(!is.na(mh.ratio) && runif(1) < mh.ratio){
        MC_chain[i, 5+time+12+ndept+(1:nstrain)]<- proposedB
        likelihoodcurrent<- likelihoodproposed
        grad_current<- grad_proposed
      }
      else{
        MC_chain[i, 5+time+12+ndept+(1:nstrain)]<- MC_chain[i-1, 5+time+12+ndept+(1:nstrain)]
      }

      proposedGs<- abs(rnorm(2,mean=c(MC_chain[i-1,1], MC_chain[i-1, 2]), sd=c(sdGs, sdGs)))
      if(proposedGs[1]>1) proposedGs[1]=2-proposedGs[1]
      if(proposedGs[2]>1) proposedGs[2]=2-proposedGs[2]

      priorcurrentGs<- sum(dbeta(MC_chain[i-1,1:2], shape1 = c(2,2), shape2 = c(2,2), log=TRUE))
      priorproposedGs<- sum(dbeta(proposedGs, shape1 = c(2,2), shape2 = c(2,2), log=TRUE))

      copulaTPM<- JointTransitionMatrix_copula(G(proposedGs[1], proposedGs[2]), K=nstrain, MC_chain[i-1, ncol(MC_chain)])
      Allquantities<- copulagradmultstrainLoglikelihood2_cpp(y=y, e_it=e_it, nstrain=nstrain,  r=MC_chain[i, 5+(1:time)], s=MC_chain[i, 5+time+(1:12)], u=MC_chain[i, 5+time+12+(1:ndept)], Gamma=copulaTPM, B=MC_chain[i, 5+time+12+ndept+(1:nstrain)], Bits=Bits, a_k=MC_chain[i-1, 5+time+12+ndept+nstrain+(1:nstrain)], Model=Model,Q_r=Q_r,Q_s = Q_s,Q_u=Q_u)
      grad_proposed <- list(grad_r=as.numeric(Allquantities$grad_r), grad_s=as.numeric(Allquantities$grad_s), grad_u=as.numeric(Allquantities$grad_u), cov_r=Allquantities$cov_r, cov_s=Allquantities$cov_s)

      likelihoodproposed<- Allquantities$loglike

      mh.ratio<- exp(likelihoodproposed + priorproposedGs
                     - likelihoodcurrent - priorcurrentGs)

      #print(mh.ratio)

      if(!is.na(mh.ratio) && runif(1) < mh.ratio){
        MC_chain[i, 1:2]<- proposedGs
        likelihoodcurrent<- likelihoodproposed
        grad_current<- grad_proposed
      }
      else{
        MC_chain[i, 1:2]<- MC_chain[i-1,1:2]
      }
      proposedcopPs<- rnorm(n_copparams,mean=log(MC_chain[i-1, ncol(MC_chain)]), sd=rep(0.7, n_copparams))

      #proposedcopPs<- rnorm(1, mean=MC_chain[i-1, ncol(MC_chain)], sd=0.4)
      #if(proposedcopPs< -1 || proposedcopPs> 1) proposedcopPs = MC_chain[i-1, ncol(MC_chain)]

      #priorcurrentcopPs<- dgamma(MC_chain[i-1, ncol(MC_chain)], shape = 1, rate = 0.001, log=TRUE)
      #priorproposedcopPs<- dgamma(proposedcopPs, shape = 1, rate = 0.001, log=TRUE)

      proposalcurrentcop<- dnorm(log(MC_chain[i-1, ncol(MC_chain)]), mean=proposedcopPs, sd=rep(0.7, n_copparams), log = T) + MC_chain[i-1, ncol(MC_chain)]
      proposalproposedcop<- dnorm(proposedcopPs, mean=log(MC_chain[i-1, ncol(MC_chain)]), sd=rep(0.7, n_copparams), log = T) + exp(proposedcopPs)

      copulaTPM<- JointTransitionMatrix_copula(G(MC_chain[i, 1],MC_chain[i, 2]), K=nstrain, exp(proposedcopPs))
      #copulaTPM<- JointTransitionMatrix_copula(G(MC_chain[i, 1],MC_chain[i, 2]), K=nstrain, proposedcopPs)

      Allquantities<- copulagradmultstrainLoglikelihood2_cpp(y=y, e_it=e_it, nstrain=nstrain,  r=MC_chain[i, 5+(1:time)], s=MC_chain[i, 5+time+(1:12)], u=MC_chain[i, 5+time+12+(1:ndept)], Gamma=copulaTPM, B=MC_chain[i, 5+time+12+ndept+(1:nstrain)], Bits=Bits, a_k=MC_chain[i-1, 5+time+12+ndept+nstrain+(1:nstrain)], Model=Model,Q_r=Q_r,Q_s = Q_s,Q_u=Q_u)
      grad_proposed <- list(grad_r=as.numeric(Allquantities$grad_r), grad_s=as.numeric(Allquantities$grad_s), grad_u=as.numeric(Allquantities$grad_u), cov_r=Allquantities$cov_r, cov_s=Allquantities$cov_s)

      likelihoodproposed<- Allquantities$loglike

      mh.ratio<- exp(likelihoodproposed - likelihoodcurrent + proposalcurrentcop - proposalproposedcop)
                     #+ priorproposedcopPs - priorcurrentcopPs)

      #print(mh.ratio)

      if(!is.na(mh.ratio) && runif(1) < mh.ratio){
        MC_chain[i, ncol(MC_chain)]<- exp(proposedcopPs)
        #MC_chain[i, ncol(MC_chain)]<- proposedcopPs
        likelihoodcurrent<- likelihoodproposed
        grad_current<- grad_proposed
      }
      else{
        MC_chain[i, ncol(MC_chain)]<- MC_chain[i-1, ncol(MC_chain)]
      }
    }

    #Pseudo-Gibbs A_k's update
    propShape<- 0.01+SumYk_vec
    propRate<- as.numeric(Allquantities$poisMean4GibbsUpdate) + 0.01/exp(-15)
    proposedAks<- log(rgamma(nstrain, shape = propShape, rate = propRate))
    copulaTPM<- JointTransitionMatrix_copula(G(MC_chain[i, 1],MC_chain[i, 2]), K=nstrain, MC_chain[i, ncol(MC_chain)])
    Allquantities<- copulagradmultstrainLoglikelihood2_cpp(y=y, e_it=e_it, nstrain=nstrain,  r=MC_chain[i, 5+(1:time)], s=MC_chain[i, 5+time+(1:12)], u=MC_chain[i, 5+time+12+(1:ndept)], Gamma=copulaTPM, B=MC_chain[i, 5+time+12+ndept+(1:nstrain)], Bits=Bits, a_k=proposedAks, Model=Model,Q_r=Q_r,Q_s = Q_s,Q_u=Q_u)
    grad_proposed <- list(grad_r=as.numeric(Allquantities$grad_r), grad_s=as.numeric(Allquantities$grad_s), grad_u=as.numeric(Allquantities$grad_u), cov_r=Allquantities$cov_r, cov_s=Allquantities$cov_s)
    likelihoodproposed<- Allquantities$loglike

    proposalcurrentAks <- sum(dgamma(exp(MC_chain[i-1, 5+time+12+ndept+nstrain+(1:nstrain)]),
                                 shape=propShape, rate=propRate, log=TRUE) + MC_chain[i-1, 5+time+12+ndept+nstrain+(1:nstrain)])
    proposalproposedAks <- sum(dgamma(exp(proposedAks),
                                  shape=propShape, rate=propRate, log=TRUE) + proposedAks)

    priorcurrentAks <- sum(dgamma(exp(MC_chain[i-1, 5+time+12+ndept+nstrain+(1:nstrain)]),
                              shape=0.01, rate=0.01/exp(-15), log=TRUE) + MC_chain[i-1, 5+time+12+ndept+nstrain+(1:nstrain)])
    priorproposedAks <- sum(dgamma(exp(proposedAks),
                               shape=0.01, rate=0.01/exp(-15), log=TRUE) +  proposedAks)

    mh.ratio<- exp(likelihoodproposed - likelihoodcurrent + proposalcurrentAks - proposalproposedAks
                   + priorproposedAks - priorcurrentAks)

    #print(mh.ratio)

    if(!is.na(mh.ratio) && runif(1) < mh.ratio){
      MC_chain[i, 5+time+12+ndept+nstrain+(1:nstrain)]<- proposedAks
      likelihoodcurrent<- likelihoodproposed
      grad_current<- grad_proposed
      #pseudoGibbsAcceptance<- pseudoGibbsAcceptance + 1
    }
    else{
      MC_chain[i, 5+time+12+ndept+nstrain+(1:nstrain)]<- MC_chain[i-1, 5+time+12+ndept+nstrain+(1:nstrain)]
    }

    if(i %% 1000 == 0) cat("Iteration:", i, "\n")
  }
  colnames(MC_chain) <- paste(c("G12", "G21", "kappa_r", "kappa_s", "kappa_u", paste("r", 1:time, sep=""), paste("s", 1:12, sep=""), paste("u", 1:ndept, sep=""), paste("B", 1:nstrain, sep=""), paste("a_k", 1:nstrain, sep=""), paste("copulaParam", 1:n_copparams, sep="")))
  MC_chain<- as.data.frame(MC_chain)
  end_time <- Sys.time()
  time_taken<- end_time - start_time
  print(time_taken)
  #print(pseudoGibbsAcceptance/num_iteration)
  return(MC_chain)
}



#Riemann Manifold Langevin updates - with copula modelling
CopulaCPPmultMMALAInference_per_strain<- function(y, e_it, Model, adjmat, step_sizes, num_iteration = 15000, sdGs=0.05, sdBs=0.01, sdAs=0.01, n_copparams=1){
  start_time <- Sys.time()
  ndept <- nrow(e_it)
  time <- ncol(e_it)
  nstrain<- dim(y)[3]
  Bits<- encodeBits(nstrain)

  R<- -1 * adjmat
  diag(R)<- -rowSums(R, na.rm = T)
  rankdef<- nrow(R)-qr(R)$rank

  RW1PrecMat<- matrix(0, nrow=12, ncol=12)
  RW1PrecMat[1, ]<- c(2,-1, rep(0, 12-3), -1)
  RW1PrecMat[2, ]<- c(-1,2,-1, rep(0, 12-3))
  RW1PrecMat[3, ]<- c(0, -1,2,-1, rep(0, 12-4))
  RW1PrecMat[(12-1), ]<- c(rep(0, 12-3), -1,2,-1)
  RW1PrecMat[12, ]<- c(-1, rep(0, 12-3), -1, 2)
  for(i in 3:(12-3)){
    RW1PrecMat[i+1, ((i):(i+2))]<- c(-1,2,-1)
  }

  RW2PrecMat<- matrix(0, nrow=time, ncol=time)
  RW2PrecMat[1,(1:3)]<- c(1,-2,1)
  RW2PrecMat[2,(1:4)]<- c(-2,5,-4,1)
  RW2PrecMat[3,(1:5)]<- c(1,-4,6,-4,1)
  RW2PrecMat[(time-1),((time-3):time)]<- c(1,-4,5,-2)
  RW2PrecMat[time,((time-2):time)]<- c(1,-2,1)
  for(i in 3:(time-3)){
    RW2PrecMat[i+1, ((i-1):(i+3))]<- c(1,-4,6,-4,1)
  }

  SumYk_vec<- numeric(nstrain)
  SumYk_vec[1]<- sum(y[,,1])
  sumY<- y[,,1]
  for(k in 2:nstrain){
    sumY<- sumY + y[,,k]
    SumYk_vec[k]<- sum(y[,,k])
  }

  crudeResults<- DetectOutbreaks:::crudeEst(sumY, e_it)
  crudeR<- crudeResults[[1]] - mean(crudeResults[[1]])
  crudeS<- crudeResults[[2]]
  crudeU<- crudeResults[[3]]
  crudeU<- ifelse(is.nan(crudeU), mean(crudeU[is.finite(crudeU)]), crudeU)
  crudeblock<- floor(time/12)
  crudeblock<- ((crudeblock*12)-11):(crudeblock*12)

  MC_chain<- matrix(NA, nrow=num_iteration, ncol=2*nstrain+3+time+12+ndept+nstrain+nstrain+n_copparams)
  initGs<- runif(2*nstrain)
  MC_chain[1, ]<- c(initGs, 1/var(crudeR), 1/var(crudeS), 1/var(crudeU), crudeR, crudeS[crudeblock-12], crudeU, rep(0, nstrain), rep(mean(crudeResults[[1]]), nstrain), rep(0.5, n_copparams))

  #flag missing data for inference
  y<- ifelse(is.na(y), -1, y)
  yflat<- as.numeric(aperm(y, c(2,1,3)))

  Q_r<- MC_chain[1,2*nstrain+1] * RW2PrecMat
  Q_s<- MC_chain[1,2*nstrain+2] * RW1PrecMat
  Q_u<- MC_chain[1,2*nstrain+3] * R

  #Compute gradients
  Gamma_lists<- BuildGamma_list(MC_chain[1, 1:(2*nstrain)])
  copulaTPM<- JointTransitionMatrix_copula_per_strain(Gamma_lists, MC_chain[1, ncol(MC_chain)])
  Allquantities<- copulagradmultstrainLoglikelihood2_cpp(y=y, e_it=e_it, nstrain=nstrain,  r=MC_chain[1, 2*nstrain+3+(1:time)], s=MC_chain[1, 2*nstrain+3+time+(1:12)], u=MC_chain[1, 2*nstrain+3+time+12+(1:ndept)], Gamma=copulaTPM, B=MC_chain[1, 2*nstrain+3+time+12+ndept+(1:nstrain)], Bits=Bits, a_k=MC_chain[1, 2*nstrain+3+time+12+ndept+nstrain+(1:nstrain)], Model=Model,Q_r=Q_r,Q_s = Q_s,Q_u=Q_u)
  likelihoodcurrent<- Allquantities$loglike
  priorcurrentRcomps<- randomwalk2(MC_chain[1, 2*nstrain+3+(1:time)], MC_chain[1, 2*nstrain+1])
  priorcurrentScomps<- seasonalComp2(MC_chain[1, 2*nstrain+3+time+(1:12)], MC_chain[1, 2*nstrain+2], RW1PrecMat)
  priorcurrentUcomps<- logIGMRF1(MC_chain[1, 2*nstrain+3+time+12+(1:ndept)], MC_chain[1, 2*nstrain+3], R, rankdef)

  grad_current <- list(grad_r=as.numeric(Allquantities$grad_r), grad_s=as.numeric(Allquantities$grad_s), grad_u=as.numeric(Allquantities$grad_u), cov_r=Allquantities$cov_r, cov_s=Allquantities$cov_s)

  for (i in 2:num_iteration) {

#    MC_chain[i,2*nstrain+1]<- rgamma(1, shape = 1 + (time-2)/2, rate = 0.0001 + (t(MC_chain[i-1, 2*nstrain+3+(1:time)]) %*% RW2PrecMat %*% MC_chain[i-1, 2*nstrain+3+(1:time)])/2)
#    MC_chain[i,2*nstrain+2]<- rgamma(1, shape = 1 + 11/2, rate = 0.001 + (t(MC_chain[i-1, 2*nstrain+3+time+(1:12)]) %*% RW1PrecMat %*% MC_chain[i-1, 2*nstrain+3+time+(1:12)])/2)
#    MC_chain[i,2*nstrain+3]<- rgamma(1, shape = 1 + (ndept-1)/2, rate = 0.01 + (t(MC_chain[i-1, 2*nstrain+3+time+12+(1:ndept)]) %*% R %*% MC_chain[i-1, 2*nstrain+3+time+12+(1:ndept)])/2)

#    Q_r<- MC_chain[i,2*nstrain+1] * RW2PrecMat
#    Q_s<- MC_chain[i,2*nstrain+2] * RW1PrecMat
#    Q_u<- MC_chain[i,2*nstrain+3] * R

#    current_r <- MC_chain[i-1, 2*nstrain+3+(1:time)]
#    current_s <- MC_chain[i-1, 2*nstrain+3+time+(1:12)]
#    current_u <- MC_chain[i-1, 2*nstrain+3+time+12+(1:ndept)]

    #Update s
#    Mmatcs<- as.numeric(grad_current$cov_s %*% grad_current$grad_s)

#    eps_s <- rnorm(12)
#    proposedScomps <- as.numeric(MC_chain[i-1, 2*nstrain+3+time+(1:12)] + 0.5 * step_sizes$s^2 * Mmatcs + step_sizes$s * chol(grad_current$cov_s) %*% eps_s)
#    proposedScomps<- proposedScomps - mean(proposedScomps)

#    Gamma_lists<- BuildGamma_list(MC_chain[i-1,1:(2*nstrain)])
#    copulaTPM<- JointTransitionMatrix_copula_per_strain(Gamma_lists, MC_chain[i-1, ncol(MC_chain)])

#    Allquantities<- copulagradmultstrainLoglikelihood2_cpp(y=y, e_it=e_it, nstrain=nstrain,  r=current_r, s=proposedScomps, u=current_u, Gamma=copulaTPM, B=MC_chain[i-1, 2*nstrain+3+time+12+ndept+(1:nstrain)], Bits=Bits, a_k=MC_chain[i-1, 2*nstrain+3+time+12+ndept+nstrain+(1:nstrain)], Model=Model,Q_r=Q_r,Q_s = Q_s,Q_u=Q_u)
#    grad_proposed <- list(grad_r=as.numeric(Allquantities$grad_r), grad_s=as.numeric(Allquantities$grad_s), grad_u=as.numeric(Allquantities$grad_u), cov_r=Allquantities$cov_r, cov_s=Allquantities$cov_s)

#    Mmatps<- as.numeric(grad_proposed$cov_s %*% grad_proposed$grad_s)

#    likelihoodproposed<- Allquantities$loglike
#    q_prop <- mvnfast::dmvn(proposedScomps, mu = MC_chain[i-1, 2*nstrain+3+time+(1:12)] + 0.5 * step_sizes$s^2 * Mmatcs, sigma = grad_current$cov_s * step_sizes$s^2, log = TRUE)
#    q_curr <- mvnfast::dmvn(MC_chain[i-1, 2*nstrain+3+time+(1:12)], mu = proposedScomps + 0.5 * step_sizes$s^2 * Mmatps, sigma = grad_proposed$cov_s * step_sizes$s^2, log = TRUE)
#    priorproposedScomps<- seasonalComp2(proposedScomps, MC_chain[i, 2*nstrain+2], RW1PrecMat)

#    log_alpha_s <- likelihoodproposed + priorproposedScomps + q_curr - likelihoodcurrent - priorcurrentScomps - q_prop

    MC_chain[i, 2*nstrain+3+time+(1:12)]<- copulamultmod1nstrain2[["s"]]

#    if (is.finite(log_alpha_s) && log(runif(1)) < log_alpha_s){
#      MC_chain[i, 2*nstrain+3+time+(1:12)]<- proposedScomps
#      likelihoodcurrent<- likelihoodproposed
#      priorcurrentScomps<- priorproposedScomps
#      grad_current<- grad_proposed
#    }else{
#      MC_chain[i, 2*nstrain+3+time+(1:12)]<- MC_chain[i-1, 2*nstrain+3+time+(1:12)]
#    }

    #Update r
#    eps_r <- rnorm(time)
#    Mmatrc<- as.numeric(grad_current$cov_r %*% grad_current$grad_r)

#    proposedRcomps <- as.numeric(current_r + 0.5 * step_sizes$r^2 * Mmatrc + step_sizes$r * chol(grad_current$cov_r) %*% eps_r)
#    proposedRcomps<- proposedRcomps - mean(proposedRcomps)

#    Allquantities<- copulagradmultstrainLoglikelihood2_cpp(y=y, e_it=e_it, nstrain=nstrain,  r=proposedRcomps, s=MC_chain[i, 2*nstrain+3+time+(1:12)], u=current_u, Gamma=copulaTPM, B=MC_chain[i-1, 2*nstrain+3+time+12+ndept+(1:nstrain)], Bits=Bits, a_k=MC_chain[i-1, 2*nstrain+3+time+12+ndept+nstrain+(1:nstrain)], Model=Model,Q_r=Q_r,Q_s = Q_s,Q_u=Q_u)
#    grad_proposed <- list(grad_r=as.numeric(Allquantities$grad_r), grad_s=as.numeric(Allquantities$grad_s), grad_u=as.numeric(Allquantities$grad_u), cov_r=Allquantities$cov_r, cov_s=Allquantities$cov_s)

#    Mmatrp<- as.numeric(grad_proposed$cov_r %*% grad_proposed$grad_r)

#    q_prop <- mvnfast::dmvn(proposedRcomps, mu = current_r + 0.5 * step_sizes$r^2 * Mmatrc, sigma = grad_current$cov_r * step_sizes$r^2, log = TRUE)
#    q_curr <- mvnfast::dmvn(current_r, mu = proposedRcomps + 0.5 * step_sizes$r^2 * Mmatrp, sigma = grad_proposed$cov_r * step_sizes$r^2, log = TRUE)

#    likelihoodproposed<- Allquantities$loglike
#    priorproposedRcomps <- randomwalk2(proposedRcomps, MC_chain[i, 2*nstrain+1])

#    log_alpha_r <- likelihoodproposed + priorproposedRcomps + q_curr - likelihoodcurrent - priorcurrentRcomps - q_prop

    MC_chain[i, 2*nstrain+3+(1:time)] <- copulamultmod1nstrain2[["r"]]

#    if (is.finite(log_alpha_r) && log(runif(1)) < log_alpha_r){
#      MC_chain[i, 2*nstrain+3+(1:time)] <- proposedRcomps
#      likelihoodcurrent<- likelihoodproposed
#      priorcurrentRcomps<- priorproposedRcomps
#      grad_current<- grad_proposed
#    }else{
#      MC_chain[i, 2*nstrain+3+(1:time)]<- MC_chain[i-1, 2*nstrain+3+(1:time)]
 #   }

    #Update u
#    eps_u <- rnorm(ndept)

#    proposedUcomps <- as.numeric(MC_chain[i-1, 2*nstrain+3+time+12+(1:ndept)] + 0.5 * step_sizes$u^2 * grad_current$grad_u + step_sizes$u * eps_u)
#    proposedUcomps<- proposedUcomps - mean(proposedUcomps)

#    Allquantities<- copulagradmultstrainLoglikelihood2_cpp(y=y, e_it=e_it, nstrain=nstrain,  r=MC_chain[i, 2*nstrain+3+(1:time)], s=MC_chain[i, 2*nstrain+3+time+(1:12)], u=proposedUcomps, Gamma=copulaTPM, B=MC_chain[i-1, 2*nstrain+3+time+12+ndept+(1:nstrain)], Bits=Bits, a_k=MC_chain[i-1, 2*nstrain+3+time+12+ndept+nstrain+(1:nstrain)], Model=Model,Q_r=Q_r,Q_s = Q_s,Q_u=Q_u)
#    grad_proposed <- list(grad_r=as.numeric(Allquantities$grad_r), grad_s=as.numeric(Allquantities$grad_s), grad_u=as.numeric(Allquantities$grad_u), cov_r=Allquantities$cov_r, cov_s=Allquantities$cov_s)

#    likelihoodproposed<- Allquantities$loglike

#    q_prop <- sum(dnorm(proposedUcomps, mean = MC_chain[i-1, 2*nstrain+3+time+12+(1:ndept)] + 0.5 * step_sizes$u^2 * grad_current$grad_u, sd = step_sizes$u, log = TRUE))
#    q_curr <- sum(dnorm(MC_chain[i-1, 2*nstrain+3+time+12+(1:ndept)], mean = proposedUcomps + 0.5 * step_sizes$u^2 * grad_proposed$grad_u, sd = step_sizes$u, log = TRUE))

#    priorproposedUcomps<- logIGMRF1(proposedUcomps, MC_chain[i, 2*nstrain+3], R, rankdef)
#    priorcurrentUcomps<- logIGMRF1(MC_chain[i-1, 2*nstrain+3+time+12+(1:ndept)], MC_chain[i, 2*nstrain+3], R, rankdef)

 #   log_alpha_u <- likelihoodproposed + priorproposedUcomps + q_curr - likelihoodcurrent - priorcurrentUcomps - q_prop

    MC_chain[i, 2*nstrain+3+time+12+(1:ndept)]<- copulamultmod1nstrain2[["u"]]

#    if (is.finite(log_alpha_u) && log(runif(1)) < log_alpha_u){
#      MC_chain[i, 2*nstrain+3+time+12+(1:ndept)]<- proposedUcomps
#      likelihoodcurrent<- likelihoodproposed
      #priorcurrentUcomps<- priorproposedUcomps
#      grad_current<- grad_proposed
#    }else{
#      MC_chain[i, 2*nstrain+3+time+12+(1:ndept)]<- MC_chain[i-1, 2*nstrain+3+time+12+(1:ndept)]
#    }

#    if(Model == 0){
#      MC_chain[i, 2*nstrain+3+time+12+ndept+(1:nstrain)]<-  MC_chain[i-1, 2*nstrain+3+time+12+ndept+(1:nstrain)]
#      MC_chain[i, 1:(2*nstrain)]<- MC_chain[i-1, 1:(2*nstrain)]
#      MC_chain[i, ncol(MC_chain)]<- MC_chain[i-1, ncol(MC_chain)]
#    }else{
#      proposedB <- abs(rnorm(nstrain, mean = MC_chain[i-1, 2*nstrain+3+time+12+ndept+(1:nstrain)], sd = rep(sdBs, nstrain)))
#      priorcurrentB<- sum(dgamma(MC_chain[i-1, 2*nstrain+3+time+12+ndept+(1:nstrain)], shape = rep(2, nstrain), rate = rep(2,nstrain), log=TRUE))
#      priorproposedB<- sum(dgamma(proposedB, shape = rep(2, nstrain), rate = rep(2, nstrain), log=TRUE))

#      Allquantities<- copulagradmultstrainLoglikelihood2_cpp(y=y, e_it=e_it, nstrain=nstrain,  r=MC_chain[i, 2*nstrain+3+(1:time)], s=MC_chain[i, 2*nstrain+3+time+(1:12)], u=MC_chain[i, 2*nstrain+3+time+12+(1:ndept)], Gamma=copulaTPM, B=proposedB, Bits=Bits, a_k=MC_chain[i-1, 2*nstrain+3+time+12+ndept+nstrain+(1:nstrain)], Model=Model,Q_r=Q_r,Q_s = Q_s,Q_u=Q_u)
#      grad_proposed <- list(grad_r=as.numeric(Allquantities$grad_r), grad_s=as.numeric(Allquantities$grad_s), grad_u=as.numeric(Allquantities$grad_u), cov_r=Allquantities$cov_r, cov_s=Allquantities$cov_s)

#      likelihoodproposed<- Allquantities$loglike

#      mh.ratio<- exp(likelihoodproposed + priorproposedB
#                     - likelihoodcurrent - priorcurrentB)

      #print(paste("mh.ratioB = ", mh.ratio))

      MC_chain[i, 2*nstrain+3+time+12+ndept+(1:nstrain)]<- copulamultmod1nstrain2[["B"]]

#      if(!is.na(mh.ratio) && runif(1) < mh.ratio){
#        MC_chain[i, 2*nstrain+3+time+12+ndept+(1:nstrain)]<- proposedB
#        likelihoodcurrent<- likelihoodproposed
#        grad_current<- grad_proposed
#      }
#      else{
#        MC_chain[i, 2*nstrain+3+time+12+ndept+(1:nstrain)]<- MC_chain[i-1, 2*nstrain+3+time+12+ndept+(1:nstrain)]
 #     }

#      proposedGs<- abs(rnorm(2*nstrain,mean=MC_chain[i-1,1:(2*nstrain)], sd=rep(sdGs, 2*nstrain)))
#      proposedGs<- ifelse(proposedGs<1, proposedGs, 2-proposedGs)

#      priorcurrentGs<- sum(dbeta(MC_chain[i-1,1:(2*nstrain)], shape1 = rep(2,2*nstrain), shape2 = rep(2,2*nstrain), log=TRUE))
#      priorproposedGs<- sum(dbeta(proposedGs, shape1 = rep(2,2*nstrain), shape2 = rep(2,2*nstrain), log=TRUE))

#      Gamma_lists<- BuildGamma_list(proposedGs)
#      copulaTPM<- JointTransitionMatrix_copula_per_strain(Gamma_lists, MC_chain[i-1, ncol(MC_chain)])

#      Allquantities<- copulagradmultstrainLoglikelihood2_cpp(y=y, e_it=e_it, nstrain=nstrain,  r=MC_chain[i, 2*nstrain+3+(1:time)], s=MC_chain[i, 2*nstrain+3+time+(1:12)], u=MC_chain[i, 2*nstrain+3+time+12+(1:ndept)], Gamma=copulaTPM, B=MC_chain[i, 2*nstrain+3+time+12+ndept+(1:nstrain)], Bits=Bits, a_k=MC_chain[i-1, 2*nstrain+3+time+12+ndept+nstrain+(1:nstrain)], Model=Model,Q_r=Q_r,Q_s = Q_s,Q_u=Q_u)
#      grad_proposed <- list(grad_r=as.numeric(Allquantities$grad_r), grad_s=as.numeric(Allquantities$grad_s), grad_u=as.numeric(Allquantities$grad_u), cov_r=Allquantities$cov_r, cov_s=Allquantities$cov_s)

#      likelihoodproposed<- Allquantities$loglike

#      mh.ratio<- exp(likelihoodproposed + priorproposedGs
#                     - likelihoodcurrent - priorcurrentGs)

      #print(mh.ratio)

      MC_chain[i, 1:(2*nstrain)]<- copulamultmod1nstrain2[["T.prob"]]

#      if(!is.na(mh.ratio) && runif(1) < mh.ratio){
#        MC_chain[i, 1:(2*nstrain)]<- proposedGs
#        likelihoodcurrent<- likelihoodproposed
#        grad_current<- grad_proposed
#      }
#      else{
#        MC_chain[i, 1:(2*nstrain)]<- MC_chain[i-1,1:(2*nstrain)]
 #     }

      proposedcopPs<- rnorm(n_copparams,mean=log(MC_chain[i-1, ncol(MC_chain)]), sd=rep(0.7, n_copparams))

#      proposedcopPs<- rnorm(1, mean=MC_chain[i-1, ncol(MC_chain)], sd=0.7)
#      if(proposedcopPs< -1 || proposedcopPs> 1) proposedcopPs = MC_chain[i-1, ncol(MC_chain)]

      priorcurrentcopPs<- dlnorm(MC_chain[i-1, ncol(MC_chain)], meanlog = exp(proposedcopPs), sdlog = 500, log=TRUE)
      priorproposedcopPs<- dlnorm(exp(proposedcopPs), meanlog = MC_chain[i-1, ncol(MC_chain)], sdlog = 500, log=TRUE)

#      priorcurrentcopPs<- sum(dunif(MC_chain[i-1, (ncol(MC_chain)-n_copparams):ncol(MC_chain)], min = -1, max = 1, log=TRUE))
#      priorproposedcopPs<- sum(dunif(proposedcopPs, min = -1, max = 1, log=TRUE))

      proposalcurrentcop<- dnorm(log(MC_chain[i-1, ncol(MC_chain)]), mean=proposedcopPs, sd=rep(0.7, n_copparams), log = T) - log(MC_chain[i-1, ncol(MC_chain)])
      proposalproposedcop<- dnorm(proposedcopPs, mean=log(MC_chain[i-1, ncol(MC_chain)]), sd=rep(0.7, n_copparams), log = T) - proposedcopPs

      Gamma_lists<- BuildGamma_list(MC_chain[i, 1:(2*nstrain)])
      copulaTPM<- JointTransitionMatrix_copula_per_strain(Gamma_lists, exp(proposedcopPs))
 #     copulaTPM<- JointTransitionMatrix_copula_per_strain(Gamma_lists, proposedcopPs)

      Allquantities<- copulagradmultstrainLoglikelihood2_cpp(y=y, e_it=e_it, nstrain=nstrain,  r=MC_chain[i, 2*nstrain+3+(1:time)], s=MC_chain[i, 2*nstrain+3+time+(1:12)], u=MC_chain[i, 2*nstrain+3+time+12+(1:ndept)], Gamma=copulaTPM, B=MC_chain[i, 2*nstrain+3+time+12+ndept+(1:nstrain)], Bits=Bits, a_k=MC_chain[i-1, 2*nstrain+3+time+12+ndept+nstrain+(1:nstrain)], Model=Model,Q_r=Q_r,Q_s = Q_s,Q_u=Q_u)
#      grad_proposed <- list(grad_r=as.numeric(Allquantities$grad_r), grad_s=as.numeric(Allquantities$grad_s), grad_u=as.numeric(Allquantities$grad_u), cov_r=Allquantities$cov_r, cov_s=Allquantities$cov_s)

      likelihoodproposed<- Allquantities$loglike

      mh.ratio<- exp(likelihoodproposed - likelihoodcurrent + proposalcurrentcop - proposalproposedcop
      + priorproposedcopPs
      - priorcurrentcopPs)

      #print(mh.ratio)

      if(!is.na(mh.ratio) && runif(1) < mh.ratio){
        MC_chain[i, ncol(MC_chain)]<- exp(proposedcopPs)
#        MC_chain[i, ncol(MC_chain)]<- proposedcopPs
        likelihoodcurrent<- likelihoodproposed
 #       grad_current<- grad_proposed
      }
      else{
        MC_chain[i, ncol(MC_chain)]<- MC_chain[i-1, ncol(MC_chain)]
      }
#    }

    #Gibbs A_k's update
    MC_chain[i, 2*nstrain+3+time+12+ndept+nstrain+(1:nstrain)]<- copulamultmod1nstrain2[["a_k"]]
    #MC_chain[i, 2*nstrain+3+time+12+ndept+nstrain+(1:nstrain)]<- log(rgamma(nstrain, shape = 0.01+SumYk_vec, rate = as.numeric(Allquantities$poisMean4GibbsUpdate) + 0.01/exp(-15)))

    if(i %% 1000 == 0) cat("Iteration:", i, "\n")
  }
  colnames(MC_chain) <- paste(c(paste0(rep(c("G12", "G21"), nstrain), "Strain", rep(1:nstrain,each=2)), "kappa_r", "kappa_s", "kappa_u", paste("r", 1:time, sep=""), paste("s", 1:12, sep=""), paste("u", 1:ndept, sep=""), paste("B", 1:nstrain, sep=""), paste("a_k", 1:nstrain, sep=""), paste("copulaParam", 1:n_copparams, sep="")))
  MC_chain<- as.data.frame(MC_chain)
  end_time <- Sys.time()
  time_taken<- end_time - start_time
  print(time_taken)
  return(MC_chain)
}


#Riemann Manifold Langevin updates
FinalCPPmultMMALAInference<- function(y, e_it, Model, adjmat, step_sizes, num_iteration = 15000, sdGs=0.05, sdBs=0.03, sdAs=0.03, sdCops=0.03) {
  start_time <- Sys.time()
  ndept <- nrow(e_it)
  time <- ncol(e_it)
  nstrain<- dim(y)[3]
  Bits<- encodeBits(nstrain)
  n_copParams<- (nstrain*(nstrain-1))/2

  R<- -1 * adjmat
  diag(R)<- -rowSums(R, na.rm = T)
  rankdef<- nrow(R)-qr(R)$rank

  RW1PrecMat<- matrix(0, nrow=12, ncol=12)
  RW1PrecMat[1, ]<- c(2,-1, rep(0, 12-3), -1)
  RW1PrecMat[2, ]<- c(-1,2,-1, rep(0, 12-3))
  RW1PrecMat[3, ]<- c(0, -1,2,-1, rep(0, 12-4))
  RW1PrecMat[(12-1), ]<- c(rep(0, 12-3), -1,2,-1)
  RW1PrecMat[12, ]<- c(-1, rep(0, 12-3), -1, 2)
  for(i in 3:(12-3)){
    RW1PrecMat[i+1, ((i):(i+2))]<- c(-1,2,-1)
  }

  RW2PrecMat<- matrix(0, nrow=time, ncol=time)
  RW2PrecMat[1,(1:3)]<- c(1,-2,1)
  RW2PrecMat[2,(1:4)]<- c(-2,5,-4,1)
  RW2PrecMat[3,(1:5)]<- c(1,-4,6,-4,1)
  RW2PrecMat[(time-1),((time-3):time)]<- c(1,-4,5,-2)
  RW2PrecMat[time,((time-2):time)]<- c(1,-2,1)
  for(i in 3:(time-3)){
    RW2PrecMat[i+1, ((i-1):(i+3))]<- c(1,-4,6,-4,1)
  }

  SumYk_vec<- numeric(nstrain)
  SumYk_vec[1]<- sum(y[,,1])
  sumY<- y[,,1]
  for(k in 2:nstrain){
    sumY<- sumY + y[,,k]
    SumYk_vec[k]<- sum(y[,,k])
  }

  crudeResults<- DetectOutbreaks:::crudeEst(sumY, e_it)
  crudeR<- crudeResults[[1]] - mean(crudeResults[[1]])
  crudeS<- crudeResults[[2]]
  crudeU<- crudeResults[[3]]
  crudeU<- ifelse(is.nan(crudeU), mean(crudeU[is.finite(crudeU)]), crudeU)
  crudeblock<- floor(time/12)
  crudeblock<- ((crudeblock*12)-11):(crudeblock*12)

  MC_chain<- matrix(NA, nrow=num_iteration, ncol=5+time+12+ndept+nstrain+nstrain+n_copParams+1)
  initG12<- runif(1)
  initG21<- runif(1)
  initstateD<- stationarydist(G(initG12, initG21))[2]
  MC_chain[1,]<- c(initG12, initG21, 1/var(crudeR), 1/var(crudeS), 1/var(crudeU), crudeR, crudeS[crudeblock-12], crudeU, rep(0, nstrain), rep(mean(crudeResults[[1]]), nstrain), rep(0, n_copParams), initstateD)

  #flag missing data for inference
  y<- ifelse(is.na(y), -1, y)
  yflat<- as.numeric(aperm(y, c(2,1,3)))

  Q_r<- MC_chain[1,3] * RW2PrecMat
  Q_s<- MC_chain[1,4] * RW1PrecMat
  Q_u<- MC_chain[1,5] * R

  #Compute gradients
  JointTPM<- JointTransitionMatrix_copula_cpp(G(MC_chain[1,1],MC_chain[1,2]), nstrain, MC_chain[1,(ncol(MC_chain)-n_copParams):(ncol(MC_chain)-1)])
  JointTPM<- ifelse(JointTPM<=0,1e-6,JointTPM)
  JointTPM<- ifelse(JointTPM>=1,1-1e-6,JointTPM)

  Allquantities<- FFBSgradmultstrainLoglikelihood2_cpp(y=y, e_it=e_it, nstrain=nstrain,  r=MC_chain[1, 5+(1:time)], s=MC_chain[1, 5+time+(1:12)], u=MC_chain[1, 5+time+12+(1:ndept)], Gamma=JointTPM, B=MC_chain[1, 5+time+12+ndept+(1:nstrain)], Bits=Bits, a_k=MC_chain[1, 5+time+12+ndept+nstrain+(1:nstrain)], Model=Model,Q_r=Q_r,Q_s = Q_s,Q_u=Q_u)
  likelihoodcurrent<- Allquantities$loglike
  priorcurrentRcomps<- randomwalk2(MC_chain[1, 5+(1:time)], MC_chain[1, 3])
  priorcurrentScomps<- seasonalComp2(MC_chain[1, 5+time+(1:12)], MC_chain[1, 4], RW1PrecMat)
  priorcurrentUcomps<- logIGMRF1(MC_chain[1, 5+time+12+(1:ndept)], MC_chain[1, 5], R, rankdef)

  grad_current <- list(grad_r=as.numeric(Allquantities$grad_r), grad_s=as.numeric(Allquantities$grad_s), grad_u=as.numeric(Allquantities$grad_u), cov_r=Allquantities$cov_r, cov_s=Allquantities$cov_s)

  for (i in 2:num_iteration) {

    MC_chain[i,3]<- rgamma(1, shape = 1 + (time-2)/2, rate = 0.0001 + (t(MC_chain[i-1, 5+(1:time)]) %*% RW2PrecMat %*% MC_chain[i-1, 5+(1:time)])/2)
    MC_chain[i,4]<- rgamma(1, shape = 1 + 11/2, rate = 0.001 + (t(MC_chain[i-1, 5+time+(1:12)]) %*% RW1PrecMat %*% MC_chain[i-1, 5+time+(1:12)])/2)
    MC_chain[i,5]<- rgamma(1, shape = 1 + (ndept-1)/2, rate = 0.01 + (t(MC_chain[i-1, 5+time+12+(1:ndept)]) %*% R %*% MC_chain[i-1, 5+time+12+(1:ndept)])/2)

    Q_r<- MC_chain[i,3] * RW2PrecMat
    Q_s<- MC_chain[i,4] * RW1PrecMat
    Q_u<- MC_chain[i,5] * R

    current_r <- MC_chain[i-1, 5+(1:time)]
    current_s <- MC_chain[i-1, 5+time+(1:12)]
    current_u <- MC_chain[i-1, 5+time+12+(1:ndept)]

    #Update s
    Mmatcs<- as.numeric(grad_current$cov_s %*% grad_current$grad_s)

    eps_s <- rnorm(12)
    proposedScomps <- as.numeric(MC_chain[i-1, 5+time+(1:12)] + 0.5 * step_sizes$s^2 * Mmatcs + step_sizes$s * chol(grad_current$cov_s) %*% eps_s)
    proposedScomps<- proposedScomps - mean(proposedScomps)

    Allquantities<- FFBSgradmultstrainLoglikelihood2_cpp(y=y, e_it=e_it, nstrain=nstrain,  r=current_r, s=proposedScomps, u=current_u, Gamma=JointTPM, B=MC_chain[i-1, 5+time+12+ndept+(1:nstrain)], Bits=Bits, a_k=MC_chain[i-1, 5+time+12+ndept+nstrain+(1:nstrain)], Model=Model,Q_r=Q_r,Q_s = Q_s,Q_u=Q_u)
    grad_proposed <- list(grad_r=as.numeric(Allquantities$grad_r), grad_s=as.numeric(Allquantities$grad_s), grad_u=as.numeric(Allquantities$grad_u), cov_r=Allquantities$cov_r, cov_s=Allquantities$cov_s)

    Mmatps<- as.numeric(grad_proposed$cov_s %*% grad_proposed$grad_s)

    likelihoodproposed<- Allquantities$loglike
    q_prop <- mvnfast::dmvn(proposedScomps, mu = MC_chain[i-1, 5+time+(1:12)] + 0.5 * step_sizes$s^2 * Mmatcs, sigma = grad_current$cov_s * step_sizes$s^2, log = TRUE)
    q_curr <- mvnfast::dmvn(MC_chain[i-1, 5+time+(1:12)], mu = proposedScomps + 0.5 * step_sizes$s^2 * Mmatps, sigma = grad_proposed$cov_s * step_sizes$s^2, log = TRUE)
    priorproposedScomps<- seasonalComp2(proposedScomps, MC_chain[i, 4], RW1PrecMat)

    log_alpha_s <- likelihoodproposed + priorproposedScomps + q_curr - likelihoodcurrent - priorcurrentScomps - q_prop
    if (is.finite(log_alpha_s) && log(runif(1)) < log_alpha_s){
      MC_chain[i, 5+time+(1:12)]<- proposedScomps
      likelihoodcurrent<- likelihoodproposed
      priorcurrentScomps<- priorproposedScomps
      grad_current<- grad_proposed
    }else{
      MC_chain[i, 5+time+(1:12)]<- MC_chain[i-1, 5+time+(1:12)]
    }

    #Update r
    eps_r <- rnorm(time)
    Mmatrc<- as.numeric(grad_current$cov_r %*% grad_current$grad_r)

    proposedRcomps <- as.numeric(current_r + 0.5 * step_sizes$r^2 * Mmatrc + step_sizes$r * chol(grad_current$cov_r) %*% eps_r)
    proposedRcomps<- proposedRcomps - mean(proposedRcomps)

    Allquantities<- FFBSgradmultstrainLoglikelihood2_cpp(y=y, e_it=e_it, nstrain=nstrain,  r=proposedRcomps, s=MC_chain[i, 5+time+(1:12)], u=current_u, Gamma=JointTPM, B=MC_chain[i-1, 5+time+12+ndept+(1:nstrain)], Bits=Bits, a_k=MC_chain[i-1, 5+time+12+ndept+nstrain+(1:nstrain)], Model=Model,Q_r=Q_r,Q_s = Q_s,Q_u=Q_u)
    grad_proposed <- list(grad_r=as.numeric(Allquantities$grad_r), grad_s=as.numeric(Allquantities$grad_s), grad_u=as.numeric(Allquantities$grad_u), cov_r=Allquantities$cov_r, cov_s=Allquantities$cov_s)

    Mmatrp<- as.numeric(grad_proposed$cov_r %*% grad_proposed$grad_r)

    q_prop <- mvnfast::dmvn(proposedRcomps, mu = current_r + 0.5 * step_sizes$r^2 * Mmatrc, sigma = grad_current$cov_r * step_sizes$r^2, log = TRUE)
    q_curr <- mvnfast::dmvn(current_r, mu = proposedRcomps + 0.5 * step_sizes$r^2 * Mmatrp, sigma = grad_proposed$cov_r * step_sizes$r^2, log = TRUE)

    likelihoodproposed<- Allquantities$loglike
    priorproposedRcomps <- randomwalk2(proposedRcomps, MC_chain[i, 3])

    log_alpha_r <- likelihoodproposed + priorproposedRcomps + q_curr - likelihoodcurrent - priorcurrentRcomps - q_prop

    if (is.finite(log_alpha_r) && log(runif(1)) < log_alpha_r){
      MC_chain[i, 5+(1:time)] <- proposedRcomps
      likelihoodcurrent<- likelihoodproposed
      priorcurrentRcomps<- priorproposedRcomps
      grad_current<- grad_proposed
    }else{
      MC_chain[i, 5+(1:time)]<- MC_chain[i-1, 5+(1:time)]
    }

    #Update u
    eps_u <- rnorm(ndept)

    proposedUcomps <- as.numeric(MC_chain[i-1, 5+time+12+(1:ndept)] + 0.5 * step_sizes$u^2 * grad_current$grad_u + step_sizes$u * eps_u)
    proposedUcomps<- proposedUcomps - mean(proposedUcomps)

    Allquantities<- FFBSgradmultstrainLoglikelihood2_cpp(y=y, e_it=e_it, nstrain=nstrain,  r=MC_chain[i, 5+(1:time)], s=MC_chain[i, 5+time+(1:12)], u=proposedUcomps, Gamma=JointTPM, B=MC_chain[i-1, 5+time+12+ndept+(1:nstrain)], Bits=Bits, a_k=MC_chain[i-1, 5+time+12+ndept+nstrain+(1:nstrain)], Model=Model,Q_r=Q_r,Q_s = Q_s,Q_u=Q_u)
    grad_proposed <- list(grad_r=as.numeric(Allquantities$grad_r), grad_s=as.numeric(Allquantities$grad_s), grad_u=as.numeric(Allquantities$grad_u), cov_r=Allquantities$cov_r, cov_s=Allquantities$cov_s)

    likelihoodproposed<- Allquantities$loglike

    q_prop <- sum(dnorm(proposedUcomps, mean = MC_chain[i-1, 5+time+12+(1:ndept)] + 0.5 * step_sizes$u^2 * grad_current$grad_u, sd = step_sizes$u, log = TRUE))
    q_curr <- sum(dnorm(MC_chain[i-1, 5+time+12+(1:ndept)], mean = proposedUcomps + 0.5 * step_sizes$u^2 * grad_proposed$grad_u, sd = step_sizes$u, log = TRUE))

    priorproposedUcomps<- logIGMRF1(proposedUcomps, MC_chain[i, 5], R, rankdef)
    priorcurrentUcomps<- logIGMRF1(MC_chain[i-1, 5+time+12+(1:ndept)], MC_chain[i, 5], R, rankdef)

    log_alpha_u <- likelihoodproposed + priorproposedUcomps + q_curr - likelihoodcurrent - priorcurrentUcomps - q_prop
    if (is.finite(log_alpha_u) && log(runif(1)) < log_alpha_u){
      MC_chain[i, 5+time+12+(1:ndept)]<- proposedUcomps
      likelihoodcurrent<- likelihoodproposed
      #priorcurrentUcomps<- priorproposedUcomps
      grad_current<- grad_proposed
    }else{
      MC_chain[i, 5+time+12+(1:ndept)]<- MC_chain[i-1, 5+time+12+(1:ndept)]
    }

    if(Model == 0){
      MC_chain[i, 5+time+12+ndept+(1:nstrain)]<-  MC_chain[i-1, 5+time+12+ndept+(1:nstrain)]
      MC_chain[i, 1:2]<- MC_chain[i-1, 1:2]
    }else{
      proposedB <- abs(rnorm(nstrain, mean = MC_chain[i-1, 5+time+12+ndept+(1:nstrain)], sd = rep(sdBs, nstrain)))
      priorcurrentB<- sum(dgamma(MC_chain[i-1, 5+time+12+ndept+(1:nstrain)], shape = rep(2, nstrain), rate = rep(2,nstrain), log=TRUE))
      priorproposedB<- sum(dgamma(proposedB, shape = rep(2, nstrain), rate = rep(2, nstrain), log=TRUE))

      Allquantities<- FFBSgradmultstrainLoglikelihood2_cpp(y=y, e_it=e_it, nstrain=nstrain,  r=MC_chain[i, 5+(1:time)], s=MC_chain[i, 5+time+(1:12)], u=MC_chain[i, 5+time+12+(1:ndept)], Gamma=JointTPM, B=proposedB, Bits=Bits, a_k=MC_chain[i-1, 5+time+12+ndept+nstrain+(1:nstrain)], Model=Model,Q_r=Q_r,Q_s = Q_s,Q_u=Q_u)
      grad_proposed <- list(grad_r=as.numeric(Allquantities$grad_r), grad_s=as.numeric(Allquantities$grad_s), grad_u=as.numeric(Allquantities$grad_u), cov_r=Allquantities$cov_r, cov_s=Allquantities$cov_s)

      likelihoodproposed<- Allquantities$loglike

      mh.ratio<- exp(likelihoodproposed + priorproposedB
                     - likelihoodcurrent - priorcurrentB)

      #print(paste("mh.ratioB = ", mh.ratio))

      if(!is.na(mh.ratio) && runif(1) < mh.ratio){
        MC_chain[i, 5+time+12+ndept+(1:nstrain)]<- proposedB
        likelihoodcurrent<- likelihoodproposed
        grad_current<- grad_proposed
      }
      else{
        MC_chain[i, 5+time+12+ndept+(1:nstrain)]<- MC_chain[i-1, 5+time+12+ndept+(1:nstrain)]
      }

      proposedGs<- abs(rnorm(2,mean=c(MC_chain[i-1,1], MC_chain[i-1, 2]), sd=c(sdGs, sdGs)))
      if(proposedGs[1]>1) proposedGs[1]=2-proposedGs[1]
      if(proposedGs[2]>1) proposedGs[2]=2-proposedGs[2]

      priorcurrentGs<- sum(dbeta(MC_chain[i-1,1:2], shape1 = c(2,2), shape2 = c(2,2), log=TRUE))
      priorproposedGs<- sum(dbeta(proposedGs, shape1 = c(2,2), shape2 = c(2,2), log=TRUE))

      proposedcopPs<- rnorm(n_copParams, mean=MC_chain[i-1, (ncol(MC_chain)-n_copParams):(ncol(MC_chain)-1)], sd=rep(sdCops, n_copParams))

      if(nstrain>3){
        JointTPM1<- ParallelJointTransitionMatrix_copula_cpp2(G(proposedGs[1],proposedGs[2]), nstrain, proposedcopPs)
      }else{
        JointTPM1<- JointTransitionMatrix_copula_cpp(G(proposedGs[1],proposedGs[2]), nstrain, proposedcopPs)
      }
      JointTPM1<- ifelse(JointTPM1<=0,1e-6,JointTPM1)
      JointTPM1<- ifelse(JointTPM1>=1,1-1e-6,JointTPM1)

      if(any(!is.finite(JointTPM1))){
        MC_chain[i, (ncol(MC_chain)-n_copParams):(ncol(MC_chain)-1)]<- MC_chain[i-1, (ncol(MC_chain)-n_copParams):(ncol(MC_chain)-1)]
        MC_chain[i, 1:2]<- MC_chain[i-1,1:2]
      }else{

      Allquantities<- FFBSgradmultstrainLoglikelihood2_cpp(y=y, e_it=e_it, nstrain=nstrain,  r=MC_chain[i, 5+(1:time)], s=MC_chain[i, 5+time+(1:12)], u=MC_chain[i, 5+time+12+(1:ndept)], Gamma=JointTPM1, B=MC_chain[i, 5+time+12+ndept+(1:nstrain)], Bits=Bits, a_k=MC_chain[i-1, 5+time+12+ndept+nstrain+(1:nstrain)], Model=Model,Q_r=Q_r,Q_s = Q_s,Q_u=Q_u)
      grad_proposed <- list(grad_r=as.numeric(Allquantities$grad_r), grad_s=as.numeric(Allquantities$grad_s), grad_u=as.numeric(Allquantities$grad_u), cov_r=Allquantities$cov_r, cov_s=Allquantities$cov_s)

      likelihoodproposed<- Allquantities$loglike

      mh.ratio<- exp(likelihoodproposed + priorproposedGs
                     - likelihoodcurrent - priorcurrentGs)

      #print(mh.ratio)

      if(!is.na(mh.ratio) && runif(1) < mh.ratio){
        MC_chain[i, 1:2]<- proposedGs
        MC_chain[i, (ncol(MC_chain)-n_copParams):(ncol(MC_chain)-1)]<- proposedcopPs
        likelihoodcurrent<- likelihoodproposed
        grad_current<- grad_proposed
        JointTPM<- JointTPM1
      }
      else{
        MC_chain[i, 1:2]<- MC_chain[i-1,1:2]
        MC_chain[i, (ncol(MC_chain)-n_copParams):(ncol(MC_chain)-1)]<- MC_chain[i-1, (ncol(MC_chain)-n_copParams):(ncol(MC_chain)-1)]
        }
      }




#      proposedcopPs<- rnorm(n_copParams, mean=MC_chain[i-1, (ncol(MC_chain)-n_copParams):(ncol(MC_chain)-1)], sd=rep(sdCops, n_copParams))
#      if(any(proposedcopPs< -1) || any(proposedcopPs> 1)){
#        MC_chain[i, (ncol(MC_chain)-n_copParams):(ncol(MC_chain)-1)]<- MC_chain[i-1, (ncol(MC_chain)-n_copParams):(ncol(MC_chain)-1)]
#      }else{

#      if(nstrain>3){
#        JointTPM<- ParallelJointTransitionMatrix_copula_cpp2(G(MC_chain[i, 1],MC_chain[i, 2]), K=nstrain, proposedcopPs)
#      }else{
#        JointTPM<- JointTransitionMatrix_copula_cpp(G(MC_chain[i, 1],MC_chain[i, 2]), K=nstrain, proposedcopPs)
#      }
#      JointTPM<- ifelse(JointTPM<=0,1e-6,JointTPM)
#      JointTPM<- ifelse(JointTPM>=1,1-1e-6,JointTPM)

#      Allquantities<- FFBSgradmultstrainLoglikelihood2_cpp(y=y, e_it=e_it, nstrain=nstrain,  r=MC_chain[i, 5+(1:time)], s=MC_chain[i, 5+time+(1:12)], u=MC_chain[i, 5+time+12+(1:ndept)], Gamma=JointTPM, B=MC_chain[i, 5+time+12+ndept+(1:nstrain)], Bits=Bits, a_k=MC_chain[i-1, 5+time+12+ndept+nstrain+(1:nstrain)], Model=Model,Q_r=Q_r,Q_s = Q_s,Q_u=Q_u)
#      grad_proposed <- list(grad_r=as.numeric(Allquantities$grad_r), grad_s=as.numeric(Allquantities$grad_s), grad_u=as.numeric(Allquantities$grad_u), cov_r=Allquantities$cov_r, cov_s=Allquantities$cov_s)

#      likelihoodproposed<- Allquantities$loglike

#      mh.ratio<- exp(likelihoodproposed - likelihoodcurrent) #+ proposalcurrentcop - proposalproposedcop)
      #+ priorproposedcopPs - priorcurrentcopPs)

      #print(mh.ratio)

#      if(!is.na(mh.ratio) && runif(1) < mh.ratio){
#        MC_chain[i, (ncol(MC_chain)-n_copParams):(ncol(MC_chain)-1)]<- proposedcopPs
#        likelihoodcurrent<- likelihoodproposed
#        grad_current<- grad_proposed
#      }
#      else{
#        MC_chain[i, (ncol(MC_chain)-n_copParams):(ncol(MC_chain)-1)]<- MC_chain[i-1, (ncol(MC_chain)-n_copParams):(ncol(MC_chain)-1)]
#        }
#      }
    }
    MC_chain[i, ncol(MC_chain)]<- stationarydist(G(MC_chain[i, 1], MC_chain[i, 2]))[2]

    #Gibbs A_k's update
    if(all(is.finite(as.numeric(Allquantities$poisMean4GibbsUpdate)))){
    MC_chain[i, 5+time+12+ndept+nstrain+(1:nstrain)]<- log(rgamma(nstrain, shape = 0.01+SumYk_vec, rate = as.numeric(Allquantities$poisMean4GibbsUpdate) + 0.01/exp(-15)))
    }else{
      MC_chain[i, 5+time+12+ndept+nstrain+(1:nstrain)]<- MC_chain[i-1, 5+time+12+ndept+nstrain+(1:nstrain)]
    }

    if(i %% 1000 == 0) cat("Iteration:", i, "\n")
  }
  colnames(MC_chain) <- paste(c("G12", "G21", "kappa_r", "kappa_s", "kappa_u", paste("r", 1:time, sep=""), paste("s", 1:12, sep=""), paste("u", 1:ndept, sep=""), paste("B", 1:nstrain, sep=""), paste("a_k", 1:nstrain, sep=""), paste("copulaParam", 1:n_copParams, sep =""), "StationaryDistribution"))
  MC_chain<- as.data.frame(MC_chain)
  end_time <- Sys.time()
  time_taken<- end_time - start_time
  print(time_taken)
  return(MC_chain)
}
