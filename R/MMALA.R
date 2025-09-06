#Gradients, inverse-Fisher's, loglikelihood
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

    return(list(loglike = loglike, grad_r = grad_r, grad_s = grad_s, grad_u = grad_u, cov_r=cov_r, cov_s=cov_s))
  }else{
    # outputs to build
    grad_r <- numeric(time)
    grad_s <- numeric(12)
    grad_u <- numeric(ndept)

    Fisher_r <- matrix(0, nrow = time, ncol = time)
    Fisher_s <- matrix(0, nrow = 12,   ncol = 12)
    Fisher_u <- matrix(0, nrow = ndept, ncol = ndept)

    loglike_total <- 0

    JointTPM <- JointTransitionMatrix(gamma = Gamma, K = nstrain)
    logJointTPM <- log(JointTPM)

    E_lambda_tk <- array(0, dim = c(ndept,time,nstrain))   #Expected Poisson mean

    for(i in 1:ndept){

      logEmissions <- matrix(NA, nrow = time, ncol = nstate)
      lambda_array  <- array(0, dim = c(time, nstate, nstrain))

      for(t in 1:time){
        month_index <- (t-1) %% 12 + 1
        for(n in 1:nstate){
          for(k in 1:nstrain){
            lambda_array[t,n,k] <- e_it[i,t] * exp(a_k[k] + r[t] + s[month_index] + u[i] + as.numeric(B %*% Bits[n, ]))
          }
          logEmissions[t,n] <- sum(dpois(y[i,t,], lambda = lambda_array[t,n,], log = TRUE))
        }
      }

      # forward pass
      loginit <- log(stationarydist(JointTPM))
      logalpha <- matrix(-Inf, nrow = time, ncol = nstate)
      logalpha[1, ] <- loginit + logEmissions[1, ]
      for(t in 2:time){
        for(n in 1:nstate){
          logalpha[t,n] <- logSumExp_cpp(logalpha[t-1, ] + logJointTPM[, n]) + logEmissions[t,n]
        }
      }

      loglik_i <- logSumExp_cpp(logalpha[time, ])
      loglike_total <- loglike_total + loglik_i

      # backward pass
      logbeta <- matrix(-Inf, nrow = time, ncol = nstate)
      logbeta[time, ] <- 0
      for(t in seq(time-1, 1, by = -1)){
        for(m in 1:nstate){
          logbeta[t,m] <- logSumExp_cpp( logJointTPM[m, ] + logEmissions[t+1, ] + logbeta[t+1, ] )
        }
      }

      #Marginal posterior probabilities P_s
      logP_s <- matrix(NA, nrow = time, ncol = nstate)
      for(t in 1:time){
        logP_s[t, ] <- logalpha[t, ] + logbeta[t, ] - loglik_i
      }
      P_s<- exp(logP_s)

      #Expected Poisson mean
      for(t in 1:time){
        for(k in 1:nstrain){
          E_lambda_tk[i,t,k] <- sum(P_s[t, ] * lambda_array[t, , k])
        }
      }
    }

    poisMean<- matrix(0, nrow = ndept, ncol = time)
    delta<- matrix(0, nrow = ndept, ncol = time)
    for(k in 1:nstrain){
      delta<- delta + (y[,,k] - E_lambda_tk[,,k])
      poisMean<- poisMean + E_lambda_tk[,,k]
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

    return(list(loglike = loglike_total, grad_r = grad_r, grad_s = grad_s, grad_u = grad_u, cov_r=cov_r, cov_s=cov_s))
  }
}

#Riemann Manifold Langevin updates
multMMALAInference<- function(y, e_it, Model, adjmat, step_sizes, num_iteration = 15000, independentChains=0) {
  start_time <- Sys.time()
  ndept <- nrow(e_it)
  time <- ncol(e_it)
  nstrain<- dim(y)[3]
  Bits<- encodeBits(nstrain)

  R<- -1 * adjmat
  diag(R)<- -rowSums(R, na.rm = T)

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

  sumY<- y[,,1]
  for(k in 2:nstrain){
    sumY<- sumY + y[,,k]
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
      MC_chain[i, 5+time+(1:12)]<- proposedScomps - mean(proposedScomps)
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

    log_alpha_u <- likelihoodproposed + priorproposedUcomps + q_curr - likelihoodcurrent - priorcurrentUcomps - q_prop
    if (is.finite(log_alpha_u) && log(runif(1)) < log_alpha_u){
      MC_chain[i, 5+time+12+(1:ndept)]<- proposedUcomps
      likelihoodcurrent<- likelihoodproposed
      priorcurrentUcomps<- priorproposedUcomps
      grad_current<- grad_proposed
    }else{
      MC_chain[i, 5+time+12+(1:ndept)]<- MC_chain[i-1, 5+time+12+(1:ndept)]
    }

    if(Model == 0){
      MC_chain[i, 5+time+12+ndept+(1:nstrain)]<-  MC_chain[i-1, 5+time+12+ndept+(1:nstrain)]
      MC_chain[i, 1:2]<- MC_chain[i-1, 1:2]
    }else{
      proposedB <- abs(rnorm(nstrain, mean = MC_chain[i-1, 5+time+12+ndept+(1:nstrain)], sd = rep(0.03, nstrain)))
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

    proposedGs<- abs(rnorm(2,mean=c(MC_chain[i-1,1], MC_chain[i-1, 2]), sd=c(0.1, 0.1)))
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

    #Random-wal Ak's update
    proposeda_k <- rnorm(nstrain, mean = MC_chain[i-1, 5+time+12+ndept+nstrain+(1:nstrain)], sd = rep(0.03, nstrain))

    Allquantities<- gradmultstrainLoglikelihood2(y=y, e_it=e_it, nstrain=nstrain,  r=MC_chain[i, 5+(1:time)], s=MC_chain[i, 5+time+(1:12)], u=MC_chain[i, 5+time+12+(1:ndept)], Gamma=G(MC_chain[i, 1],MC_chain[i, 2]), B=MC_chain[i, 5+time+12+ndept+(1:nstrain)], Bits=Bits, a_k=proposeda_k, Model=Model,Q_r=Q_r,Q_s = Q_s,Q_u=Q_u)
    grad_proposed <- list(grad_r=Allquantities$grad_r, grad_s=Allquantities$grad_s, grad_u=Allquantities$grad_u, cov_r=Allquantities$cov_r, cov_s=Allquantities$cov_s)

    likelihoodproposed<- Allquantities$loglike

    mh.ratio<- exp(likelihoodproposed - likelihoodcurrent)

    if(!is.na(mh.ratio) && runif(1) < mh.ratio){
      MC_chain[i, 5+time+12+ndept+nstrain+(1:nstrain)]<- proposeda_k
      likelihoodcurrent<- likelihoodproposed
      grad_current<- grad_proposed
    }
    else{
      MC_chain[i, 5+time+12+ndept+nstrain+(1:nstrain)]<- MC_chain[i-1, 5+time+12+ndept+nstrain+(1:nstrain)]
    }
    if(i %% 1000 == 0) cat("Iteration:", i, "\n")
  }
  colnames(MC_chain) <- paste(c("G12", "G21", "kappa_r", "kappa_s", "kappa_u", paste("r", 1:time, sep=""), paste("s", 1:12, sep=""), paste("u", 1:ndept, sep=""), paste("B", 1:nstrain, sep=""), paste("a_k", 1:nstrain, sep=""), "StationaryDistribution"))
  MC_chain<- as.data.frame(MC_chain)
  end_time <- Sys.time()
  time_taken<- end_time - start_time
  print(time_taken)
  return(MC_chain)
}


#multMmalaRes2<- multMMALAInference(y=multmod1nstrain5[[1]], e_it = multmod1nstrain5[[2]], Model = 1, adjmat = sim_adjmat, step_sizes = list("r"=0.3,"s"=0.3,"u"=0.025), num_iteration = 20000)
#set.seed(111);newmultmod1nstrain5<- Multstrain.simulate(Model = 1, time=60, adj.matrix = sim_adjmat, nstrain=5)
#set.seed(111);newmultmod1nstrain3<- Multstrain.simulate(Model = 1, time=60, adj.matrix = sim_adjmat, nstrain=3, B=c(1.2,0.95,1.1))
