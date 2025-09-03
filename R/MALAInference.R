log_posterior_grad <- function(y, e_it, r, s, u, Q_r, Q_s, Q_u) {
  ndept <- nrow(y)
  time  <- ncol(y)
  month_indexes <- ((1:time - 1) %% 12) + 1

  r_mat <- matrix(rep(r, each = ndept), nrow = ndept)
  s_mat <- matrix(rep(s[month_indexes], each = ndept), nrow = ndept)
  u_mat <- matrix(rep(u, times = time), nrow = ndept)

  log_risk <- r_mat + s_mat + u_mat

  poisMean<- e_it * exp(log_risk)
  delta <- y - poisMean

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

  return(list(grad_r = grad_r, grad_s = grad_s, grad_u = grad_u, cov_r=cov_r, cov_s=cov_s))
}

MALAInference<- function(y, e_it, Model, adjmat, step_sizes, num_iteration = 15000) {
  start_time <- Sys.time()
  ndept <- nrow(y)
  time <- ncol(y)

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

  design_matrix_func <- get(paste0("DesignMatrixModel", Model))
  z_it <- design_matrix_func(original.y, adjmat)[[1]]
  z_it2 <- design_matrix_func(original.y, adjmat)[[2]]

  crudeResults<- crudeEst(y, e_it)
  crudeR<- crudeResults[[1]]
  crudeS<- crudeResults[[2]]
  crudeU<- crudeResults[[3]]
  crudeU<- ifelse(is.nan(crudeU), mean(crudeU[is.finite(crudeU)]), crudeU)
  crudeblock<- floor(time/12)
  crudeblock<- ((crudeblock*12)-11):(crudeblock*12)

  MC_chain<- matrix(NA, nrow=num_iteration, ncol=5+time+12+ndept+2+1)
  initG12<- runif(1)
  initG21<- runif(1)
  initstateD<- state_dist_cpp(initG12, initG21)[2]
  MC_chain[1,]<- c(initG12, initG21, 1/var(crudeR), 1/var(crudeS), 1/var(crudeU), crudeR, crudeS[crudeblock], crudeU, rep(0, 2), initstateD)

  zigmaR<- diag(rep(0.1, time), nrow = time, ncol = time)
  optconstantR<- 2.38^2/(time-2)
  lambdaR<- 1

  likelihoodcurrent<- GeneralLoglikelihood_cpp2(y, MC_chain[1, 5+(1:time)], MC_chain[1, 5+time+(1:12)], MC_chain[1, 5+time+12+(1:ndept)], G(MC_chain[1, 1], MC_chain[1, 2]), e_it, MC_chain[1, 5+time+12+ndept+(1:2)], Model, z_it, z_it2)
  priorcurrentRcomps<- randomwalk2(MC_chain[1, 5+(1:time)], MC_chain[1, 3])
  priorcurrentScomps<- seasonalComp2(MC_chain[1, 5+time+(1:12)], MC_chain[1, 4], RW1PrecMat)
  priorcurrentUcomps<- logIGMRF1(MC_chain[1, 5+time+12+(1:ndept)], MC_chain[1, 5], R, rankdef)

  Q_r<- MC_chain[1,3] * RW2PrecMat
  Q_s<- MC_chain[1,4] * RW1PrecMat
  Q_u<- MC_chain[1,5] * R

  #Compute gradients
  grad_current <- log_posterior_grad(y=y, e_it=e_it, r=MC_chain[1, 5+(1:time)], s=MC_chain[1, 5+time+(1:12)], u=MC_chain[1, 5+time+12+(1:ndept)],
                                     Q_r=Q_r, Q_s=Q_s, Q_u=Q_u)


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

    grad_proposed <- log_posterior_grad(y=y, e_it=e_it, r=MC_chain[i-1, 5+(1:time)],
                                        s=proposedScomps, u=current_u, Q_r=Q_r, Q_s=Q_s, Q_u=Q_u)

    Mmatps<- as.numeric(grad_proposed$cov_s %*% grad_proposed$grad_s)

    likelihoodproposed<- GeneralLoglikelihood_cpp2(y, MC_chain[i-1, 5+(1:time)], proposedScomps, MC_chain[i-1, 5+time+12+(1:ndept)], G(MC_chain[i-1, 1], MC_chain[i-1, 2]), e_it, MC_chain[i-1, 5+time+12+ndept+(1:2)], Model, z_it, z_it2)
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

    grad_proposed <- log_posterior_grad(y=y, e_it=e_it, r=proposedRcomps,
                                        s=MC_chain[i, 5+time+(1:12)], u=current_u, Q_r=Q_r, Q_s=Q_s, Q_u=Q_u)

    Mmatrp<- as.numeric(grad_proposed$cov_r %*% grad_proposed$grad_r)

    q_prop <- mvnfast::dmvn(proposedRcomps, mu = current_r + 0.5 * step_sizes$r^2 * Mmatrc, sigma = grad_current$cov_r * step_sizes$r^2, log = TRUE)
    q_curr <- mvnfast::dmvn(current_r, mu = proposedRcomps + 0.5 * step_sizes$r^2 * Mmatrp, sigma = grad_proposed$cov_r * step_sizes$r^2, log = TRUE)

    likelihoodproposed<- GeneralLoglikelihood_cpp2(y, proposedRcomps, MC_chain[i, 5+time+(1:12)], MC_chain[i-1, 5+time+12+(1:ndept)], G(MC_chain[i-1, 1], MC_chain[i-1, 2]), e_it, MC_chain[i-1, 5+time+12+ndept+(1:2)], Model, z_it, z_it2)
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

    grad_proposed <- log_posterior_grad(y=y, e_it=e_it, r=MC_chain[i, 5+(1:time)],
                                        s=MC_chain[i, 5+time+(1:12)], u=proposedUcomps, Q_r=Q_r, Q_s=Q_s, Q_u=Q_u)


    likelihoodproposed<- GeneralLoglikelihood_cpp2(y, MC_chain[i, 5+(1:time)], MC_chain[i, 5+time+(1:12)], proposedUcomps, G(MC_chain[i-1, 1], MC_chain[i-1, 2]), e_it, MC_chain[i-1, 5+time+12+ndept+(1:2)], Model, z_it, z_it2)

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
    proposedB <- c(0, 0)
    priorcurrentB<- -99999
    priorproposedB<- -99999
  }else if(Model %in% c(1,2,4,5,7)) {
    proposedB <- abs(rnorm(1, mean = MC_chain[i-1, 5+time+12+ndept+1], sd = 0.1))
    proposedB <- c(proposedB, 0)
    priorcurrentB<- dgamma(MC_chain[i-1, 5+time+12+ndept+1], shape = 2, rate = 2, log=TRUE)
    priorproposedB<- dgamma(proposedB[1], shape = 2, rate = 2, log=TRUE)
  }else if(Model %in% c(3,6)){
    proposedB <- abs(rnorm(2, mean = MC_chain[i-1, 5+time+12+ndept+(1:2)], sd = c(0.1, 0.1)))
    priorcurrentB<- sum(dgamma(MC_chain[i-1, 5+time+12+ndept+(1:2)], shape = c(2, 2), rate = c(2,2), log=TRUE))
    priorproposedB<- sum(dgamma(proposedB, shape = c(2, 2), rate = c(2, 2), log=TRUE))
  }

  likelihoodproposed<- GeneralLoglikelihood_cpp2(y, MC_chain[i, 5+(1:time)], MC_chain[i, 5+time+(1:12)], MC_chain[i, 5+time+12+(1:ndept)], G(MC_chain[i-1, 1], MC_chain[i-1,2]), e_it, proposedB, Model,z_it, z_it2)

  mh.ratio<- exp(likelihoodproposed + priorproposedB
                 - likelihoodcurrent - priorcurrentB)

  if(!is.na(mh.ratio) && runif(1) < mh.ratio){
    MC_chain[i, 5+time+12+ndept+(1:2)]<- proposedB
    likelihoodcurrent<- likelihoodproposed
  }
  else{
    MC_chain[i, 5+time+12+ndept+(1:2)]<- MC_chain[i-1, 5+time+12+ndept+(1:2)]
  }

  proposedGs<- abs(rnorm(2,mean=c(MC_chain[i-1,1], MC_chain[i-1, 2]), sd=c(0.1, 0.1)))
  if(proposedGs[1]>1) proposedGs[1]=2-proposedGs[1]
  if(proposedGs[2]>1) proposedGs[2]=2-proposedGs[2]

  priorcurrentGs<- sum(dbeta(MC_chain[i-1,1:2], shape1 = c(2,2), shape2 = c(2,2), log=TRUE))
  priorproposedGs<- sum(dbeta(proposedGs, shape1 = c(2,2), shape2 = c(2,2), log=TRUE))

  likelihoodproposed<- GeneralLoglikelihood_cpp2(y, MC_chain[i, 5+(1:time)], MC_chain[i, 5+time+(1:12)], MC_chain[i, 5+time+12+(1:ndept)], G(proposedGs[1], proposedGs[2]),e_it, MC_chain[i, 5+time+12+ndept+(1:2)], Model,z_it, z_it2)

  mh.ratio<- exp(likelihoodproposed + priorproposedGs
                 - likelihoodcurrent - priorcurrentGs)

  if(!is.na(mh.ratio) && runif(1) < mh.ratio){
    MC_chain[i, 1:2]<- proposedGs
    likelihoodcurrent<- likelihoodproposed
  }
  else{
    MC_chain[i, 1:2]<- MC_chain[i-1,1:2]
  }
  MC_chain[i, 5+time+12+ndept+2+1]<- state_dist_cpp(MC_chain[i, 1], MC_chain[i, 2])[2]
  if(i %% 1000 == 0) cat("Iteration:", i, "\n")
}
  colnames(MC_chain) <- paste(c("G12", "G21", "kappa_r", "kappa_s", "kappa_u", paste("r", 1:time, sep=""), paste("s", 1:12, sep=""), paste("u", 1:ndept, sep=""), paste("B", 1:2, sep=""), "StationaryDistribution"))
  MC_chain<- as.data.frame(MC_chain)
  end_time <- Sys.time()
  time_taken<- end_time - start_time
  print(time_taken)
  return(MC_chain)
}

#malaRes<- MALAInference(y=mixedpop0[[1]], e_it = mixedpop0[[2]], Model = 0, adjmat = sim_adjmat, step_sizes = list("r"=0.3,"s"=0.3,"u"=0.025), num_iteration = 5000)
#mcmc.plot(malaRes)

multlog_posterior_grad <- function(y, e_it, r, s, u, Q_r, Q_s, Q_u, a_k) {
  ndept <- nrow(e_it)
  time  <- ncol(e_it)
  nstrain<- dim(y)[3]
  month_indexes <- ((1:time - 1) %% 12) + 1

  r_mat <- matrix(rep(r, each = ndept), nrow = ndept)
  s_mat <- matrix(rep(s[month_indexes], each = ndept), nrow = ndept)
  u_mat <- matrix(rep(u, times = time), nrow = ndept)

  log_risk <- r_mat + s_mat + u_mat

  poisMean<- matrix(0, nrow = ndept, ncol = time)
  delta<- matrix(0, nrow = ndept, ncol = time)
  for(k in 1:nstrain){
    delta<- delta + (y[,,k] - (e_it * exp(log_risk + a_k[k])))
    poisMean<- poisMean + e_it * exp(log_risk + a_k[k])
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

  return(list(grad_r = grad_r, grad_s = grad_s, grad_u = grad_u, cov_r=cov_r, cov_s=cov_s))
}

multMALAInference<- function(y, e_it, Model, adjmat, step_sizes, num_iteration = 15000, independentChains=0) {
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

  crudeResults<- crudeEst(sumY, e_it)
  crudeR<- crudeResults[[1]] - mean(crudeResults[[1]])
  crudeS<- crudeResults[[2]]
  crudeU<- crudeResults[[3]]
  crudeU<- ifelse(is.nan(crudeU), mean(crudeU[is.finite(crudeU)]), crudeU)
  crudeblock<- floor(time/12)
  crudeblock<- ((crudeblock*12)-11):(crudeblock*12)

  MC_chain<- matrix(NA, nrow=num_iteration, ncol=5+time+12+ndept+nstrain+nstrain+1)
  initG12<- runif(1)
  initG21<- runif(1)
  initstateD<- state_dist_cpp(initG12, initG21)[2]
  MC_chain[1,]<- c(initG12, initG21, 1/var(crudeR), 1/var(crudeS), 1/var(crudeU), crudeR, crudeS[crudeblock-12], crudeU, rep(0, nstrain), rep(mean(crudeResults[[1]]), nstrain), initstateD)

  #flag missing data for inference
  y<- ifelse(is.na(y), -1, y)
  yflat<- as.numeric(aperm(y, c(2,1,3)))

  likelihoodcurrent<- multGeneralLoglikelihood_cpp2(y=yflat,ndept=ndept,time=time,nstrain=nstrain,a_k=MC_chain[1, 5+time+12+ndept+nstrain+(1:nstrain)], r=MC_chain[1,5+(1:time)],s=MC_chain[1,5+time+(1:12)],u=MC_chain[1,5+time+12+(1:ndept)],Gamma=G(MC_chain[1,1],MC_chain[1,2]),e_it=e_it, B=MC_chain[1, 5+time+12+ndept+(1:nstrain)], Bits=Bits, model=Model,independentChains=independentChains)
  priorcurrentRcomps<- randomwalk2(MC_chain[1, 5+(1:time)], MC_chain[1, 3])
  priorcurrentScomps<- seasonalComp2(MC_chain[1, 5+time+(1:12)], MC_chain[1, 4], RW1PrecMat)
  priorcurrentUcomps<- logIGMRF1(MC_chain[1, 5+time+12+(1:ndept)], MC_chain[1, 5], R, rankdef)

  Q_r<- MC_chain[1,3] * RW2PrecMat
  Q_s<- MC_chain[1,4] * RW1PrecMat
  Q_u<- MC_chain[1,5] * R

  #Compute gradients
  grad_current <- multlog_posterior_grad(y=y, e_it=e_it, r=MC_chain[1, 5+(1:time)], s=MC_chain[1, 5+time+(1:12)], u=MC_chain[1, 5+time+12+(1:ndept)],
                                     Q_r=Q_r, Q_s=Q_s, Q_u=Q_u, a_k=MC_chain[1, 5+time+12+ndept+nstrain+(1:nstrain)])


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

    grad_proposed <- multlog_posterior_grad(y=y, e_it=e_it, r=MC_chain[i-1, 5+(1:time)],
                                        s=proposedScomps, u=current_u, Q_r=Q_r, Q_s=Q_s, Q_u=Q_u, a_k=MC_chain[i-1, 5+time+12+ndept+nstrain+(1:nstrain)])

    Mmatps<- as.numeric(grad_proposed$cov_s %*% grad_proposed$grad_s)

    likelihoodproposed<- multGeneralLoglikelihood_cpp2(y=yflat,ndept=ndept,time=time,nstrain=nstrain,a_k=MC_chain[i-1, 5+time+12+ndept+nstrain+(1:nstrain)],r=MC_chain[i-1,5+(1:time)],s=proposedScomps,u=MC_chain[i-1,5+time+12+(1:ndept)],Gamma=G(MC_chain[i-1,1], MC_chain[i-1,2]),e_it=e_it, B=MC_chain[i-1, 5+time+12+ndept+(1:nstrain)], Bits=Bits,model=Model,independentChains=independentChains)
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

    grad_proposed <- multlog_posterior_grad(y=y, e_it=e_it, r=proposedRcomps,
                                        s=MC_chain[i, 5+time+(1:12)], u=current_u, Q_r=Q_r, Q_s=Q_s, Q_u=Q_u, a_k=MC_chain[i-1, 5+time+12+ndept+nstrain+(1:nstrain)])

    Mmatrp<- as.numeric(grad_proposed$cov_r %*% grad_proposed$grad_r)

    q_prop <- mvnfast::dmvn(proposedRcomps, mu = current_r + 0.5 * step_sizes$r^2 * Mmatrc, sigma = grad_current$cov_r * step_sizes$r^2, log = TRUE)
    q_curr <- mvnfast::dmvn(current_r, mu = proposedRcomps + 0.5 * step_sizes$r^2 * Mmatrp, sigma = grad_proposed$cov_r * step_sizes$r^2, log = TRUE)

    likelihoodproposed<- multGeneralLoglikelihood_cpp2(y=yflat,ndept=ndept,time=time,nstrain=nstrain,a_k=MC_chain[i-1, 5+time+12+ndept+nstrain+(1:nstrain)], r=proposedRcomps, s=MC_chain[i, 5+time+(1:12)], u=MC_chain[i-1, 5+time+12+(1:ndept)], Gamma=G(MC_chain[i-1, 1], MC_chain[i-1, 2]), e_it=e_it, B=MC_chain[i-1, 5+time+12+ndept+(1:nstrain)], Bits=Bits, model=Model, independentChains = independentChains)
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

    grad_proposed <- multlog_posterior_grad(y=y, e_it=e_it, r=MC_chain[i, 5+(1:time)],
                                        s=MC_chain[i, 5+time+(1:12)], u=proposedUcomps, Q_r=Q_r, Q_s=Q_s, Q_u=Q_u, a_k=MC_chain[i-1, 5+time+12+ndept+nstrain+(1:nstrain)])


    likelihoodproposed<- multGeneralLoglikelihood_cpp2(y=yflat,ndept=ndept,time=time,nstrain=nstrain,a_k=MC_chain[i-1, 5+time+12+ndept+nstrain+(1:nstrain)],r=MC_chain[i,5+(1:time)],s=MC_chain[i,5+time+(1:12)],u=proposedUcomps,Gamma=G(MC_chain[i-1,1],MC_chain[i-1,2]),e_it=e_it, B=MC_chain[i-1, 5+time+12+ndept+(1:nstrain)], Bits=Bits, model=Model,independentChains=independentChains)

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
      proposedB <- rep(0, nstrain)
      priorcurrentB<- -99999
      priorproposedB<- -99999
    }else{
      proposedB <- abs(rnorm(nstrain, mean = MC_chain[i-1, 5+time+12+ndept+(1:nstrain)], sd = rep(0.1, nstrain)))
      priorcurrentB<- sum(dgamma(MC_chain[i-1, 5+time+12+ndept+(1:nstrain)], shape = rep(2, nstrain), rate = rep(2,nstrain), log=TRUE))
      priorproposedB<- sum(dgamma(proposedB, shape = rep(2, nstrain), rate = rep(2, nstrain), log=TRUE))
    }

    likelihoodcurrent<- multGeneralLoglikelihood_cpp2(y=yflat,ndept=ndept,time=time,nstrain=nstrain,a_k=MC_chain[i-1, 5+time+12+ndept+nstrain+(1:nstrain)], r=MC_chain[i, 5+(1:time)], s=MC_chain[i, 5+time+(1:12)], u=MC_chain[i, 5+time+12+(1:ndept)], Gamma=G(MC_chain[i-1, 1],MC_chain[i-1,2]), e_it=e_it, B=MC_chain[i-1, 5+time+12+ndept+(1:nstrain)], Bits=Bits, model=Model, independentChains=independentChains)
    likelihoodproposed<- multGeneralLoglikelihood_cpp2(y=yflat,ndept=ndept,time=time,nstrain=nstrain,a_k=MC_chain[i-1, 5+time+12+ndept+nstrain+(1:nstrain)], r=MC_chain[i, 5+(1:time)], s=MC_chain[i, 5+time+(1:12)], u=MC_chain[i, 5+time+12+(1:ndept)], Gamma=G(MC_chain[i-1, 1], MC_chain[i-1,2]), e_it=e_it, B=proposedB, Bits=Bits, model=Model,independentChains=independentChains)

    mh.ratio<- exp(likelihoodproposed + priorproposedB
                   - likelihoodcurrent - priorcurrentB)

    #print(paste("mh.ratioB = ", mh.ratio))

    if(!is.na(mh.ratio) && runif(1) < mh.ratio){
      MC_chain[i, 5+time+12+ndept+(1:nstrain)]<- proposedB
    }
    else{
      MC_chain[i, 5+time+12+ndept+(1:nstrain)]<- MC_chain[i-1, 5+time+12+ndept+(1:nstrain)]
    }

    proposedGs<- abs(rnorm(2,mean=c(MC_chain[i-1,1], MC_chain[i-1, 2]), sd=c(0.1, 0.1)))
    if(proposedGs[1]>1) proposedGs[1]=2-proposedGs[1]
    if(proposedGs[2]>1) proposedGs[2]=2-proposedGs[2]

    priorcurrentGs<- sum(dbeta(MC_chain[i-1,1:2], shape1 = c(2,2), shape2 = c(2,2), log=TRUE))
    priorproposedGs<- sum(dbeta(proposedGs, shape1 = c(2,2), shape2 = c(2,2), log=TRUE))

    likelihoodcurrent<- multGeneralLoglikelihood_cpp2(y=yflat,ndept=ndept,time=time,nstrain=nstrain,a_k=MC_chain[i-1, 5+time+12+ndept+nstrain+(1:nstrain)], r=MC_chain[i,5+(1:time)], s=MC_chain[i,5+time+(1:12)], u=MC_chain[i, 5+time+12+(1:ndept)], Gamma=G(MC_chain[i-1,1],MC_chain[i-1,2]),e_it=e_it, B=MC_chain[i, 5+time+12+ndept+(1:nstrain)], Bits=Bits, model=Model, independentChains=independentChains)
    likelihoodproposed<- multGeneralLoglikelihood_cpp2(y=yflat,ndept=ndept,time=time,nstrain=nstrain,a_k=MC_chain[i-1, 5+time+12+ndept+nstrain+(1:nstrain)], r=MC_chain[i, 5+(1:time)], s=MC_chain[i, 5+time+(1:12)], u=MC_chain[i, 5+time+12+(1:ndept)], Gamma=G(proposedGs[1], proposedGs[2]),e_it=e_it, B=MC_chain[i, 5+time+12+ndept+(1:nstrain)], Bits=Bits, model=Model,independentChains=independentChains)

    mh.ratio<- exp(likelihoodproposed + priorproposedGs
                   - likelihoodcurrent - priorcurrentGs)

    #print(mh.ratio)

    if(!is.na(mh.ratio) && runif(1) < mh.ratio){
      MC_chain[i, 1:2]<- proposedGs
    }
    else{
      MC_chain[i, 1:2]<- MC_chain[i-1,1:2]
    }
    MC_chain[i, 5+time+12+ndept+nstrain+nstrain+1]<- state_dist_cpp(MC_chain[i, 1], MC_chain[i, 2])[2]

    #Random-wal Ak's update
    proposeda_k <- rnorm(nstrain, mean = MC_chain[i-1, 5+time+12+ndept+nstrain+(1:nstrain)], sd = rep(0.03, nstrain))

    likelihoodcurrent<- multGeneralLoglikelihood_cpp2(y=yflat,ndept=ndept,time=time,nstrain=nstrain,a_k=MC_chain[i-1, 5+time+12+ndept+nstrain+(1:nstrain)], r=MC_chain[i, 5+(1:time)], s=MC_chain[i, 5+time+(1:12)], u=MC_chain[i, 5+time+12+(1:ndept)], Gamma=G(MC_chain[i, 1],MC_chain[i,2]), e_it=e_it, B=MC_chain[i, 5+time+12+ndept+(1:nstrain)], Bits=Bits, model=Model, independentChains=independentChains)
    likelihoodproposed<- multGeneralLoglikelihood_cpp2(y=yflat,ndept=ndept,time=time,nstrain=nstrain,a_k=proposeda_k, r=MC_chain[i, 5+(1:time)], s=MC_chain[i, 5+time+(1:12)], u=MC_chain[i, 5+time+12+(1:ndept)], Gamma=G(MC_chain[i, 1], MC_chain[i,2]), e_it=e_it, B=MC_chain[i, 5+time+12+ndept+(1:nstrain)], Bits=Bits, model=Model,independentChains=independentChains)

    mh.ratio<- exp(likelihoodproposed - likelihoodcurrent)

    if(!is.na(mh.ratio) && runif(1) < mh.ratio){
      MC_chain[i, 5+time+12+ndept+nstrain+(1:nstrain)]<- proposeda_k
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

#multmalaRes<- multMALAInference(y=multmod0nstrain5[[1]], e_it = multmod0nstrain5[[2]], Model = 0, adjmat = sim_adjmat, step_sizes = list("r"=0.3,"s"=0.3,"u"=0.025), num_iteration = 5000)
