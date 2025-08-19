log_posterior_grad <- function(y, e_it, r, s, u, Q_r, Q_u, kappa_r, kappa_s, kappa_u) {
  ndept <- nrow(y)
  time <- ncol(y)
  # Compute eta = log(lambda_it)
  eta <- matrix(0, ndept, time)
  for (i in 1:ndept) {
    for (t in 1:time) {
      month_index <- (t - 1) %% 12 + 1
      eta[i, t] <- r[t] + s[month_index] + u[i] #+ x[i, t] * beta
    }
  }

  # Compute delta = y - e * exp(eta)
  lambda <- exp(eta)
  delta <- y - e_it * lambda

  # Temporal trend r
  grad_r <- colSums(delta) - kappa_r * as.numeric(Q_r %*% r)

  # Seasonal s
  grad_s <- numeric(12)
  for (month_index in 1:12) {
    t_idx <- which(((1:time - 1) %% 12 + 1) == month_index)
    grad_s[month_index] <- sum(delta[, t_idx]) - kappa_s * s[j]
  }

  # Spatial u
  grad_u <- rowSums(delta) - kappa_u * as.numeric(Q_u %*% u)

  return(list(grad_r = grad_r, grad_s = grad_s, grad_u = grad_u))
}

MALAInference<- function(y, e_it, Model, adjmat, step_sizes, num_iteration = 1000) {
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
  MC_chain[1,]<- c(initG12, initG21, runif(1, 0, 1000), runif(1, 0, 1000), runif(1, 0, 30), crudeR, crudeS[crudeblock], crudeU, rep(0, 2), initstateD)

  zigmaR<- diag(rep(0.1, time), nrow = time, ncol = time)
  optconstantR<- 2.38^2/(time-2)
  lambdaR<- 1

  likelihoodcurrent<- GeneralLoglikelihood_cpp2(y, MC_chain[1, 5+(1:time)], MC_chain[1, 5+time+(1:12)], MC_chain[1, 5+time+12+(1:ndept)], G(MC_chain[1, 1], MC_chain[1, 2]), e_it, MC_chain[1, 5+time+12+ndept+(1:2)], Model, z_it, z_it2)
  priorcurrentRcomps<- randomwalk2(MC_chain[1, 5+(1:time)], MC_chain[1, 3])
  priorcurrentScomps<- seasonalComp2(MC_chain[1, 5+time+(1:12)], MC_chain[1, 4], RW1PrecMat)
  priorcurrentUcomps<- logIGMRF1(MC_chain[1, 5+time+12+(1:ndept)], MC_chain[1, 5], R, rankdef)

  for (i in 2:num_iteration) {

    MC_chain[i,3]<- rgamma(1, shape = 1 + (time-2)/2, rate = 0.0001 + (t(MC_chain[i-1, 5+(1:time)]) %*% RW2PrecMat %*% MC_chain[i-1, 5+(1:time)])/2)
    MC_chain[i,4]<- rgamma(1, shape = 1 + 11/2, rate = (time / 12.0) * 0.001 + (t(MC_chain[i-1, 5+time+(1:12)]) %*% RW1PrecMat %*% MC_chain[i-1, 5+time+(1:12)])/2)
    MC_chain[i,5]<- rgamma(1, shape = 1 + (ndept-1)/2, rate = 0.01 + (t(MC_chain[i-1, 5+time+12+(1:ndept)]) %*% R %*% MC_chain[i-1, 5+time+12+(1:ndept)])/2)

    Q_r<- MC_chain[i,3] * RW2PrecMat
    Q_s<- MC_chain[i,4] * RW1PrecMat
    Q_u<- MC_chain[i,5] * R

    #Compute gradients
    grad <- log_posterior_grad(y=y, e_it=e_it, r=MC_chain[i-1, 5+(1:time)], s=MC_chain[i-1, 5+time+(1:12)], u=MC_chain[i-1, 5+time+12+(1:ndept)],
                               Q_r=Q_r, Q_u=Q_u,kappa_r=MC_chain[i,3], kappa_s=MC_chain[i,4], kappa_u=MC_chain[i,5])

    #Update r
    eps_r <- rnorm(time)
    proposedRcomps<- MC_chain[i-1, 5+(1:time)] + 0.5 * step_sizes$r^2 * grad$grad_r + step_sizes$r * eps_r
    #likelihoodcurrent<- GeneralLoglikelihood_cpp2(y, MC_chain[i-1, 5+(1:time)], MC_chain[i-1, 5+time+(1:12)], MC_chain[i-1, 5+time+12+(1:ndept)], G(MC_chain[i-1, 1], MC_chain[i-1, 2]), e_it, MC_chain[i-1, 5+time+12+ndept+(1:2)], Model, z_it, z_it2)
    likelihoodproposed<- GeneralLoglikelihood_cpp2(y, proposedRcomps, MC_chain[i-1, 5+time+(1:12)], MC_chain[i-1, 5+time+12+(1:ndept)], G(MC_chain[i-1, 1], MC_chain[i-1, 2]), e_it, MC_chain[i-1, 5+time+12+ndept+(1:2)], Model, z_it, z_it2)
    q_prop <- sum(dnorm(proposedRcomps, mean = MC_chain[i-1, 5+(1:time)] + 0.5 * step_sizes$r^2 * grad$grad_r, sd = step_sizes$r, log = TRUE))
    q_curr <- sum(dnorm(MC_chain[i-1, 5+(1:time)], mean = proposedRcomps + 0.5 * step_sizes$r^2 * grad$grad_r, sd = step_sizes$r, log = TRUE))
    #priorcurrentRcomps<- randomwalk2(MC_chain[i-1, 5+(1:time)], MC_chain[i, 3])
    priorproposedRcomps<- randomwalk2(proposedRcomps, MC_chain[i, 3])

    log_alpha_r <- likelihoodproposed + priorproposedRcomps + q_curr - likelihoodcurrent - priorcurrentRcomps - q_prop
    if (is.finite(log_alpha_r) && log(runif(1)) < log_alpha_r){
      MC_chain[i, 5+(1:time)] <- proposedRcomps
      likelihoodcurrent<- likelihoodproposed
      priorcurrentRcomps<- priorproposedRcomps
      #cat("accepted")
    }else{
      MC_chain[i, 5+(1:time)]<- MC_chain[i-1, 5+(1:time)]
    }

     #Update s
    eps_s <- rnorm(12)
    proposedScomps <- MC_chain[i-1, 5+time+(1:12)] + 0.5 * step_sizes$s^2 * grad$grad_s + step_sizes$s * eps_s
    #likelihoodcurrent<- GeneralLoglikelihood_cpp2(y, MC_chain[i, 5+(1:time)], MC_chain[i-1, 5+time+(1:12)], MC_chain[i-1, 5+time+12+(1:ndept)], G(MC_chain[i-1, 1], MC_chain[i-1, 2]), e_it, MC_chain[i-1, 5+time+12+ndept+(1:2)], Model, z_it, z_it2)
    likelihoodproposed<- GeneralLoglikelihood_cpp2(y, MC_chain[i, 5+(1:time)], proposedScomps, MC_chain[i-1, 5+time+12+(1:ndept)], G(MC_chain[i-1, 1], MC_chain[i-1, 2]), e_it, MC_chain[i-1, 5+time+12+ndept+(1:2)], Model, z_it, z_it2)
    q_prop <- sum(dnorm(proposedScomps, mean = MC_chain[i-1, 5+time+(1:12)] + 0.5 * step_sizes$s^2 * grad$grad_s, sd = step_sizes$s, log = TRUE))
    q_curr <- sum(dnorm(MC_chain[i-1, 5+time+(1:12)], mean = proposedScomps + 0.5 * step_sizes$s^2 * grad$grad_s, sd = step_sizes$s, log = TRUE))
    #priorcurrentScomps<- seasonalComp2(MC_chain[i-1, 5+time+(1:12)], MC_chain[i, 4], RW1PrecMat)
    priorproposedScomps<- seasonalComp2(proposedScomps, MC_chain[i, 4], RW1PrecMat)

    log_alpha_s <- likelihoodproposed + priorproposedScomps + q_curr - likelihoodcurrent - priorcurrentScomps - q_prop
    if (is.finite(log_alpha_s) && log(runif(1)) < log_alpha_s){
      MC_chain[i, 5+time+(1:12)]<- proposedScomps
      likelihoodcurrent<- likelihoodproposed
      priorcurrentScomps<- priorproposedScomps
    }else{
      MC_chain[i, 5+time+(1:12)]<- MC_chain[i-1, 5+time+(1:12)]
    }

    #Update u
    eps_u <- rnorm(ndept)
    proposedUcomps <- MC_chain[i-1, 5+time+12+(1:ndept)] + 0.5 * step_sizes$u^2 * grad$grad_u + step_sizes$u * eps_u
    #likelihoodcurrent<- GeneralLoglikelihood_cpp2(y, MC_chain[i, 5+(1:time)], MC_chain[i, 5+time+(1:12)], MC_chain[i-1, 5+time+12+(1:ndept)], G(MC_chain[i-1, 1], MC_chain[i-1, 2]), e_it, MC_chain[i-1, 5+time+12+ndept+(1:2)], Model, z_it, z_it2)
    likelihoodproposed<- GeneralLoglikelihood_cpp2(y, MC_chain[i, 5+(1:time)], MC_chain[i, 5+time+(1:12)], proposedUcomps, G(MC_chain[i-1, 1], MC_chain[i-1, 2]), e_it, MC_chain[i-1, 5+time+12+ndept+(1:2)], Model, z_it, z_it2)
    q_prop <- sum(dnorm(proposedUcomps, mean = MC_chain[i-1, 5+time+12+(1:ndept)] + 0.5 * step_sizes$u^2 * grad$grad_u, sd = step_sizes$u, log = TRUE))
    q_curr <- sum(dnorm(MC_chain[i-1, 5+time+12+(1:ndept)], mean = proposedUcomps + 0.5 * step_sizes$u^2 * grad$grad_u, sd = step_sizes$u, log = TRUE))
    #priorcurrentUcomps<- logIGMRF1(MC_chain[i-1, 5+time+12+(1:ndept)], MC_chain[i, 5], R, rankdef)
    priorproposedUcomps<- logIGMRF1(proposedUcomps, MC_chain[i, 5], R, rankdef)

    log_alpha_u <- likelihoodproposed + priorproposedUcomps + q_curr - likelihoodcurrent - priorcurrentUcomps - q_prop
    if (is.finite(log_alpha_u) && log(runif(1)) < log_alpha_u){
      MC_chain[i, 5+time+12+(1:ndept)]<- proposedUcomps
      likelihoodcurrent<- likelihoodproposed
      priorcurrentUcomps<- priorproposedUcomps
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

  #likelihoodcurrent<- GeneralLoglikelihood_cpp2(y, MC_chain[i, 5+(1:time)], MC_chain[i, 5+time+(1:12)], MC_chain[i, 5+time+12+(1:ndept)], G(MC_chain[i-1, 1],MC_chain[i-1,2]), e_it, MC_chain[i-1, 5+time+12+ndept+(1:2)], Model, z_it, z_it2)
  likelihoodproposed<- GeneralLoglikelihood_cpp2(y, MC_chain[i, 5+(1:time)], MC_chain[i, 5+time+(1:12)], MC_chain[i, 5+time+12+(1:ndept)], G(MC_chain[i-1, 1], MC_chain[i-1,2]), e_it, proposedB, Model,z_it, z_it2)

  mh.ratio<- exp(likelihoodproposed + priorproposedB
                 - likelihoodcurrent - priorcurrentB)

  #print(paste("mh.ratioB = ", mh.ratio))

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

  #likelihoodcurrent<- GeneralLoglikelihood_cpp2(y, MC_chain[i,5+(1:time)], MC_chain[i,5+time+(1:12)], MC_chain[i, 5+time+12+(1:ndept)], G(MC_chain[i-1,1],MC_chain[i-1,2]),e_it, MC_chain[i, 5+time+12+ndept+(1:2)], Model, z_it, z_it2)
  likelihoodproposed<- GeneralLoglikelihood_cpp2(y, MC_chain[i, 5+(1:time)], MC_chain[i, 5+time+(1:12)], MC_chain[i, 5+time+12+(1:ndept)], G(proposedGs[1], proposedGs[2]),e_it, MC_chain[i, 5+time+12+ndept+(1:2)], Model,z_it, z_it2)

  mh.ratio<- exp(likelihoodproposed + priorproposedGs
                 - likelihoodcurrent - priorcurrentGs)

  #print(mh.ratio)

  if(!is.na(mh.ratio) && runif(1) < mh.ratio){
    MC_chain[i, 1:2]<- proposedGs
    likelihoodcurrent<- likelihoodproposed
  }
  else{
    MC_chain[i, 1:2]<- MC_chain[i-1,1:2]
  }
  MC_chain[i, 5+time+12+ndept+2+1]<- state_dist_cpp(MC_chain[i, 1], MC_chain[i, 2])[2]
  #Adapting zigmaR
  if(i==5){
    epsilonR<- 0.007
    XnR<- MC_chain[1:i, 5+(1:time)]
    XnbarR <- colMeans(XnR)
    zigmaR <- cov(XnR) + epsilonR * diag(rep(1, time))
    zigmaR<- optconstantR * zigmaR
  } else if (i > 5){

    ### Using random walk after 5 conditional prior proposals
    proposedRcomps<- mvnfast::rmvn(1, mu = MC_chain[i, 5+(1:time)], sigma = zigmaR)

    priorcurrentRcomps<- randomwalk2(MC_chain[i, 5+(1:time)], MC_chain[i, 3])
    priorproposedRcomps<- randomwalk2(proposedRcomps, MC_chain[i, 3])

    proposalproposedRcomps<- mvnfast::dmvn(proposedRcomps, mu = MC_chain[i, 5+(1:time)], sigma = zigmaR, log = TRUE)
    proposalcurrentRcomps<- mvnfast::dmvn(MC_chain[i, 5+(1:time)], mu = proposedRcomps, sigma = zigmaR, log = TRUE)

    likelihoodcurrent<- GeneralLoglikelihood_cpp2(y,MC_chain[i,5+(1:time)],MC_chain[i, 5+time+(1:12)],MC_chain[i,5+time+12+(1:ndept)],G(MC_chain[i,1],MC_chain[i,2]),e_it, MC_chain[i, 5+time+12+ndept+(1:2)], Model,z_it, z_it2)
    likelihoodproposed<- GeneralLoglikelihood_cpp2(y,proposedRcomps,MC_chain[i, 5+time+(1:12)],MC_chain[i,5+time+12+(1:ndept)],G(MC_chain[i,1],MC_chain[i,2]),e_it, MC_chain[i, 5+time+12+ndept+(1:2)], Model,z_it, z_it2)

    mh.ratioR<- exp(likelihoodproposed + priorproposedRcomps + proposalcurrentRcomps
                    - likelihoodcurrent - priorcurrentRcomps - proposalproposedRcomps)

    #print(mh.ratioR)

    if(!is.na(mh.ratioR) && runif(1) < mh.ratioR){
      MC_chain[i,5+(1:time)]<- proposedRcomps
    }
    else{
      MC_chain[i,5+(1:time)]<- MC_chain[i,5+(1:time)]
    }

    XnbarPrevR <- XnbarR
    XnbarR <- (i*XnbarR + MC_chain[i, 5+(1:time)])/(i+1)
    zigmaR <- ((i-1)*zigmaR + tcrossprod(MC_chain[i, 5+(1:time)]) + i*tcrossprod(XnbarPrevR) - (i+1)*tcrossprod(XnbarR) + epsilonR*diag(rep(1,time)))/i
    #Robbins Munro tuning
    lambdaR<- lambdaR * exp((2/max(1, i-5)) * (min(mh.ratioR, 1) - 0.234))
    zigmaR<- lambdaR* optconstantR * zigmaR
    #print(zigmaR)
  }
  if(i %% 1000 == 0) cat("Iteration:", i, "\n")
}
  colnames(MC_chain) <- paste(c("G12", "G21", "kappa_r", "kappa_s", "kappa_u", paste("r", 1:time, sep=""), paste("s", 1:12, sep=""), paste("u", 1:ndept, sep=""), paste("B", 1:2, sep=""), "StationaryDistribution"))
  MC_chain<- as.data.frame(MC_chain)
  end_time <- Sys.time()
  time_taken<- end_time - start_time
  print(time_taken)
  return(MC_chain)
}

#malaRes<- MALAInference(y=mixedpop0[[1]], e_it = mixedpop0[[2]], Model = 0, adjmat = sim_adjmat, step_sizes = list("r"=0.000004,"s"=0.03,"u"=0.025), num_iteration = 70000)
#mcmc.plot(malaRes)



