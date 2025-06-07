ModelEvidence.Bridge<- function(y, e_it, adjmat, Model,inf.object, num_samples = 50000){
  logPy<- DetectOutbreaks::ModelEvidence(y=y, e_it = e_it, adjmat = adjmat, Model = Model,inf.object = inf.object, num_samples = num_samples)
  print(logPy)
  R<- -1 * adjmat
  diag(R)<- -rowSums(R, na.rm = T)
  rankdef<- nrow(R)-qr(R)$rank

  design_matrix_func <- get(paste0("DesignMatrixModel", Model))
  z_it <- design_matrix_func(y, adjmat)[[1]]
  z_it2 <- design_matrix_func(y, adjmat)[[2]]

  y<- ifelse(is.na(y), -1, y)

  time<- ncol(y)
  ndept<- nrow(y)

  SMat<- matrix(0, nrow=12, ncol=12)
  SMat[1, ]<- c(2,-1, rep(0, 12-3), -1)
  SMat[2, ]<- c(-1,2,-1, rep(0, 12-3))
  SMat[3, ]<- c(0, -1,2,-1, rep(0, 12-4))
  SMat[(12-1), ]<- c(rep(0, 12-3), -1,2,-1)
  SMat[12, ]<- c(-1, rep(0, 12-3), -1, 2)
  for(i in 3:(12-3)){
    SMat[i+1, ((i):(i+2))]<- c(-1,2,-1)
  }

  #R MCMC samples
  if(is.data.frame(inf.object)){
    if(Model == 0){
      inf.object<- inf.object[,-(c(1,2,(ncol(inf.object)-2):ncol(inf.object)))]
    }else if(Model %in% c(1, 2, 4, 5, 7)){
      inf.object<- inf.object[,-((ncol(inf.object)-1):ncol(inf.object))]
    }else if(Model %in% c(3, 6)){
      inf.object<- inf.object[, -ncol(inf.object)]
    }
  }else{
    #Stan HMC samples
    cond<- checkB(inf.object)
    if(cond){
      inf.object<- inf.object$draws(variables = c("G12", "G21", "kappa_r", "kappa_s", "kappa_u", "r", "s", "uconstrained", "B"))
    }else{
      inf.object<- inf.object$draws(variables = c("G12", "G21", "kappa_r", "kappa_s", "kappa_u", "r", "s", "uconstrained"))
    }
    inf.object<- inf.object[,1,]
    inf.object<- as.data.frame(inf.object)
    if(Model == 0){
      inf.object<- inf.object[,-(1:2)]
    }
  }

  post.s<- nrow(inf.object)/2
  s1<- post.s/(num_samples+post.s)
  s2<- num_samples/(num_samples+post.s)

  mu<- colMeans(inf.object[1:post.s, ])
  varcov<- cov(inf.object[1:post.s, ])
  cond<- check_cholesky(varcov)

  logDensProposal<- numeric(num_samples)
  logDensPosterior<- numeric(post.s)

  if(cond){
    theta <- matrix(nrow = 1, ncol = length(mu))
    class(theta) <- "numeric"
    for(i in 1:num_samples){
      mvnfast::rmvt(n=1, mu = mu, sigma = varcov, df = 3, A = theta)
      if(Model==0){
        if(theta[1] <= 0 || theta[2] <= 0 || theta[3]<= 0){
          logDensProposal[i]<- -Inf
        }else{
          Gammas<- c(0.5, 0.5)
          G_mat<- matrix(c(1-Gammas[1], Gammas[1], Gammas[2], 1-Gammas[2]), nrow = length(Gammas), byrow = T)
          kappaR<- theta[1]
          kappaS<- theta[2]
          kappaU<- theta[3]
          r<- theta[3+(1:time)]
          s<- theta[3+time+(1:12)]
          u<- theta[3+time+12+(1:ndept)]
          ARcoeff<- 0
          logDensProposal[i]<- GeneralLoglikelihood_cpp2(y=y, r=r, s=s, u=u, G_mat, e_it, B = ARcoeff, model = Model, z_it, z_it2)
          + sum(dgamma(c(kappaR, kappaS), shape = c(1, 1), rate = c(0.0001, 0.001), log = TRUE)) +
            dgamma(kappaU, shape = 1, rate = 0.01, log = TRUE) +
            randomwalk2(r, kappaR) +
            seasonalComp2(s, kappaS, SMat) +
            logIGMRF1(u, kappaU, R, rankdef) -
           mvnfast::dmvt(theta, mu = mu, sigma = varcov, df = 3, log = TRUE)
        }
      }else if(Model %in% c(1,2,3,4,5,6,7)){
        if(theta[1] <= 0 || theta[1] >= 1 ||theta[2] <= 0 || theta[2] >= 1 || theta[3] <= 0 || theta[4] <= 0 || theta[5]<= 0){
          logDensProposal[i]<- -Inf
        }else{
          Gammas<- theta[1:2]
          G_mat<- matrix(c(1-Gammas[1], Gammas[1], Gammas[2], 1-Gammas[2]), nrow = length(Gammas), byrow = T)
          kappaR<- theta[3]
          kappaS<- theta[4]
          kappaU<- theta[5]
          r<- theta[5+(1:time)]
          s<- theta[5+time+(1:12)]
          u<- theta[5+time+12+(1:ndept)]
          if(Model %in% c(3, 6)){
            ARcoeff<- theta[5+time+12+ndept+(1:2)]
          }else{
            ARcoeff<- theta[5+time+12+ndept+1]
          }
          logDensProposal[i]<- GeneralLoglikelihood_cpp2(y=y, r=r, s=s, u=u, G_mat, e_it, B = ARcoeff, model = Model, z_it, z_it2)
          + sum(dbeta(Gammas, shape1 = c(2, 2), shape2 = c(2, 2), log = TRUE)) +
            sum(dgamma(c(kappaR, kappaS), shape = c(1, 1), rate = c(0.0001, 0.001), log = TRUE)) +
            dgamma(kappaU, shape = 1, rate = 0.01, log = TRUE) +
            randomwalk2(r, kappaR) +
            seasonalComp2(s, kappaS, SMat) +
            logIGMRF1(u, kappaU, R, rankdef) +
            sum(dgamma(ARcoeff, shape = rep(2, length(ARcoeff)), rate = rep(2, length(ARcoeff)), log = TRUE)) -
           mvnfast::dmvt(theta, mu = mu, sigma = varcov, df = 3, log = TRUE)
        }
      }
    }
    for(j in (post.s+1):nrow(inf.object)){
      theta<- as.numeric(inf.object[j, ])
      if(Model==0){
        Gammas<- c(0.5, 0.5)
        G_mat<- matrix(c(1-Gammas[1], Gammas[1], Gammas[2], 1-Gammas[2]), nrow = length(Gammas), byrow = T)
        kappaR<- theta[1]
        kappaS<- theta[2]
        kappaU<- theta[3]
        r<- theta[3+(1:time)]
        s<- theta[3+time+(1:12)]
        u<- theta[3+time+12+(1:ndept)]
        ARcoeff<- 0
      }else{
        Gammas<- theta[1:2]
        G_mat<- matrix(c(1-Gammas[1], Gammas[1], Gammas[2], 1-Gammas[2]), nrow = length(Gammas), byrow = T)
        kappaR<- theta[3]
        kappaS<- theta[4]
        kappaU<- theta[5]
        r<- theta[5+(1:time)]
        s<- theta[5+time+(1:12)]
        u<- theta[5+time+12+(1:ndept)]
        if(Model %in% c(3, 6)){
          ARcoeff<- theta[5+time+12+ndept+(1:2)]
        }else{
          ARcoeff<- theta[5+time+12+ndept+1]
        }
      }
      logDensPosterior[j-post.s]<- GeneralLoglikelihood_cpp2(y=y, r=r, s=s, u=u, G_mat, e_it, B = ARcoeff, model = Model, z_it, z_it2)
      + sum(dbeta(Gammas, shape1 = c(2, 2), shape2 = c(2, 2), log = TRUE)) +
        sum(dgamma(c(kappaR, kappaS), shape = c(1, 1), rate = c(0.0001, 0.001), log = TRUE)) +
        dgamma(kappaU, shape = 1, rate = 0.01, log = TRUE) +
        randomwalk2(r, kappaR) +
        seasonalComp2(s, kappaS, SMat) +
        logIGMRF1(u, kappaU, R, rankdef) +
        sum(dgamma(ARcoeff, shape = rep(2, length(ARcoeff)), rate = rep(2, length(ARcoeff)), log = TRUE)) -
       mvnfast::dmvt(theta, mu = mu, sigma = varcov, df = 3, log = TRUE)
    }
  }else{
    for(i in 1:num_samples){
      theta<- mvtnorm::rmvt(n=1, delta = mu, sigma = varcov, df = 3)
      if(Model==0){
        if(theta[1] <= 0 || theta[2] <= 0 || theta[3]<= 0){
          logDensProposal[i]<- -Inf
        }else{
          Gammas<- c(0.5, 0.5)
          G_mat<- matrix(c(1-Gammas[1], Gammas[1], Gammas[2], 1-Gammas[2]), nrow = length(Gammas), byrow = T)
          kappaR<- theta[1]
          kappaS<- theta[2]
          kappaU<- theta[3]
          r<- theta[3+(1:time)]
          s<- theta[3+time+(1:12)]
          u<- theta[3+time+12+(1:ndept)]
          ARcoeff<- 0
          logDensProposal[i]<- GeneralLoglikelihood_cpp2(y=y, r=r, s=s, u=u, G_mat, e_it, B = ARcoeff, model = Model, z_it, z_it2)
          + sum(dgamma(c(kappaR, kappaS), shape = c(1, 1), rate = c(0.0001, 0.001), log = TRUE)) +
            dgamma(kappaU, shape = 1, rate = 0.01, log = TRUE) +
            randomwalk2(r, kappaR) +
            seasonalComp2(s, kappaS, SMat) +
            logIGMRF1(u, kappaU, R, rankdef) -
           mvtnorm::dmvt(theta, delta = mu, sigma = varcov, df = 3, log = TRUE)
        }
      }else if(Model %in% c(1,2,3,4,5,6,7)){
        if(theta[1] <= 0 || theta[1] >= 1 ||theta[2] <= 0 || theta[2] >= 1 || theta[3] <= 0 || theta[4] <= 0 || theta[5]<= 0){
          logDensProposal[i]<- -Inf
        }else{
          Gammas<- theta[1:2]
          G_mat<- matrix(c(1-Gammas[1], Gammas[1], Gammas[2], 1-Gammas[2]), nrow = length(Gammas), byrow = T)
          kappaR<- theta[3]
          kappaS<- theta[4]
          kappaU<- theta[5]
          r<- theta[5+(1:time)]
          s<- theta[5+time+(1:12)]
          u<- theta[5+time+12+(1:ndept)]
          if(Model %in% c(3, 6)){
            ARcoeff<- theta[5+time+12+ndept+(1:2)]
          }else{
            ARcoeff<- theta[5+time+12+ndept+1]
          }
          logDensProposal[i]<- GeneralLoglikelihood_cpp2(y=y, r=r, s=s, u=u, G_mat, e_it, B = ARcoeff, model = Model, z_it, z_it2)
          + sum(dbeta(Gammas, shape1 = c(2, 2), shape2 = c(2, 2), log = TRUE)) +
            sum(dgamma(c(kappaR, kappaS), shape = c(1, 1), rate = c(0.0001, 0.001), log = TRUE)) +
            dgamma(kappaU, shape = 1, rate = 0.01, log = TRUE) +
            randomwalk2(r, kappaR) +
            seasonalComp2(s, kappaS, SMat) +
            logIGMRF1(u, kappaU, R, rankdef) +
            sum(dgamma(ARcoeff, shape = rep(2, length(ARcoeff)), rate = rep(2, length(ARcoeff)), log = TRUE)) -
           mvtnorm::dmvt(theta, delta = mu, sigma = varcov, df = 3, log = TRUE)
        }
      }
    }
    for(j in (post.s+1):nrow(inf.object)){
      theta<- as.numeric(inf.object[j, ])
      if(Model==0){
        Gammas<- c(0.5, 0.5)
        G_mat<- matrix(c(1-Gammas[1], Gammas[1], Gammas[2], 1-Gammas[2]), nrow = length(Gammas), byrow = T)
        kappaR<- theta[1]
        kappaS<- theta[2]
        kappaU<- theta[3]
        r<- theta[3+(1:time)]
        s<- theta[3+time+(1:12)]
        u<- theta[3+time+12+(1:ndept)]
        ARcoeff<- 0
      }else{
        Gammas<- theta[1:2]
        G_mat<- matrix(c(1-Gammas[1], Gammas[1], Gammas[2], 1-Gammas[2]), nrow = length(Gammas), byrow = T)
        kappaR<- theta[3]
        kappaS<- theta[4]
        kappaU<- theta[5]
        r<- theta[5+(1:time)]
        s<- theta[5+time+(1:12)]
        u<- theta[5+time+12+(1:ndept)]
        if(Model %in% c(3, 6)){
          ARcoeff<- theta[5+time+12+ndept+(1:2)]
        }else{
          ARcoeff<- theta[5+time+12+ndept+1]
        }
      }
      logDensPosterior[j-post.s]<- GeneralLoglikelihood_cpp2(y=y, r=r, s=s, u=u, G_mat, e_it, B = ARcoeff, model = Model, z_it, z_it2)
      + sum(dbeta(Gammas, shape1 = c(2, 2), shape2 = c(2, 2), log = TRUE)) +
        sum(dgamma(c(kappaR, kappaS), shape = c(1, 1), rate = c(0.0001, 0.001), log = TRUE)) +
        dgamma(kappaU, shape = 1, rate = 0.01, log = TRUE) +
        randomwalk2(r, kappaR) +
        seasonalComp2(s, kappaS, SMat) +
        logIGMRF1(u, kappaU, R, rankdef) +
        sum(dgamma(ARcoeff, shape = rep(2, length(ARcoeff)), rate = rep(2, length(ARcoeff)), log = TRUE)) -
       mvtnorm::dmvt(theta, delta = mu, sigma = varcov, df = 3, log = TRUE)
    }
  }

  #print(logDensProposal)
  #print(logDensPosterior)
  log_w<- logDensProposal[logDensProposal != -Inf] - logSumExp_cpp(logDensProposal[logDensProposal != -Inf])
  propESS<- exp(-logSumExp_cpp(2 * log_w))
  postESS<- exp((2*logSumExp_cpp(logDensPosterior))-(logSumExp_cpp(2*logDensPosterior)))
  print(paste("ESS from proposal samples is ", propESS))
  print(paste("ESS from posterior samples is ", postESS))

  numerator<- numeric(num_samples)
  denominator<- numeric(post.s)
  delta<- 0.03

  while(delta > 0.01){
    for(i in 1:num_samples){
      numerator[i]<- logDensProposal[i] - logSumExp_cpp(c(log(s1)+logDensProposal[i], log(s2) + logPy))
    }
    for(j in (post.s+1):nrow(inf.object)){
      denominator[j-post.s]<- log(1) - logSumExp_cpp(c(log(s1)+logDensPosterior[j-post.s], log(s2) + logPy))
    }
    numerator<- numerator[numerator != -Inf]
    denominator<- denominator[denominator != -Inf]
    currentlogPy<- (-log(length(numerator))+logSumExp_cpp(numerator))-(-log(post.s)+logSumExp_cpp(denominator))
    delta<- abs(logPy-currentlogPy)
    print(currentlogPy)
    logPy<- currentlogPy
  }
  #return(list(logDensProposal, logDensPosterior))
  return(currentlogPy)
}



#library(doParallel)

# Register parallel backend
#cores <- detectCores() - 1  # Use one less than available cores
#registerDoParallel(cores)
# Perform parallel computation
#results <- foreach(i = 1:10, .combine = 'c', .export = c("y", "e_it", "sim_adjmat", "Mod0"), .packages = 'DetectOutbreaks') %dopar% {
  # Simulate a computation-heavy task
#  Sys.sleep(1)  # Sleep for 1 second
#  DetectOutbreaks::ModelEvidence(y=y, e_it = e_it, adjmat = sim_adjmat, Model = 0,inf.object = Mod0, num_samples = 500)
#}

# Stop the parallel backend
#stopImplicitCluster()

#print(results)

ParallelModelEvidence.Bridge<- function(y, e_it, adjmat, Model,inf.object, num_samples = 50000){

  logPy<- DetectOutbreaks::ModelEvidence(y=y, e_it = e_it, adjmat = adjmat, Model = Model,inf.object = inf.object, num_samples = num_samples)
  print(logPy)
  R<- -1 * adjmat
  diag(R)<- -rowSums(R, na.rm = T)
  rankdef<- nrow(R)-qr(R)$rank

  design_matrix_func <- get(paste0("DesignMatrixModel", Model))
  z_it <- design_matrix_func(y, adjmat)[[1]]
  z_it2 <- design_matrix_func(y, adjmat)[[2]]

  y<- ifelse(is.na(y), -1, y)

  time<- ncol(y)
  ndept<- nrow(y)

  SMat<- matrix(0, nrow=12, ncol=12)
  SMat[1, ]<- c(2,-1, rep(0, 12-3), -1)
  SMat[2, ]<- c(-1,2,-1, rep(0, 12-3))
  SMat[3, ]<- c(0, -1,2,-1, rep(0, 12-4))
  SMat[(12-1), ]<- c(rep(0, 12-3), -1,2,-1)
  SMat[12, ]<- c(-1, rep(0, 12-3), -1, 2)
  for(i in 3:(12-3)){
    SMat[i+1, ((i):(i+2))]<- c(-1,2,-1)
  }

  #R MCMC samples
  if(is.data.frame(inf.object)){
    if(Model == 0){
      inf.object<- inf.object[,-(c(1,2,(ncol(inf.object)-2):ncol(inf.object)))]
    }else if(Model %in% c(1, 2, 4, 5, 7)){
      inf.object<- inf.object[,-((ncol(inf.object)-1):ncol(inf.object))]
    }else if(Model %in% c(3, 6)){
      inf.object<- inf.object[, -ncol(inf.object)]
    }
  }else{
    #Stan HMC samples
    cond<- checkB(inf.object)
    if(cond){
      inf.object<- inf.object$draws(variables = c("G12", "G21", "kappa_r", "kappa_s", "kappa_u", "r", "s", "uconstrained", "B"))
    }else{
      inf.object<- inf.object$draws(variables = c("G12", "G21", "kappa_r", "kappa_s", "kappa_u", "r", "s", "uconstrained"))
    }
    inf.object<- inf.object[,1,]
    inf.object<- as.data.frame(inf.object)
    if(Model == 0){
      inf.object<- inf.object[,-(1:2)]
    }
  }

  post.s<- nrow(inf.object)/2
  s1<- post.s/(num_samples+post.s)
  s2<- num_samples/(num_samples+post.s)

  mu<- colMeans(inf.object[1:post.s, ])
  varcov<- cov(inf.object[1:post.s, ])
  cond<- check_cholesky(varcov)

  logDensProposal<- numeric(num_samples)
  logDensPosterior<- numeric(post.s)

  library(doParallel)
  cores<- detectCores()-1
  cl <- makeCluster(cores)
  clusterEvalQ(cl, devtools::load_all("C:/Users/Matthew Adeoye/Documents/PhD Statistics/Year 3/SpatMet"))
  registerDoParallel(cl)

  if(cond){
    logDensProposal<- foreach(i = 1:num_samples, .combine = c, .export = 'logI', .packages = 'mvnfast') %dopar% {
      theta<- mvnfast::rmvt(n = 1, mu = mu, sigma = varcov, df = 3)
      logI(theta=theta, Model=Model, speed=1, y=y, e_it=e_it, z_it=z_it, z_it2=z_it2,
                                SMat=SMat, R=R, rankdef=rankdef, mu=mu, varcov=varcov)
     }

    logDensPosterior<- foreach(j = (post.s+1):nrow(inf.object), .combine = 'c', .export = 'logI', .packages = 'mvnfast') %dopar% {
        theta.post <- as.numeric(inf.object[j, ])
         logI(theta=theta.post, Model=Model, speed=1, y=y, e_it=e_it, z_it=z_it, z_it2=z_it2,
             SMat=SMat, R=R, rankdef=rankdef, mu=mu, varcov=varcov)
    }
  }else{
    logDensProposal<- foreach(i = 1:num_samples, .combine = c, .export = 'logI', .packages = 'mvtnorm') %dopar% {
      theta<- mvtnorm::rmvt(n=1, delta = mu, sigma = varcov, df = 3)
      logI(theta=theta, Model=Model,speed=0, y=y, e_it=e_it, z_it=z_it, z_it2=z_it2, SMat=SMat, R=R, rankdef=rankdef, mu=mu, varcov=varcov)
    }
    logDensPosterior<- foreach(j = (post.s+1):nrow(inf.object), .combine = c, .export = 'logI', .packages = 'mvtnorm') %dopar% {
      theta.post<- as.numeric(inf.object[j, ])
      logI(theta=theta.post, Model=Model,speed=0, y=y, e_it=e_it, z_it=z_it, z_it2=z_it2, SMat=SMat, R=R, rankdef=rankdef, mu=mu, varcov=varcov)
    }
  }
  stopCluster(cl)

  numerator<- numeric(num_samples)
  denominator<- numeric(post.s)

  delta<- 0.03

  while(delta > 0.01){
    for(i in 1:num_samples){
      numerator[i]<- logDensProposal[i] - logSumExp_cpp(c(log(s1)+logDensProposal[i], log(s2) + logPy))
    }
    for(j in (post.s+1):nrow(inf.object)){
      theta.post<- as.numeric(inf.object[j, ])
      denominator[j-post.s]<- log(1) - logSumExp_cpp(c(log(s1)+logDensPosterior[j-post.s], log(s2) + logPy))
    }
    numerator<- numerator[numerator != -Inf]
    denominator<- denominator[denominator != -Inf]
    currentlogPy<- (log(1)-log(num_samples)+logSumExp_cpp(numerator))-(log(1)-log(post.s)+logSumExp_cpp(denominator))
    delta<- abs(logPy-currentlogPy)
    print(currentlogPy)
    logPy<- currentlogPy
  }
  return(currentlogPy)
}


ModelEvidence.IS<- function(y, e_it, adjmat, Model, inf.object, num_samples = 50000){

  R<- -1 * adjmat
  diag(R)<- -rowSums(R, na.rm = T)
  rankdef<- nrow(R)-qr(R)$rank

  design_matrix_func <- get(paste0("DesignMatrixModel", Model))
  z_it <- design_matrix_func(y, adjmat)[[1]]
  z_it2 <- design_matrix_func(y, adjmat)[[2]]

  y<- ifelse(is.na(y), -1, y)

  time<- ncol(y)
  ndept<- nrow(y)

  SMat<- matrix(0, nrow=12, ncol=12)
  SMat[1, ]<- c(2,-1, rep(0, 12-3), -1)
  SMat[2, ]<- c(-1,2,-1, rep(0, 12-3))
  SMat[3, ]<- c(0, -1,2,-1, rep(0, 12-4))
  SMat[(12-1), ]<- c(rep(0, 12-3), -1,2,-1)
  SMat[12, ]<- c(-1, rep(0, 12-3), -1, 2)
  for(i in 3:(12-3)){
    SMat[i+1, ((i):(i+2))]<- c(-1,2,-1)
  }

  #R MCMC samples
  if(is.data.frame(inf.object)){
    if(Model == 0){
      inf.object<- inf.object[,-(c(1,2,(ncol(inf.object)-2):ncol(inf.object)))]
    }else if(Model %in% c(1, 2, 4, 5, 7)){
      inf.object<- inf.object[,-((ncol(inf.object)-1):ncol(inf.object))]
    }else if(Model %in% c(3, 6)){
      inf.object<- inf.object[, -ncol(inf.object)]
    }
  }else{
    #Stan HMC samples
    cond<- checkB(inf.object)
    if(cond){
      inf.object<- inf.object$draws(variables = c("G12", "G21", "kappa_r", "kappa_s", "kappa_u", "r", "s", "uconstrained", "B"))
    }else{
      inf.object<- inf.object$draws(variables = c("G12", "G21", "kappa_r", "kappa_s", "kappa_u", "r", "s", "uconstrained"))
    }
    inf.object<- inf.object[,1,]
    inf.object<- as.data.frame(inf.object)
    if(Model == 0){
      inf.object<- inf.object[,-(1:2)]
    }
  }

  mu<- colMeans(inf.object)
  varcov<- cov(inf.object)
  cond<- check_cholesky(varcov)
  a<- numeric(num_samples)

  if(cond){
    theta <- matrix(nrow = 1, ncol = length(mu))
    class(theta) <- "numeric"
    for(i in 1:num_samples){
      mvnfast::rmvt(n=1, mu = mu, sigma = varcov, df = 3, A = theta)
      if(Model==0){
        if(theta[1] <= 0 || theta[2] <= 0 || theta[3]<= 0){
          a[i]<- -Inf
        }else{
          Gammas<- c(0.5, 0.5)
          G_mat<- matrix(c(1-Gammas[1], Gammas[1], Gammas[2], 1-Gammas[2]), nrow = length(Gammas), byrow = T)
          kappaR<- theta[1]
          kappaS<- theta[2]
          kappaU<- theta[3]
          r<- theta[3+(1:time)]
          s<- theta[3+time+(1:12)]
          u<- theta[3+time+12+(1:ndept)]
          ARcoeff<- 0
          a[i]<- GeneralLoglikelihood_cpp2(y=y, r=r, s=s, u=u, G_mat, e_it, B = ARcoeff, model = Model, z_it, z_it2)
          + sum(dgamma(c(kappaR, kappaS), shape = c(1, 1), rate = c(0.0001, 0.001), log = TRUE)) +
            dgamma(kappaU, shape = 1, rate = 0.01, log = TRUE) +
            randomwalk2(r, kappaR) +
            seasonalComp2(s, kappaS, SMat) +
            logIGMRF1(u, kappaU, R, rankdef)
          - mvnfast::dmvt(theta, mu = mu, sigma = varcov, df = 3, log = TRUE)
        }
      }else if(Model %in% c(1,2,3,4,5,6,7)){
        if(theta[1] <= 0 || theta[1] >= 1 ||theta[2] <= 0 || theta[2] >= 1 || theta[3] <= 0 || theta[4] <= 0 || theta[5]<= 0){
          a[i]<- -Inf
        }else{
          Gammas<- theta[1:2]
          G_mat<- matrix(c(1-Gammas[1], Gammas[1], Gammas[2], 1-Gammas[2]), nrow = length(Gammas), byrow = T)
          kappaR<- theta[3]
          kappaS<- theta[4]
          kappaU<- theta[5]
          r<- theta[5+(1:time)]
          s<- theta[5+time+(1:12)]
          u<- theta[5+time+12+(1:ndept)]
          if(Model %in% c(3, 6)){
            ARcoeff<- theta[5+time+12+ndept+(1:2)]
          }else{
            ARcoeff<- theta[5+time+12+ndept+1]
          }
          a[i]<- GeneralLoglikelihood_cpp2(y=y, r=r, s=s, u=u, G_mat, e_it, B = ARcoeff, model = Model, z_it, z_it2)
          + sum(dbeta(Gammas, shape1 = c(2, 2), shape2 = c(2, 2), log = TRUE)) +
            sum(dgamma(c(kappaR, kappaS), shape = c(1, 1), rate = c(0.0001, 0.001), log = TRUE)) +
            dgamma(kappaU, shape = 1, rate = 0.01, log = TRUE) +
            randomwalk2(r, kappaR) +
            seasonalComp2(s, kappaS, SMat) +
            logIGMRF1(u, kappaU, R, rankdef) +
            sum(dgamma(ARcoeff, shape = rep(2, length(ARcoeff)), rate = rep(2, length(ARcoeff)), log = TRUE))
          - mvnfast::dmvt(theta, mu = mu, sigma = varcov, df = 3, log = TRUE)
        }
      }
    }
  }else{
    for(i in 1:num_samples){
      theta<- mvtnorm::rmvt(n=1, delta = mu, sigma = varcov, df = 3)
      if(Model==0){
        if(theta[1] <= 0 || theta[2] <= 0 || theta[3]<= 0){
          a[i]<- -Inf
        }else{
          Gammas<- c(0.5, 0.5)
          G_mat<- matrix(c(1-Gammas[1], Gammas[1], Gammas[2], 1-Gammas[2]), nrow = length(Gammas), byrow = T)
          kappaR<- theta[1]
          kappaS<- theta[2]
          kappaU<- theta[3]
          r<- theta[3+(1:time)]
          s<- theta[3+time+(1:12)]
          u<- theta[3+time+12+(1:ndept)]
          ARcoeff<- 0
          a[i]<- GeneralLoglikelihood_cpp2(y=y, r=r, s=s, u=u, G_mat, e_it, B = ARcoeff, model = Model, z_it, z_it2)
          + sum(dgamma(c(kappaR, kappaS), shape = c(1, 1), rate = c(0.0001, 0.001), log = TRUE)) +
            dgamma(kappaU, shape = 1, rate = 0.01, log = TRUE) +
            randomwalk2(r, kappaR) +
            seasonalComp2(s, kappaS, SMat) +
            logIGMRF1(u, kappaU, R, rankdef)
          - mvtnorm::dmvt(theta, delta = mu, sigma = varcov, df = 3, log = TRUE)
        }
      }else if(Model %in% c(1,2,3,4,5,6,7)){
        if(theta[1] <= 0 || theta[1] >= 1 ||theta[2] <= 0 || theta[2] >= 1 || theta[3] <= 0 || theta[4] <= 0 || theta[5]<= 0){
          a[i]<- -Inf
        }else{
          Gammas<- theta[1:2]
          G_mat<- matrix(c(1-Gammas[1], Gammas[1], Gammas[2], 1-Gammas[2]), nrow = length(Gammas), byrow = T)
          kappaR<- theta[3]
          kappaS<- theta[4]
          kappaU<- theta[5]
          r<- theta[5+(1:time)]
          s<- theta[5+time+(1:12)]
          u<- theta[5+time+12+(1:ndept)]
          if(Model %in% c(3, 6)){
            ARcoeff<- theta[5+time+12+ndept+(1:2)]
          }else{
            ARcoeff<- theta[5+time+12+ndept+1]
          }
          a[i]<- GeneralLoglikelihood_cpp2(y=y, r=r, s=s, u=u, G_mat, e_it, B = ARcoeff, model = Model, z_it, z_it2)
          + sum(dbeta(Gammas, shape1 = c(2, 2), shape2 = c(2, 2), log = TRUE)) +
            sum(dgamma(c(kappaR, kappaS), shape = c(1, 1), rate = c(0.0001, 0.001), log = TRUE)) +
            dgamma(kappaU, shape = 1, rate = 0.01, log = TRUE) +
            randomwalk2(r, kappaR) +
            seasonalComp2(s, kappaS, SMat) +
            logIGMRF1(u, kappaU, R, rankdef) +
            sum(dgamma(ARcoeff, shape = rep(2, length(ARcoeff)), rate = rep(2, length(ARcoeff)), log = TRUE))
          - mvtnorm::dmvt(theta, delta = mu, sigma = varcov, df = 3, log = TRUE)
        }
      }
    }
  }
  #print(a)
  validDensities<- a[a != -Inf]
  log_w<- validDensities - logSumExp_cpp(validDensities)
  nESS<- exp(-logSumExp_cpp(2 * log_w))
  uESS<- exp((2*logSumExp_cpp(validDensities))-(logSumExp_cpp(2*validDensities)))
  print(paste("n.ESS is ", nESS))
  print(paste("un.ESS is ", uESS))
  MarginalLikelihood<- log(1) - log(length(validDensities)) + logSumExp_cpp(validDensities)
  return(MarginalLikelihood)
}


#Alternative implementation 1
ModelEvidence.Bridge1<- function(y, e_it, adjmat, Model,inf.object, num_samples = 50000){
  rt<- 0.0000001
  R<- -1 * adjmat
  diag(R)<- -rowSums(R, na.rm = T)
  rankdef<- nrow(R)-qr(R)$rank

  design_matrix_func <- get(paste0("DesignMatrixModel", Model))
  z_it <- design_matrix_func(y, adjmat)[[1]]
  z_it2 <- design_matrix_func(y, adjmat)[[2]]

  y<- ifelse(is.na(y), -1, y)

  time<- ncol(y)
  ndept<- nrow(y)

  SMat<- matrix(0, nrow=12, ncol=12)
  SMat[1, ]<- c(2,-1, rep(0, 12-3), -1)
  SMat[2, ]<- c(-1,2,-1, rep(0, 12-3))
  SMat[3, ]<- c(0, -1,2,-1, rep(0, 12-4))
  SMat[(12-1), ]<- c(rep(0, 12-3), -1,2,-1)
  SMat[12, ]<- c(-1, rep(0, 12-3), -1, 2)
  for(i in 3:(12-3)){
    SMat[i+1, ((i):(i+2))]<- c(-1,2,-1)
  }

  #R MCMC samples
  if(is.data.frame(inf.object)){
    if(Model == 0){
      inf.object<- inf.object[,-(c(1,2,(ncol(inf.object)-2):ncol(inf.object)))]
    }else if(Model %in% c(1, 2, 4, 5, 7)){
      inf.object<- inf.object[,-((ncol(inf.object)-1):ncol(inf.object))]
    }else if(Model %in% c(3, 6)){
      inf.object<- inf.object[, -ncol(inf.object)]
    }
  }else{
    #Stan HMC samples
    cond<- checkB(inf.object)
    if(cond){
      inf.object<- inf.object$draws(variables = c("G12", "G21", "kappa_r", "kappa_s", "kappa_u", "r", "s", "uconstrained", "B"))
    }else{
      inf.object<- inf.object$draws(variables = c("G12", "G21", "kappa_r", "kappa_s", "kappa_u", "r", "s", "uconstrained"))
    }
    inf.object<- inf.object[,1,]
    inf.object<- as.data.frame(inf.object)
    if(Model == 0){
      inf.object<- inf.object[,-(1:2)]
    }
  }

  post.s<- nrow(inf.object)/2
  s1<- post.s/(num_samples+post.s)
  s2<- num_samples/(num_samples+post.s)

  mu<- colMeans(inf.object[1:post.s, ])
  varcov<- cov(inf.object[1:post.s, ])
  cond<- check_cholesky(varcov)

  logDensProposal<- numeric(num_samples)
  logDensPosterior<- numeric(post.s)

  if(cond){
    theta <- matrix(nrow = 1, ncol = length(mu))
    class(theta) <- "numeric"
    for(i in 1:num_samples){
      mvnfast::rmvt(n=1, mu = mu, sigma = varcov, df = 3, A = theta)
      if(Model==0){
        if(theta[1] <= 0 || theta[2] <= 0 || theta[3]<= 0){
          logDensProposal[i]<- -Inf
        }else{
          Gammas<- c(0.5, 0.5)
          G_mat<- matrix(c(1-Gammas[1], Gammas[1], Gammas[2], 1-Gammas[2]), nrow = length(Gammas), byrow = T)
          kappaR<- theta[1]
          kappaS<- theta[2]
          kappaU<- theta[3]
          r<- theta[3+(1:time)]
          s<- theta[3+time+(1:12)]
          u<- theta[3+time+12+(1:ndept)]
          ARcoeff<- 0
          logDensProposal[i]<- GeneralLoglikelihood_cpp2(y=y, r=r, s=s, u=u, G_mat, e_it, B = ARcoeff, model = Model, z_it, z_it2)
          + sum(dgamma(c(kappaR, kappaS), shape = c(1, 1), rate = c(0.0001, 0.001), log = TRUE)) +
            dgamma(kappaU, shape = 1, rate = 0.01, log = TRUE) +
            randomwalk2(r, kappaR) +
            seasonalComp2(s, kappaS, SMat) +
            logIGMRF1(u, kappaU, R, rankdef)
          - mvnfast::dmvt(theta, mu = mu, sigma = varcov, df = 3, log = TRUE)
        }
      }else if(Model %in% c(1,2,3,4,5,6,7)){
        if(theta[1] <= 0 || theta[1] >= 1 ||theta[2] <= 0 || theta[2] >= 1 || theta[3] <= 0 || theta[4] <= 0 || theta[5]<= 0){
          logDensProposal[i]<- -Inf
        }else{
          Gammas<- theta[1:2]
          G_mat<- matrix(c(1-Gammas[1], Gammas[1], Gammas[2], 1-Gammas[2]), nrow = length(Gammas), byrow = T)
          kappaR<- theta[3]
          kappaS<- theta[4]
          kappaU<- theta[5]
          r<- theta[5+(1:time)]
          s<- theta[5+time+(1:12)]
          u<- theta[5+time+12+(1:ndept)]
          if(Model %in% c(3, 6)){
            ARcoeff<- theta[5+time+12+ndept+(1:2)]
          }else{
            ARcoeff<- theta[5+time+12+ndept+1]
          }
          logDensProposal[i]<- GeneralLoglikelihood_cpp2(y=y, r=r, s=s, u=u, G_mat, e_it, B = ARcoeff, model = Model, z_it, z_it2)
          + sum(dbeta(Gammas, shape1 = c(2, 2), shape2 = c(2, 2), log = TRUE)) +
            sum(dgamma(c(kappaR, kappaS), shape = c(1, 1), rate = c(0.0001, 0.001), log = TRUE)) +
            dgamma(kappaU, shape = 1, rate = 0.01, log = TRUE) +
            randomwalk2(r, kappaR) +
            seasonalComp2(s, kappaS, SMat) +
            logIGMRF1(u, kappaU, R, rankdef) +
            sum(dgamma(ARcoeff, shape = rep(2, length(ARcoeff)), rate = rep(2, length(ARcoeff)), log = TRUE))
          - mvnfast::dmvt(theta, mu = mu, sigma = varcov, df = 3, log = TRUE)
        }
      }
    }
    for(j in (post.s+1):nrow(inf.object)){
      theta<- as.numeric(inf.object[j, ])
      if(Model==0){
        Gammas<- c(0.5, 0.5)
        G_mat<- matrix(c(1-Gammas[1], Gammas[1], Gammas[2], 1-Gammas[2]), nrow = length(Gammas), byrow = T)
        kappaR<- theta[1]
        kappaS<- theta[2]
        kappaU<- theta[3]
        r<- theta[3+(1:time)]
        s<- theta[3+time+(1:12)]
        u<- theta[3+time+12+(1:ndept)]
        ARcoeff<- 0
      }else{
        Gammas<- theta[1:2]
        G_mat<- matrix(c(1-Gammas[1], Gammas[1], Gammas[2], 1-Gammas[2]), nrow = length(Gammas), byrow = T)
        kappaR<- theta[3]
        kappaS<- theta[4]
        kappaU<- theta[5]
        r<- theta[5+(1:time)]
        s<- theta[5+time+(1:12)]
        u<- theta[5+time+12+(1:ndept)]
        if(Model %in% c(3, 6)){
          ARcoeff<- theta[5+time+12+ndept+(1:2)]
        }else{
          ARcoeff<- theta[5+time+12+ndept+1]
        }
      }
      logDensPosterior[j-post.s]<- GeneralLoglikelihood_cpp2(y=y, r=r, s=s, u=u, G_mat, e_it, B = ARcoeff, model = Model, z_it, z_it2)
      + sum(dbeta(Gammas, shape1 = c(2, 2), shape2 = c(2, 2), log = TRUE)) +
        sum(dgamma(c(kappaR, kappaS), shape = c(1, 1), rate = c(0.0001, 0.001), log = TRUE)) +
        dgamma(kappaU, shape = 1, rate = 0.01, log = TRUE) +
        randomwalk2(r, kappaR) +
        seasonalComp2(s, kappaS, SMat) +
        logIGMRF1(u, kappaU, R, rankdef) +
        sum(dgamma(ARcoeff, shape = rep(2, length(ARcoeff)), rate = rep(2, length(ARcoeff)), log = TRUE))
      - mvnfast::dmvt(theta, mu = mu, sigma = varcov, df = 3, log = TRUE)
    }
  }else{
    for(i in 1:num_samples){
      theta<- mvtnorm::rmvt(n=1, delta = mu, sigma = varcov, df = 3)
      if(Model==0){
        if(theta[1] <= 0 || theta[2] <= 0 || theta[3]<= 0){
          logDensProposal[i]<- -Inf
        }else{
          Gammas<- c(0.5, 0.5)
          G_mat<- matrix(c(1-Gammas[1], Gammas[1], Gammas[2], 1-Gammas[2]), nrow = length(Gammas), byrow = T)
          kappaR<- theta[1]
          kappaS<- theta[2]
          kappaU<- theta[3]
          r<- theta[3+(1:time)]
          s<- theta[3+time+(1:12)]
          u<- theta[3+time+12+(1:ndept)]
          ARcoeff<- 0
          logDensProposal[i]<- GeneralLoglikelihood_cpp2(y=y, r=r, s=s, u=u, G_mat, e_it, B = ARcoeff, model = Model, z_it, z_it2)
          + sum(dgamma(c(kappaR, kappaS), shape = c(1, 1), rate = c(0.0001, 0.001), log = TRUE)) +
            dgamma(kappaU, shape = 1, rate = 0.01, log = TRUE) +
            randomwalk2(r, kappaR) +
            seasonalComp2(s, kappaS, SMat) +
            logIGMRF1(u, kappaU, R, rankdef)
          - mvtnorm::dmvt(theta, delta = mu, sigma = varcov, df = 3, log = TRUE)
        }
      }else if(Model %in% c(1,2,3,4,5,6,7)){
        if(theta[1] <= 0 || theta[1] >= 1 ||theta[2] <= 0 || theta[2] >= 1 || theta[3] <= 0 || theta[4] <= 0 || theta[5]<= 0){
          logDensProposal[i]<- -Inf
        }else{
          Gammas<- theta[1:2]
          G_mat<- matrix(c(1-Gammas[1], Gammas[1], Gammas[2], 1-Gammas[2]), nrow = length(Gammas), byrow = T)
          kappaR<- theta[3]
          kappaS<- theta[4]
          kappaU<- theta[5]
          r<- theta[5+(1:time)]
          s<- theta[5+time+(1:12)]
          u<- theta[5+time+12+(1:ndept)]
          if(Model %in% c(3, 6)){
            ARcoeff<- theta[5+time+12+ndept+(1:2)]
          }else{
            ARcoeff<- theta[5+time+12+ndept+1]
          }
          logDensProposal[i]<- GeneralLoglikelihood_cpp2(y=y, r=r, s=s, u=u, G_mat, e_it, B = ARcoeff, model = Model, z_it, z_it2)
          + sum(dbeta(Gammas, shape1 = c(2, 2), shape2 = c(2, 2), log = TRUE)) +
            sum(dgamma(c(kappaR, kappaS), shape = c(1, 1), rate = c(0.0001, 0.001), log = TRUE)) +
            dgamma(kappaU, shape = 1, rate = 0.01, log = TRUE) +
            randomwalk2(r, kappaR) +
            seasonalComp2(s, kappaS, SMat) +
            logIGMRF1(u, kappaU, R, rankdef) +
            sum(dgamma(ARcoeff, shape = rep(2, length(ARcoeff)), rate = rep(2, length(ARcoeff)), log = TRUE))
          - mvtnorm::dmvt(theta, delta = mu, sigma = varcov, df = 3, log = TRUE)
        }
      }
    }
    for(j in (post.s+1):nrow(inf.object)){
      theta<- as.numeric(inf.object[j, ])
      if(Model==0){
        Gammas<- c(0.5, 0.5)
        G_mat<- matrix(c(1-Gammas[1], Gammas[1], Gammas[2], 1-Gammas[2]), nrow = length(Gammas), byrow = T)
        kappaR<- theta[1]
        kappaS<- theta[2]
        kappaU<- theta[3]
        r<- theta[3+(1:time)]
        s<- theta[3+time+(1:12)]
        u<- theta[3+time+12+(1:ndept)]
        ARcoeff<- 0
      }else{
        Gammas<- theta[1:2]
        G_mat<- matrix(c(1-Gammas[1], Gammas[1], Gammas[2], 1-Gammas[2]), nrow = length(Gammas), byrow = T)
        kappaR<- theta[3]
        kappaS<- theta[4]
        kappaU<- theta[5]
        r<- theta[5+(1:time)]
        s<- theta[5+time+(1:12)]
        u<- theta[5+time+12+(1:ndept)]
        if(Model %in% c(3, 6)){
          ARcoeff<- theta[5+time+12+ndept+(1:2)]
        }else{
          ARcoeff<- theta[5+time+12+ndept+1]
        }
      }
      logDensPosterior[j-post.s]<- GeneralLoglikelihood_cpp2(y=y, r=r, s=s, u=u, G_mat, e_it, B = ARcoeff, model = Model, z_it, z_it2)
      + sum(dbeta(Gammas, shape1 = c(2, 2), shape2 = c(2, 2), log = TRUE)) +
        sum(dgamma(c(kappaR, kappaS), shape = c(1, 1), rate = c(0.0001, 0.001), log = TRUE)) +
        dgamma(kappaU, shape = 1, rate = 0.01, log = TRUE) +
        randomwalk2(r, kappaR) +
        seasonalComp2(s, kappaS, SMat) +
        logIGMRF1(u, kappaU, R, rankdef) +
        sum(dgamma(ARcoeff, shape = rep(2, length(ARcoeff)), rate = rep(2, length(ARcoeff)), log = TRUE))
      - mvtnorm::dmvt(theta, delta = mu, sigma = varcov, df = 3, log = TRUE)
    }
  }

  #print(logDensProposal)
  #print(logDensPosterior)
  log_w<- logDensProposal[logDensProposal != -Inf] - logSumExp_cpp(logDensProposal[logDensProposal != -Inf])
  propESS<- exp(-logSumExp_cpp(2 * log_w))
  postESS<- exp((2*logSumExp_cpp(logDensPosterior))-(logSumExp_cpp(2*logDensPosterior)))
  print(paste("ESS from proposal samples is ", propESS))
  print(paste("ESS from posterior samples is ", postESS))

  numerator<- numeric(num_samples)
  denominator<- numeric(post.s)

  medL<- median(logDensPosterior)
  delta<- 0.3

  while(delta > 0.01){
    for(i in 1:num_samples){
      numerator[i]<- exp(logDensProposal[i] - medL)/(s1 * exp(logDensProposal[i] - medL) + s2 * rt)
    }
    for(j in 1:post.s){
      denominator[j]<- 1/(s1 * exp(logDensPosterior[j] - medL) + s2 * rt)
    }
    numerator<- numerator[numerator != -Inf]
    current.rt<- (sum(numerator)/length(numerator))/(sum(denominator)/post.s)
    delta<- abs(rt-current.rt)
    print(current.rt)
    rt<- current.rt
  }
  current.rt<- log(current.rt) + medL
  return(current.rt)
}

