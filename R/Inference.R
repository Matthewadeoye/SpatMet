multstrainInfer<- function(y, e_it, nstrain, Model, adjmat, independentChains, num_iteration = 30000, Stan = TRUE, GPU = FALSE, nchains = 4,
                 iter = 4000, seed = NULL, verbose = F, ModEvid = F, OutbreakProb = F, adaptdelta = 0.90, Burn.in = 1000){
  start_time <- Sys.time()

  R<- -1 * adjmat
  diag(R)<- -rowSums(R, na.rm = T)
  rankdef<- nrow(R)-qr(R)$rank
  ndept<- nrow(e_it)
  time<- ncol(e_it)

  original.y<- y

  #flag missing data for inference
  y<- ifelse(is.na(y), -1, y)

  yflat<- as.numeric(aperm(y, c(2,1,3)))

  RW1PrecMat<- matrix(0, nrow=12, ncol=12)
  RW1PrecMat[1, ]<- c(2,-1, rep(0, 12-3), -1)
  RW1PrecMat[2, ]<- c(-1,2,-1, rep(0, 12-3))
  RW1PrecMat[3, ]<- c(0, -1,2,-1, rep(0, 12-4))
  RW1PrecMat[(12-1), ]<- c(rep(0, 12-3), -1,2,-1)
  RW1PrecMat[12, ]<- c(-1, rep(0, 12-3), -1, 2)
  for(i in 3:(12-3)){
    RW1PrecMat[i+1, ((i):(i+2))]<- c(-1,2,-1)
  }
  strs<- RW1PrecMat

  if(Model == 0){
    npar <- 0
  }else{
    npar <- nstrain
  }

  if(independentChains == 0){
    nstate<- 2^nstrain
  }else{
    nstate<- nstrain
  }
  Bits<- encodeBits(nstrain)

  if(Stan){
    initials <- list(G12 = 0.1, G21 = 0.3, u = rep(0, ndept-1), rraw = rep(0, time), sraw = rep(0, 11), kappa_u=20, kappa_r=20, kappa_s=20, B=rep(0.1, npar), a_k=rep(-10, nstrain))
    initials_list <- lapply(1:nchains, function(x) initials)
    if(GPU){
      mod <- cmdstanr::cmdstan_model(system.file("stan", "multstrain.stan", package = "SpatMet", mustWork = T), compile = F, cpp_options = list(stan_opencl = TRUE))
      if(verbose){
        mod$compile()
        fit<- mod$sample(data = list(ndept=ndept, time=time, nstate=nstate, rankdef=rankdef, nstrain=nstrain, npar=npar, y=yflat, e_it=e_it, R=R,
                                     Bits=Bits, independentChains=independentChains, SMat = strs, Model = Model),
                         init = initials_list, chains = nchains, iter_warmup = round(iter*0.25),
                         iter_sampling = round(iter*0.75), parallel_chains = nchains,
                         seed=seed, adapt_delta = adaptdelta, opencl_ids = c(0, 0))
      }else{
        invisible(capture.output(suppressMessages({
          mod$compile()
          fit<- mod$sample(data = list(ndept=ndept, time=time, nstate=nstate, rankdef=rankdef, nstrain=nstrain, npar=npar, y=yflat, e_it=e_it, R=R,
                                       Bits=Bits, independentChains=independentChains, SMat = strs, Model = Model),
                           init = initials_list, chains = nchains, iter_warmup = round(iter*0.25),
                           iter_sampling = round(iter*0.75), parallel_chains = nchains,
                           seed=seed, adapt_delta = adaptdelta)
        })))
      }
    }else{
      mod <- cmdstanr::cmdstan_model(system.file("stan", "multstrain.stan", package = "SpatMet", mustWork = T), compile = F)
      if(verbose){
        mod$compile()
        fit<- mod$sample(data = list(ndept=ndept, time=time, nstate=nstate, rankdef=rankdef, nstrain=nstrain, npar=npar, y=yflat, e_it=e_it, R=R,
                                     Bits=Bits, independentChains=independentChains, SMat = strs, Model = Model),
                         init = initials_list, chains = nchains, iter_warmup = round(iter*0.25),
                         iter_sampling = round(iter*0.75), parallel_chains = nchains,
                         seed=seed, adapt_delta = adaptdelta)
      }else{
        invisible(capture.output(suppressMessages({
          mod$compile()
          fit<- mod$sample(data = list(ndept=ndept, time=time, nstate=nstate, rankdef=rankdef, nstrain=nstrain, npar=npar, y=yflat, e_it=e_it, R=R,
                                       Bits=Bits, independentChains=independentChains, SMat = strs, Model = Model),
                           init = initials_list, chains = nchains, iter_warmup = round(iter*0.25),
                           iter_sampling = round(iter*0.75), parallel_chains = nchains,
                           seed=seed, adapt_delta = adaptdelta)
        })))
      }
    }
    if(ModEvid){
      ME<- ModelEvidence(y = y, e_it = e_it, Model = Model, adjmat = adjmat, inf.object = fit)
      print(paste0("Marginal loglikelihood is ", ME))
    }
    if(OutbreakProb){
      OutP<- OutbreakProbability(y = y, e_it = e_it, inf.object = fit, adjmat = adjmat, Model = Model, Cyclic = T)
      end_time <- Sys.time()
      time_taken<- end_time - start_time
      print(time_taken)
      return(list(fit, OutP))
    }else{
      end_time <- Sys.time()
      time_taken<- end_time - start_time
      print(time_taken)
      return(fit)
    }
  }else{

    SumYk_vec<- numeric(nstrain)
    for(k in 1:nstrain){
      SumYk_vec[k]<- sum(y[,,k])
    }

    yk<- y[,,1]
    for(k in 2:nstrain){
      yk<- yk + y[,,k]
    }
    crudeResults<- crudeEst(yk, e_it)
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

    zigmaR<- diag(rep(0.1, time), nrow = time, ncol = time)
    zigmaS<- diag(rep(0.1, 11), nrow = 11, ncol = 11)
    zigmaU<- diag(rep(0.08, ndept-1), nrow=ndept-1, ncol=ndept-1)
    optconstantR<- 2.38^2/(time-2)
    optconstantS<- 2.38^2/11
    optconstantU<- 2.38^2/(ndept-1)
    lambdaR<- 1
    lambdaS<- 1
    lambdaU<- 1

    RW2PrecMat<- matrix(0, nrow=time, ncol=time)
    RW2PrecMat[1,(1:3)]<- c(1,-2,1)
    RW2PrecMat[2,(1:4)]<- c(-2,5,-4,1)
    RW2PrecMat[3,(1:5)]<- c(1,-4,6,-4,1)
    RW2PrecMat[(time-1),((time-3):time)]<- c(1,-4,5,-2)
    RW2PrecMat[time,((time-2):time)]<- c(1,-2,1)
    for(i in 3:(time-3)){
      RW2PrecMat[i+1, ((i-1):(i+3))]<- c(1,-4,6,-4,1)
    }
    strr<- RW2PrecMat

    Rsize<- floor(time/8)

    A<- 5+(1:Rsize)
    B<- 5+((Rsize+1):(2*Rsize))
    C<- 5+((2*Rsize+1):(3*Rsize))
    D<- 5+((3*Rsize+1):(4*Rsize))
    E<- 5+((4*Rsize+1):(5*Rsize))
    f<- 5+((5*Rsize+1):(6*Rsize))
    g<- 5+((6*Rsize+1):(7*Rsize))
    H<- 5+((7*Rsize+1):time)

    Blocks<- list(A,B,C,D,E,f,g,H)

    invRconditionalcovA<- solve(RW2PrecMat[A-5, A-5])
    invRconditionalcovB<- solve(RW2PrecMat[B-5, B-5])
    invRconditionalcovC<- solve(RW2PrecMat[C-5, C-5])
    invRconditionalcovD<- solve(RW2PrecMat[D-5, D-5])
    invRconditionalcovE<- solve(RW2PrecMat[E-5, E-5])
    invRconditionalcovf<- solve(RW2PrecMat[f-5, f-5])
    invRconditionalcovg<- solve(RW2PrecMat[g-5, g-5])
    invRconditionalcovH<- solve(RW2PrecMat[H-5, H-5])


    for(i in 2:num_iteration){
      #print(i)

      #proposedkappaR<- rgamma(1, shape = 1 + (time-2)/2, rate = 0.0001 + (t(MC_chain[i-1, 5+(1:time)]) %*% strr %*% MC_chain[i-1, 5+(1:time)])/2)
      #MC_chain[i,3]<- proposedkappaR
      #print(paste("GibbskappaR = ", proposedkappaR))

      proposedkappaS<- rgamma(1, shape = 1 + 11/2, rate = 0.001 + (t(MC_chain[i-1, 5+time+(1:12)]) %*% strs %*% MC_chain[i-1, 5+time+(1:12)])/2)
      MC_chain[i,4]<- proposedkappaS
      #print(paste("GibbskappaS = ", proposedkappaS))

      proposedkappaU<- rgamma(1, shape = 1 + (ndept-1)/2, rate = 0.01 + (t(MC_chain[i-1, 5+time+12+(1:ndept)]) %*% R %*% MC_chain[i-1, 5+time+12+(1:ndept)])/2)
      MC_chain[i,5]<- proposedkappaU
      #print(paste("GibbskappaU = ", proposedkappaU))

      proposedspatcomps<- mvnfast::rmvn(1, mu=MC_chain[i-1, 5+time+12+(1:(ndept-1))], sigma = zigmaU)
      proposedspatcomps<- c(proposedspatcomps, -sum(proposedspatcomps))

      priorcurrentUcomps<- logIGMRF1(MC_chain[i-1, 5+time+12+(1:ndept)], MC_chain[i, 5], R, rankdef)
      priorproposedUcomps<- logIGMRF1(proposedspatcomps, MC_chain[i, 5], R, rankdef)

      proposalproposedcompsU<- mvnfast::dmvn(proposedspatcomps[-ndept], mu = MC_chain[i-1, 5+time+12+(1:(ndept-1))], sigma = zigmaU, log = TRUE)
      proposalcurrentcompsU<- mvnfast::dmvn(MC_chain[i-1, 5+time+12+(1:(ndept-1))], mu = proposedspatcomps[-ndept], sigma = zigmaU, log = TRUE)

      likelihoodcurrent<- multGeneralLoglikelihood_cpp2(y=yflat,ndept=ndept,time=time,nstrain=nstrain,a_k=MC_chain[i-1, 5+time+12+ndept+nstrain+(1:nstrain)], r=MC_chain[i-1,5+(1:time)],s=MC_chain[i-1,5+time+(1:12)],u=MC_chain[i-1,5+time+12+(1:ndept)],Gamma=G(MC_chain[i-1,1],MC_chain[i-1,2]),e_it=e_it, B=MC_chain[i-1, 5+time+12+ndept+(1:nstrain)], Bits=Bits, model=Model,independentChains=independentChains)
      likelihoodproposed<- multGeneralLoglikelihood_cpp2(y=yflat,ndept=ndept,time=time,nstrain=nstrain,a_k=MC_chain[i-1, 5+time+12+ndept+nstrain+(1:nstrain)],r=MC_chain[i-1,5+(1:time)],s=MC_chain[i-1,5+time+(1:12)],u=proposedspatcomps,Gamma=G(MC_chain[i-1,1],MC_chain[i-1,2]),e_it=e_it, B=MC_chain[i-1, 5+time+12+ndept+(1:nstrain)], Bits=Bits, model=Model,independentChains=independentChains)

      mh.ratioU<- exp(likelihoodproposed + priorproposedUcomps + proposalcurrentcompsU
                      - likelihoodcurrent - priorcurrentUcomps - proposalproposedcompsU)

      #print(mh.ratioU)

      if(!is.na(mh.ratioU) && runif(1) < mh.ratioU){
        MC_chain[i,5+time+12+(1:ndept)]<- proposedspatcomps
      }
      else{
        MC_chain[i,5+time+12+(1:ndept)]<- MC_chain[i-1,5+time+12+(1:ndept)]
      }

       for(j in 1:length(Blocks)){
        if(j==1){
          Bj<- length(Blocks[j])
          proposedkappaR<- rgamma(1, shape = 1 + (time-2-Bj)/2, rate = 0.0001 + (t(MC_chain[i-1, unlist(Blocks[-j])]) %*% strr[-(Blocks[[j]]-5),-(Blocks[[j]]-5)] %*% MC_chain[i-1, unlist(Blocks[-j])])/2)
          MC_chain[i,3]<- proposedkappaR

          RW2PrecMat<- MC_chain[i, 3] * strr
          RconditionalcovA<- (1/MC_chain[i, 3])* invRconditionalcovA
          RconditionalcovB<- (1/MC_chain[i, 3])* invRconditionalcovB
          RconditionalcovC<- (1/MC_chain[i, 3])* invRconditionalcovC
          RconditionalcovD<- (1/MC_chain[i, 3])* invRconditionalcovD
          RconditionalcovE<- (1/MC_chain[i, 3])* invRconditionalcovE
          Rconditionalcovf<- (1/MC_chain[i, 3])* invRconditionalcovf
          Rconditionalcovg<- (1/MC_chain[i, 3])* invRconditionalcovg
          RconditionalcovH<- (1/MC_chain[i, 3])* invRconditionalcovH

          covBlocks<- list(RconditionalcovA, RconditionalcovB, RconditionalcovC, RconditionalcovD,
                           RconditionalcovE, Rconditionalcovf, Rconditionalcovg, RconditionalcovH)

          Rconditionalmean<- -covBlocks[[j]] %*% RW2PrecMat[(Blocks[[j]])-5, (unlist(Blocks[-j]))-5] %*% MC_chain[i-1, unlist(Blocks[-j])]
          proposedRcomps<- mvnfast::rmvn(1, mu = Rconditionalmean, sigma = covBlocks[[j]])
          proposedRcomps<- c(proposedRcomps, MC_chain[i-1, unlist(Blocks[-j])])
        }
        else if(j!=1 && j!=length(Blocks)){
          Bj<- length(Blocks[j])
          proposedkappaR<- rgamma(1, shape = 1 + (time-2-Bj)/2, rate = 0.0001 + (t(MC_chain[i-1, unlist(Blocks[-j])]) %*% strr[-(Blocks[[j]]-5),-(Blocks[[j]]-5)] %*% MC_chain[i-1, unlist(Blocks[-j])])/2)
          MC_chain[i,3]<- proposedkappaR

          RW2PrecMat<- MC_chain[i, 3] * strr
          RconditionalcovA<- (1/MC_chain[i, 3])* invRconditionalcovA
          RconditionalcovB<- (1/MC_chain[i, 3])* invRconditionalcovB
          RconditionalcovC<- (1/MC_chain[i, 3])* invRconditionalcovC
          RconditionalcovD<- (1/MC_chain[i, 3])* invRconditionalcovD
          RconditionalcovE<- (1/MC_chain[i, 3])* invRconditionalcovE
          Rconditionalcovf<- (1/MC_chain[i, 3])* invRconditionalcovf
          Rconditionalcovg<- (1/MC_chain[i, 3])* invRconditionalcovg
          RconditionalcovH<- (1/MC_chain[i, 3])* invRconditionalcovH

          covBlocks<- list(RconditionalcovA, RconditionalcovB, RconditionalcovC, RconditionalcovD,
                           RconditionalcovE, Rconditionalcovf, Rconditionalcovg, RconditionalcovH)

          Rconditionalmean<- -covBlocks[[j]] %*% (RW2PrecMat[(Blocks[[j]])-5, (unlist(Blocks[1:(j-1)]))-5] %*% MC_chain[i, unlist(Blocks[1:(j-1)])] + RW2PrecMat[(Blocks[[j]])-5, (unlist(Blocks[(j+1):(length(Blocks))]))-5] %*% MC_chain[i, unlist(Blocks[(j+1):length(Blocks)])])
          proposedRcomps<- mvnfast::rmvn(1, mu = Rconditionalmean, sigma = covBlocks[[j]])
          proposedRcomps<- c(MC_chain[i, unlist(Blocks[1:(j-1)])], proposedRcomps, MC_chain[i, unlist(Blocks[(j+1):length(Blocks)])])
        }
        else if(j==length(Blocks)){
          Bj<- length(Blocks[j])
          proposedkappaR<- rgamma(1, shape = 1 + (time-2-Bj)/2, rate = 0.0001 + (t(MC_chain[i-1, unlist(Blocks[-j])]) %*% strr[-(Blocks[[j]]-5),-(Blocks[[j]]-5)] %*% MC_chain[i-1, unlist(Blocks[-j])])/2)
          MC_chain[i,3]<- proposedkappaR

          RW2PrecMat<- MC_chain[i, 3] * strr
          RconditionalcovA<- (1/MC_chain[i, 3])* invRconditionalcovA
          RconditionalcovB<- (1/MC_chain[i, 3])* invRconditionalcovB
          RconditionalcovC<- (1/MC_chain[i, 3])* invRconditionalcovC
          RconditionalcovD<- (1/MC_chain[i, 3])* invRconditionalcovD
          RconditionalcovE<- (1/MC_chain[i, 3])* invRconditionalcovE
          Rconditionalcovf<- (1/MC_chain[i, 3])* invRconditionalcovf
          Rconditionalcovg<- (1/MC_chain[i, 3])* invRconditionalcovg
          RconditionalcovH<- (1/MC_chain[i, 3])* invRconditionalcovH

          covBlocks<- list(RconditionalcovA, RconditionalcovB, RconditionalcovC, RconditionalcovD,
                           RconditionalcovE, Rconditionalcovf, Rconditionalcovg, RconditionalcovH)

          Rconditionalmean<- -covBlocks[[j]] %*% RW2PrecMat[(Blocks[[j]])-5, (unlist(Blocks[-j]))-5] %*% MC_chain[i, unlist(Blocks[-j])]
          proposedRcomps<- mvnfast::rmvn(1, mu = Rconditionalmean, sigma = covBlocks[[j]])
          proposedRcomps<- c(MC_chain[i, unlist(Blocks[-j])], proposedRcomps)
        }
        likelihoodcurrent<- multGeneralLoglikelihood_cpp2(y=yflat,ndept=ndept,time=time,nstrain=nstrain,a_k=MC_chain[i-1, 5+time+12+ndept+nstrain+(1:nstrain)], r=MC_chain[i-1, 5+(1:time)], s=MC_chain[i-1, 5+time+(1:12)], u=MC_chain[i, 5+time+12+(1:ndept)], Gamma=G(MC_chain[i-1, 1], MC_chain[i-1, 2]), e_it=e_it, B=MC_chain[i-1, 5+time+12+ndept+(1:nstrain)], Bits=Bits, model=Model, independentChains = independentChains)
        likelihoodproposed<- multGeneralLoglikelihood_cpp2(y=yflat,ndept=ndept,time=time,nstrain=nstrain,a_k=MC_chain[i-1, 5+time+12+ndept+nstrain+(1:nstrain)], r=proposedRcomps, s=MC_chain[i-1, 5+time+(1:12)], u=MC_chain[i, 5+time+12+(1:ndept)], Gamma=G(MC_chain[i-1, 1], MC_chain[i-1, 2]), e_it=e_it, B=MC_chain[i-1, 5+time+12+ndept+(1:nstrain)], Bits=Bits, model=Model, independentChains = independentChains)

        mh.ratioR<- exp(likelihoodproposed - likelihoodcurrent)

        #print(paste("mh.ratioR condpriorprop =", mh.ratioR))
        if(!is.na(mh.ratioR) && runif(1) < mh.ratioR){
          MC_chain[i, 5+(1:time)]<- proposedRcomps
        }
        else{
          if(j==1){
            MC_chain[i, 5+(1:time)]<- MC_chain[i-1, 5+(1:time)]
          }
          else if(j!=1){
            MC_chain[i, 5+(1:time)]<- MC_chain[i, 5+(1:time)]
          }
        }
      }

      proposedScomps<- mvnfast::rmvn(1, mu = MC_chain[i-1, 5+time+(1:11)], sigma = zigmaS)
      proposedScomps<- c(proposedScomps, -sum(proposedScomps))

      priorcurrentScomps<- seasonalComp2(MC_chain[i-1, 5+time+(1:12)], MC_chain[i, 4], strs)
      priorproposedScomps<- seasonalComp2(proposedScomps, MC_chain[i, 4], strs)

      proposalproposedScomps<- mvnfast::dmvn(proposedScomps[-12], mu = MC_chain[i-1, 5+time+(1:11)], sigma = zigmaS, log = TRUE)
      proposalcurrentScomps<- mvnfast::dmvn(MC_chain[i-1, 5+time+(1:11)], mu = proposedScomps[-12], sigma = zigmaS, log = TRUE)

      likelihoodcurrent<- multGeneralLoglikelihood_cpp2(y=yflat,ndept=ndept,time=time,nstrain=nstrain,a_k=MC_chain[i-1, 5+time+12+ndept+nstrain+(1:nstrain)],r=MC_chain[i,5+(1:time)],s=MC_chain[i-1,5+time+(1:12)],u=MC_chain[i,5+time+12+(1:ndept)],Gamma=G(MC_chain[i-1,1], MC_chain[i-1,2]),e_it=e_it, B=MC_chain[i-1, 5+time+12+ndept+(1:nstrain)],Bits=Bits, model=Model,independentChains=independentChains)
      likelihoodproposed<- multGeneralLoglikelihood_cpp2(y=yflat,ndept=ndept,time=time,nstrain=nstrain,a_k=MC_chain[i-1, 5+time+12+ndept+nstrain+(1:nstrain)],r=MC_chain[i,5+(1:time)],s=proposedScomps,u=MC_chain[i,5+time+12+(1:ndept)],Gamma=G(MC_chain[i-1,1], MC_chain[i-1,2]),e_it=e_it, B=MC_chain[i-1, 5+time+12+ndept+(1:nstrain)], Bits=Bits,model=Model,independentChains=independentChains)

      mh.ratioS<- exp(likelihoodproposed + priorproposedScomps + proposalcurrentScomps
                      - likelihoodcurrent - priorcurrentScomps - proposalproposedScomps)

      #print(mh.ratioS)

      if(!is.na(mh.ratioS) && runif(1) < mh.ratioS){
        MC_chain[i,5+time+(1:12)]<- proposedScomps
      }
      else{
        MC_chain[i,5+time+(1:12)]<- MC_chain[i-1,5+time+(1:12)]
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
      proposeda_k <- rnorm(nstrain, mean = MC_chain[i-1, 5+time+12+ndept+nstrain+(1:nstrain)], sd = rep(0.07, nstrain))

      likelihoodcurrent<- multGeneralLoglikelihood_cpp2(y=yflat,ndept=ndept,time=time,nstrain=nstrain,a_k=MC_chain[i-1, 5+time+12+ndept+nstrain+(1:nstrain)], r=MC_chain[i, 5+(1:time)], s=MC_chain[i, 5+time+(1:12)], u=MC_chain[i, 5+time+12+(1:ndept)], Gamma=G(MC_chain[i, 1],MC_chain[i,2]), e_it=e_it, B=MC_chain[i, 5+time+12+ndept+(1:nstrain)], Bits=Bits, model=Model, independentChains=independentChains)
      likelihoodproposed<- multGeneralLoglikelihood_cpp2(y=yflat,ndept=ndept,time=time,nstrain=nstrain,a_k=proposeda_k, r=MC_chain[i, 5+(1:time)], s=MC_chain[i, 5+time+(1:12)], u=MC_chain[i, 5+time+12+(1:ndept)], Gamma=G(MC_chain[i, 1], MC_chain[i,2]), e_it=e_it, B=MC_chain[i, 5+time+12+ndept+(1:nstrain)], Bits=Bits, model=Model,independentChains=independentChains)

      mh.ratio<- exp(likelihoodproposed - likelihoodcurrent)

      if(!is.na(mh.ratio) && runif(1) < mh.ratio){
        MC_chain[i, 5+time+12+ndept+nstrain+(1:nstrain)]<- proposeda_k
      }
      else{
        MC_chain[i, 5+time+12+ndept+nstrain+(1:nstrain)]<- MC_chain[i-1, 5+time+12+ndept+nstrain+(1:nstrain)]
      }

      #Gibbs update for a_k's
#      stateDist<- state_dist_cpp2(JointTransitionMatrix_cpp(G(MC_chain[i, 1], MC_chain[i, 2]), K=nstrain))

#      currentR<- MC_chain[i, 5+(1:time)]
#      currentS<- MC_chain[i, 5+time+(1:12)]
#      currentU<- MC_chain[i, 5+time+12+(1:ndept)]
#      currentB<- MC_chain[i, 5+time+12+ndept+(1:nstrain)]
#      poisMean<- 0
#      for(a in 1:ndept){
#        for(b in 1:time){
#          month_index<- (b-1) %% 12 + 1
#          poisMean<- poisMean + e_it[a, b] * exp(currentR[b] + currentS[month_index] + currentU[a] + rep(currentB[1], nstrain)%*%stateDist[c(1,2)] + rep(currentB[2], nstrain)%*%stateDist[c(2,4)])
#        }
#      }
#      proposedAks<- log(rgamma(nstrain, shape = 0.01+SumYk_vec, rate = poisMean + 0.01/exp(-15)))
#      if(Model==0){
#      MC_chain[i, 5+time+12+ndept+nstrain+(1:nstrain)]<- proposedAks
#      }else{
#        proposedcurrentAks<- sum(dgamma(exp(MC_chain[i-1, 5+time+12+ndept+nstrain+(1:nstrain)]), shape = 0.01+SumYk_vec, rate = poisMean + 0.01/exp(-15), log=TRUE) + MC_chain[i-1, 5+time+12+ndept+nstrain+(1:nstrain)])
#        proposedproposedAks<- sum(dgamma(exp(proposedAks), shape = 0.01+SumYk_vec, rate = poisMean + 0.01/exp(-15), log=TRUE) + proposedAks)

#        priorcurrentAks<- sum(dgamma(exp(MC_chain[i-1, 5+time+12+ndept+nstrain+(1:nstrain)]), shape = rep(0.01, nstrain), rate = rep(0.01/exp(-15),nstrain), log=TRUE))
#        priorproposedAks<- sum(dgamma(exp(proposedAks), shape = rep(0.01, nstrain), rate = rep(0.01/exp(-15),nstrain), log=TRUE))

#        likelihoodcurrent<- multGeneralLoglikelihood_cpp2(y=yflat,ndept=ndept,time=time,nstrain=nstrain,a_k=MC_chain[i-1, 5+time+12+ndept+nstrain+(1:nstrain)], r=MC_chain[i,5+(1:time)], s=MC_chain[i,5+time+(1:12)], u=MC_chain[i, 5+time+12+(1:ndept)], Gamma=G(MC_chain[i,1],MC_chain[i,2]),e_it=e_it, B=MC_chain[i, 5+time+12+ndept+(1:nstrain)], Bits=Bits, model=Model, independentChains=independentChains)
#        likelihoodproposed<- multGeneralLoglikelihood_cpp2(y=yflat,ndept=ndept,time=time,nstrain=nstrain,a_k=proposedAks, r=MC_chain[i, 5+(1:time)], s=MC_chain[i, 5+time+(1:12)], u=MC_chain[i, 5+time+12+(1:ndept)], Gamma=G(MC_chain[i,1],MC_chain[i,2]),e_it=e_it, B=MC_chain[i, 5+time+12+ndept+(1:nstrain)], Bits=Bits, model=Model,independentChains=independentChains)

#        mh.ratio<- exp(likelihoodproposed + priorproposedAks + proposedcurrentAks
#                       - likelihoodcurrent - priorcurrentAks - proposedproposedAks)

        #print(mh.ratio)

#        if(!is.na(mh.ratio) && runif(1) < mh.ratio){
#          MC_chain[i, 5+time+12+ndept+nstrain+(1:nstrain)]<- proposedAks
#        }
#        else{
#          MC_chain[i, 5+time+12+ndept+nstrain+(1:nstrain)]<- MC_chain[i-1, 5+time+12+ndept+nstrain+(1:nstrain)]
#        }
#      }

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

        likelihoodcurrent<- multGeneralLoglikelihood_cpp2(y=yflat,ndept=ndept,time=time,nstrain=nstrain,a_k=MC_chain[i, 5+time+12+ndept+nstrain+(1:nstrain)], r=MC_chain[i,5+(1:time)], s=MC_chain[i,5+time+(1:12)], u=MC_chain[i, 5+time+12+(1:ndept)], Gamma=G(MC_chain[i,1],MC_chain[i,2]),e_it=e_it, B=MC_chain[i, 5+time+12+ndept+(1:nstrain)], Bits=Bits, model=Model, independentChains=independentChains)
        likelihoodproposed<- multGeneralLoglikelihood_cpp2(y=yflat,ndept=ndept,time=time,nstrain=nstrain,a_k=MC_chain[i, 5+time+12+ndept+nstrain+(1:nstrain)], r=proposedRcomps, s=MC_chain[i,5+time+(1:12)], u=MC_chain[i, 5+time+12+(1:ndept)], Gamma=G(MC_chain[i,1],MC_chain[i,2]),e_it=e_it, B=MC_chain[i, 5+time+12+ndept+(1:nstrain)], Bits=Bits, model=Model, independentChains=independentChains)

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

      #Adapting zigmaS
      if(i==5){
        epsilonS<- 0.007
        XnS<- MC_chain[1:i, 5+time+(1:11)]
        XnbarS <- colMeans(XnS)
        zigmaS <- cov(XnS) + epsilonS*diag(rep(1, 11))
        zigmaS<- optconstantS * zigmaS
      } else if (i > 5){
        XnbarPrevS <- XnbarS
        XnbarS <- (i*XnbarS + MC_chain[i, 5+time+(1:11)])/(i+1)
        zigmaS <- ((i-1)*zigmaS + tcrossprod(MC_chain[i, 5+time+(1:11)]) + i*tcrossprod(XnbarPrevS) - (i+1)*tcrossprod(XnbarS) + epsilonS*diag(rep(1,11)))/i
        #Robbins Munro tuning
        lambdaS<- lambdaS * exp((2/max(1, i-5)) * (min(mh.ratioS, 1) - 0.234))
        zigmaS<- lambdaS* optconstantS * zigmaS
        #print(zigmaS)
      }

      #Adapting zigmaU
      if(i==5){
        epsilonU<- 0.007
        XnU<- MC_chain[1:i, 5+time+12+(1:(ndept-1))]
        XnbarU <- colMeans(XnU)
        zigmaU <- cov(XnU) + epsilonU*diag(rep(1, ndept-1))
        zigmaU<- optconstantU * zigmaU
      } else if (i > 5){
        XnbarPrevU <- XnbarU
        XnbarU <- (i*XnbarU + MC_chain[i, 5+time+12+(1:(ndept-1))])/(i+1)
        zigmaU <- ((i-1)*zigmaU + tcrossprod(MC_chain[i, 5+time+12+(1:(ndept-1))]) + i*tcrossprod(XnbarPrevU) - (i+1)*tcrossprod(XnbarU) + epsilonU*diag(rep(1,ndept-1)))/i
        #Robbins Munro tuning
        lambdaU<- lambdaU * exp((2/max(1, i-5)) * (min(mh.ratioU, 1) - 0.234))
        zigmaU<- lambdaU* optconstantU * zigmaU
        #print(zigmaU)
      }
      if(i %% 1000 == 0) cat("Iteration:", i, "\n")
    }

    colnames(MC_chain) <- paste(c("G12", "G21", "kappa_r", "kappa_s", "kappa_u", paste("r", 1:time, sep=""), paste("s", 1:12, sep=""), paste("u", 1:ndept, sep=""), paste("B", 1:nstrain, sep=""), paste("a_k", 1:nstrain, sep=""), "StationaryDistribution"))
    MC_chain<- as.data.frame(MC_chain)
    if(ModEvid){
      ME<- ModelEvidence(y = y, e_it = e_it, Model = Model, adjmat = adjmat, inf.object = MC_chain[-(1:2000), ])
      print(paste0("Marginal loglikelihood is ", ME))
    }
    if(OutbreakProb){
      OutP<- OutbreakProbability(y = y, e_it = e_it, inf.object = MC_chain, adjmat = adjmat, Model = Model, Cyclic = T)
      end_time <- Sys.time()
      time_taken<- end_time - start_time
      print(time_taken)
      return(list(MC_chain, OutP))
    }else{
      end_time <- Sys.time()
      time_taken<- end_time - start_time
      print(time_taken)
      return(MC_chain)
    }
  }
}




#MCMC Sampling with C++
CPPmultstrainInfer<- function(y, e_it, nstrain, Model, adjmat, independentChains, num_iteration = 30000, Stan = F, GPU = FALSE, nchains = 4,
                           iter = 4000, seed = NULL, verbose = F, ModEvid = F, OutbreakProb = F, adaptdelta = 0.90, Burn.in = 1000){
  start_time <- Sys.time()

  R<- -1 * adjmat
  diag(R)<- -rowSums(R, na.rm = T)
  rankdef<- nrow(R)-qr(R)$rank
  ndept<- nrow(y[,,1])
  time<- ncol(y[,,1])

  original.y<- y

  #flag missing data for inference
  y<- ifelse(is.na(y), -1, y)

  yflat<- as.numeric(aperm(y, c(2,1,3)))

  RW1PrecMat<- matrix(0, nrow=12, ncol=12)
  RW1PrecMat[1, ]<- c(2,-1, rep(0, 12-3), -1)
  RW1PrecMat[2, ]<- c(-1,2,-1, rep(0, 12-3))
  RW1PrecMat[3, ]<- c(0, -1,2,-1, rep(0, 12-4))
  RW1PrecMat[(12-1), ]<- c(rep(0, 12-3), -1,2,-1)
  RW1PrecMat[12, ]<- c(-1, rep(0, 12-3), -1, 2)
  for(i in 3:(12-3)){
    RW1PrecMat[i+1, ((i):(i+2))]<- c(-1,2,-1)
  }
  strs<- RW1PrecMat

  if(Model == 0){
    npar <- 0
  }else{
    npar <- nstrain
  }

  if(independentChains == 0){
    nstate<- 2^nstrain
  }else{
    nstate<- nstrain
  }
  Bits<- encodeBits(nstrain)

  if(Stan){
    initials <- list(G12 = 0.1, G21 = 0.3, u = rep(0, ndept-1), rraw = rep(0, time), sraw = rep(0, 11), kappa_u=20, kappa_r=20, kappa_s=20, B=rep(0.1, npar), a_k=rep(-10, nstrain))
    initials_list <- lapply(1:nchains, function(x) initials)
    if(GPU){
      mod <- cmdstanr::cmdstan_model(system.file("stan", "multstrain.stan", package = "SpatMet", mustWork = T), compile = F, cpp_options = list(stan_opencl = TRUE))
      if(verbose){
        mod$compile()
        fit<- mod$sample(data = list(ndept=ndept, time=time, nstate=nstate, rankdef=rankdef, nstrain=nstrain, npar=npar, y=yflat, e_it=e_it, R=R,
                                     Bits=Bits, independentChains=independentChains, SMat = strs, Model = Model),
                         init = initials_list, chains = nchains, iter_warmup = round(iter*0.25),
                         iter_sampling = round(iter*0.75), parallel_chains = nchains,
                         seed=seed, adapt_delta = adaptdelta, opencl_ids = c(0, 0))
      }else{
        invisible(capture.output(suppressMessages({
          mod$compile()
          fit<- mod$sample(data = list(ndept=ndept, time=time, nstate=nstate, rankdef=rankdef, nstrain=nstrain, npar=npar, y=yflat, e_it=e_it, R=R,
                                       Bits=Bits, independentChains=independentChains, SMat = strs, Model = Model),
                           init = initials_list, chains = nchains, iter_warmup = round(iter*0.25),
                           iter_sampling = round(iter*0.75), parallel_chains = nchains,
                           seed=seed, adapt_delta = adaptdelta)
        })))
      }
    }else{
      mod <- cmdstanr::cmdstan_model(system.file("stan", "multstrain.stan", package = "SpatMet", mustWork = T), compile = F)
      if(verbose){
        mod$compile()
        fit<- mod$sample(data = list(ndept=ndept, time=time, nstate=nstate, rankdef=rankdef, nstrain=nstrain, npar=npar, y=yflat, e_it=e_it, R=R,
                                     Bits=Bits, independentChains=independentChains, SMat = strs, Model = Model),
                         init = initials_list, chains = nchains, iter_warmup = round(iter*0.25),
                         iter_sampling = round(iter*0.75), parallel_chains = nchains,
                         seed=seed, adapt_delta = adaptdelta)
      }else{
        invisible(capture.output(suppressMessages({
          mod$compile()
          fit<- mod$sample(data = list(ndept=ndept, time=time, nstate=nstate, rankdef=rankdef, nstrain=nstrain, npar=npar, y=yflat, e_it=e_it, R=R,
                                       Bits=Bits, independentChains=independentChains, SMat = strs, Model = Model),
                           init = initials_list, chains = nchains, iter_warmup = round(iter*0.25),
                           iter_sampling = round(iter*0.75), parallel_chains = nchains,
                           seed=seed, adapt_delta = adaptdelta)
        })))
      }
    }
    if(ModEvid){
      ME<- ModelEvidence(y = y, e_it = e_it, Model = Model, adjmat = adjmat, inf.object = fit)
      print(paste0("Marginal loglikelihood is ", ME))
    }
    if(OutbreakProb){
      OutP<- OutbreakProbability(y = y, e_it = e_it, inf.object = fit, adjmat = adjmat, Model = Model, Cyclic = T)
      end_time <- Sys.time()
      time_taken<- end_time - start_time
      print(time_taken)
      return(list(fit, OutP))
    }else{
      end_time <- Sys.time()
      time_taken<- end_time - start_time
      print(time_taken)
      return(fit)
    }
  }else{

    yk<- y[,,1]
    for(k in 2:nstrain){
      yk<- yk + y[,,k]
    }
    crudeResults<- crudeEst(yk, e_it)
    crudeR<- crudeResults[[1]] - mean(crudeResults[[1]])
    crudeS<- crudeResults[[2]]
    crudeU<- crudeResults[[3]]
    crudeU<- ifelse(is.nan(crudeU), mean(crudeU[is.finite(crudeU)]), crudeU)
    crudeblock<- floor(time/12)
    crudeblock<- ((crudeblock*12)-11):(crudeblock*12)
    meanR<- mean(crudeResults[[1]])

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
    for(k in 1:nstrain){
      SumYk_vec[k]<- sum(y[,,k])
    }

    MC_chain<- multInfer_cpp(y=yflat, e_it=e_it, nstrain=nstrain, Model=Model, Bits=Bits, CrudeR=crudeR, CrudeS=crudeS[crudeblock-12], CrudeU=crudeU, RW2PrecMat=RW2PrecMat, RW1PrecMat=RW1PrecMat,
                  R=R, rankdef=rankdef, independentChains=independentChains, num_iteration=num_iteration, meanR = meanR, SumYk_vec = SumYk_vec)

    colnames(MC_chain) <- paste(c("G12", "G21", "kappa_r", "kappa_s", "kappa_u", paste("r", 1:time, sep=""), paste("s", 1:12, sep=""), paste("u", 1:ndept, sep=""), paste("B", 1:nstrain, sep=""), paste("a_k", 1:nstrain, sep=""), "StationaryDistribution"))
    MC_chain<- as.data.frame(MC_chain)
    end_time <- Sys.time()
    time_taken<- end_time - start_time
    print(time_taken)
    return(MC_chain)
  }
}

#chkmcmcMult<- multstrainInfer(y=multmod1[[1]], e_it=multmod1[[2]], nstrain=2, Model=0, adjmat=sim_adjmat, independentChains=0, num_iteration = 1000, Stan = FALSE, GPU = FALSE, nchains = 4,iter = 4000, seed = NULL, verbose = F, ModEvid = F, OutbreakProb = F, adaptdelta = 0.90, Burn.in = 1000)
