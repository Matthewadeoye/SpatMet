#Check valid simulations
check_sim <- function(sim.object, adjmat, Model) {
  y<- sim.object[[1]]
  e_it<- sim.object[[2]]
  r<- sim.object[[3]]
  s<- sim.object[[4]]
  u<- sim.object[[5]]
  x<- sim.object[[6]]
  ndept<- nrow(y)
  time<- ncol(y)

  lambda_it<- matrix(NA, nrow = ndept, ncol = time)
  S_naught<- matrix(NA, nrow = ndept, ncol = time)
  S_one<- matrix(NA, nrow = ndept, ncol = time)

  for(i in 1:ndept){
    for(t in 1:time){
      m<- (t - 1) %% 12 + 1
      lambda_it[i, t] = exp(r[t] + s[m] + u[i])
    }
  }

  if(Model != 0){
    for(i in 1:ndept){
      for(t in 1:time){
        if(x[i, t] == 0){
          S_naught[i, t]<- y[i, t] - e_it[i, t] * lambda_it[i, t]
        }else{
          S_one[i, t]<- y[i, t] - e_it[i, t] * lambda_it[i, t]
        }
      }
    }
    S_naught <- as.numeric(S_naught[!is.na(S_naught)])
    S_one<- as.numeric(S_one[!is.na(S_one)])

    vioSim<- data.frame(value = c(S_naught, S_one), group = factor(c(rep("S0", length(S_naught)), rep("S1", length(S_one)))))
    vioSim$group <- factor(vioSim$group, levels = unique(vioSim$group))
    Means<- c(mean(S_naught), mean(S_one))
    library(RColorBrewer)
    print(ggplot(vioSim, aes(x = group, y = value, fill = group)) +
            geom_violin(trim = FALSE, alpha = 0.7) +
            geom_point(x = 1, y = Means[1], color = "red", size = 3, shape = 18) +
            geom_point(x = 2, y = Means[2], color = "red", size = 3, shape = 18) +
            labs(title = "", x = "", y = "Value", fill = "Mean value") +
            theme_minimal() +
            scale_fill_brewer(palette = "Set2") +
            scale_x_discrete(labels = c("S0" = expression(S[O]),
                                        "S1" = expression(S[1]))) +
            scale_fill_manual(values = brewer.pal(5, "Set2"),
                              labels = c(expression(S[O]),
                                         expression(S[1])),
                              name = "Mean value"))  # Legend title
  }else{
    for(i in 1:ndept){
      for(t in 1:time){
        S_naught[i, t]<- y[i, t] - e_it[i, t] * lambda_it[i, t]
      }
    }
    S_naught <- as.numeric(S_naught[!is.na(S_naught)])
    vioSim<- data.frame(value = S_naught, group = factor(rep("S0", length(S_naught))))
    Means<- mean(S_naught)
    library(RColorBrewer)
    print(ggplot(vioSim, aes(x = group, y = value, fill = group)) +
            geom_violin(trim = FALSE, alpha = 0.7) +
            geom_point(x = 1, y = Means, color = "red", size = 3, shape = 18) +
            labs(title = "", x = "", y = "Value", fill = "Mean value") +
            theme_minimal() +
            scale_fill_brewer(palette = "Set2") +
            scale_x_discrete(labels = c("S0" = expression(S[O]))) +
            scale_fill_manual(values = brewer.pal(5, "Set2"),
                              labels = c(expression(S[O])),
                              name = "Mean value"))  # Legend title
  }
}


#Cholesky confirmation for using mvnfast
check_cholesky <- function(matrix) {
  result <- tryCatch({
    chol(matrix)
    TRUE
  }, error = function(e) {
    FALSE
  })
  return(result)
}

#Prior density for Trend components (r_t)
randomwalk2<- function(componentR, PrecisionR){
  time<- length(componentR)
  Sumres<- 0
  for(i in 3:time){
    res<- (componentR[i-2] - (2 * componentR[i-1]) + componentR[i])^2
    Sumres<- Sumres + res
  }
  return((time - 2)/2 * log(PrecisionR) - PrecisionR/2 * Sumres)
}

#Prior density for Seasonal components (s_t)
seasonalComp<- function(x, z){
  time<- length(x)
  Sumres<- 0
  for(i in 12:time){
    res<- (sum(x[(i-11):(i-0)]))^2
    Sumres<- Sumres + res
  }
  return((time - 11)/2 * log(z) - z/2 * Sumres)
}

#Cyclic RW1 for seasonal components (s_t)
seasonalComp2<- function(x, y, z) {
  n = nrow(z)
  sumC = sum(x[1:(n-1)])
  x = c(x[1:(n-1)], -sumC)
  return ((n - 1)/2 * (log(y) - log(2 * pi)) - y/2 * t(x) %*% z %*% x)
}

#Intrinsic GMRF density for spatial components (u_i)
logIGMRF1<- function(x, y, z, rankdef) {
  n = nrow(z)
  sumC = sum(x[1:(n-1)])
  x = c(x[1:(n-1)], -sumC)
  return ((n - rankdef)/2 * (log(y) - log(2 * pi)) - y/2 * t(x) %*% z %*% x)
}


#' Build transition probability matrix
#'
#' @param G12 Probability of jumping to state 2 if currently at state 1.
#' @param G21 probability of jumping to state 1 if currently in state 2.
#'
#' @return A transition probability matrix
#' @export
#'
#' @examples G(0.1, 0.2)
#'
G<- function(G12, G21){
  m<- matrix(c(1-G12,G12,G21,1-G21),2,2,byrow=T)
  return(m)
}

#Crude estimates
crudeEst<- function(y, e_it){
  ydot<- colSums(y, na.rm = T)
  edot<- colSums(e_it, na.rm = T)
  logydot<- log(ydot)
  logedot<- log(edot)
  lambdadot<- logydot-logedot
  nlambdadot<- lambdadot[lambdadot != -Inf]
  x<- 1:ncol(y)
  lambdadot<- ifelse(lambdadot== -Inf, mean(nlambdadot), lambdadot)
  success <- tryCatch({
    loess_fit <- loess(lambdadot ~ x, span = 0.5)
    TRUE
  }, error = function(e) {
    FALSE
  })

  if(success){
    loess_fit <- loess(lambdadot ~ x, span = 0.5)
    smoothed <- predict(loess_fit)
    crudeS<- lambdadot - smoothed
    crudeR<- smoothed
    crudeU<- log(rowSums(y/e_it, na.rm = T)/sum(exp(crudeR+crudeS)))
    crudeU[crudeU==-Inf]<- mean(crudeU[is.finite(crudeU)])
    crudeU<- crudeU-mean(crudeU[is.finite(crudeU)])
    #crudeU<- log(rowSums(y/e_it, na.rm = T)/sum(exp(crudeR+crudeS)))-mean(log(rowSums(y/e_it, na.rm = T)/sum(exp(crudeR+crudeS))))
    #crudeU<- rep(0, nrow(y))
  }else{
    crudeR<- rep(mean(lambdadot), ncol(y))
    crudeS<- lambdadot-mean(lambdadot)
    crudeU<- rep(0, nrow(y))
  }
  return(list(crudeR, crudeS, crudeU))
}

#Check B for figures
checkB<- function(df){
  success <- tryCatch({
    Bs <- df$draws(variables = "B")
    TRUE
  }, error = function(e) {
    FALSE
  })
  return(success)
}

#Extract posterior credible interval
posterior_interval_custom <- function(posterior_samples, prob = 0.95) {
  # Check if the input is a matrix
  if (!is.matrix(posterior_samples)) {
    stop("Input must be a matrix where rows are samples and columns are variables.")
  }

  lower <- (1 - prob) / 2
  upper <- 1 - lower

  # Apply quantile function across columns
  credible_intervals <- apply(posterior_samples, 2, quantile, probs = c(lower, upper))

  # Convert to dataframe and label rows
  credible_intervals_df <- as.data.frame(t(credible_intervals))
  colnames(credible_intervals_df) <- c("2.5%", "97.5%")

  return(credible_intervals_df)
}


#custom legend
add_legend <- function(...) {
  opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0),
              mar=c(0, 0, 0, 0), new=TRUE)
  on.exit(par(opar))
  plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
  legend(...)
}

#Infer missing data for model evidence computation
InferredData<- function(list.ye, inf.object, adjmat, Model, burn.in = 2000){
  y<- list.ye[[1]]
  e_it<- list.ye[[2]]
  time<- ncol(y)
  ndept<- nrow(y)

  y<- ifelse(is.na(y), -1, y)

  design_matrix_func <- get(paste0("DesignMatrixModel", Model))
  z_it <- design_matrix_func(y, adjmat)[[1]]
  z_it2 <- design_matrix_func(y, adjmat)[[2]]

  if(!is.data.frame(inf.object)){
    fullG12.draws<- stack(as.data.frame(inf.object$draws(variables = "G12")[,1,]))[,1]
    fullG21.draws<- stack(as.data.frame(inf.object$draws(variables = "G21")[,1,]))[,1]
    fullr.draws<- as.data.frame(inf.object$draws(variables = "r")[,1,])
    fulls.draws<- as.data.frame(inf.object$draws(variables = "s")[,1,])
    fullu.draws<- as.data.frame(inf.object$draws(variables = "uconstrained")[,1,])
    fulld.draws<- stack(as.data.frame(inf.object$draws(variables = "state1_stationary_dist")[,1,]))[,1]
  }else{
    fullG12.draws<- as.numeric(inf.object[-(1:burn.in), 1])
    fullG21.draws<- as.numeric(inf.object[-(1:burn.in), 2])
    fullr.draws<- inf.object[-(1:burn.in), 5+(1:time)]
    fulls.draws<- inf.object[-(1:burn.in), 5+time+(1:12)]
    fullu.draws<- inf.object[-(1:burn.in), 5+time+12+(1:ndept)]
    fulld.draws<- inf.object[-(1:burn.in), 5+time+12+ndept+2+1]
  }
  pred_Y<- matrix(NA, nrow = ndept, ncol = time)

  thinning<- numeric(floor(nrow(fullr.draws)/10))
  thinning[1]<- 10
  for(i in 2:length(thinning)){
    thinning[i]<- thinning[i-1] + 10
  }

  G12.draws<- fullG12.draws[thinning]
  G21.draws<- fullG21.draws[thinning]
  r.draws<- fullr.draws[thinning, ]
  s.draws<- fulls.draws[thinning, ]
  u.draws<- fullu.draws[thinning, ]
  d.draws<- fulld.draws[thinning]

  Alldata<- matrix(0, nrow = ndept, ncol = time)

  if(Model == 0){
    for(index in 1:length(thinning)){
      for(i in 1:ndept){
        for(t in 1:time){
          m<- (t - 1) %% 12 + 1
          Exlambda_it <- e_it[i, t] * exp(r.draws[index, t] + s.draws[index, m] + u.draws[index, i])
          pred_Y[i, t]<- rpois(1, Exlambda_it)
        }
      }
      Alldata<- Alldata + pred_Y
    }
  }else if(Model %in% c(1,2,4,5,7)){
    if(!is.data.frame(inf.object)){
      B.draws<- stack(as.data.frame(inf.object$draws(variables = "B")[,1,]))[,1]
    }else{
      B.draws<- as.numeric(inf.object[-(1:burn.in), 5+time+12+ndept+1])
    }
    B.draws<- B.draws[thinning]
    for(index in 1:length(thinning)){
      r<- as.numeric(r.draws[index,])
      s<- as.numeric(s.draws[index,])
      u<- as.numeric(u.draws[index,])
      Ex_Xit <- Decoding(y = y, e_it = e_it, r = r, s = s, u = u, Gamma = G(G12.draws[index], G21.draws[index]), B = B.draws[index], Model = Model, adjmat = adjmat, Cyclic = T)
      for(i in 1:ndept){
        for(t in 1:time){
          m<- (t - 1) %% 12 + 1
          P_Xit<- rbinom(1, 1, prob = Ex_Xit[i, t])
          Exlambda_it <- e_it[i, t] * exp(r.draws[index, t] + s.draws[index, m] + u.draws[index, i] + P_Xit * z_it[i, t] * B.draws[index])
          pred_Y[i, t]<- rpois(1, Exlambda_it)
        }
      }
      Alldata<- Alldata + pred_Y
    }
  }else if(Model %in% c(3,6)){
    if(!is.data.frame(inf.object)){
      B.draws<- as.data.frame(inf.object$draws(variables = "B")[,1,])
    }else{
      B.draws<- inf.object[-(1:burn.in), 5+time+12+ndept+(1:2)]
    }
    B.draws<- B.draws[thinning, ]
    for(index in 1:length(thinning)){
      r<- as.numeric(r.draws[index,])
      s<- as.numeric(s.draws[index,])
      u<- as.numeric(u.draws[index,])
      B<- as.numeric(B.draws[index,])
      Ex_Xit <- Decoding(y = y, e_it = e_it, r = r, s = s, u = u, Gamma = G(G12.draws[index], G21.draws[index]), B = B, Model = Model, adjmat = adjmat, Cyclic = T)
      sum_Xit<- sum_Xit + Ex_Xit
      for(i in 1:ndept){
        for(t in 1:time){
          m<- (t - 1) %% 12 + 1
          P_Xit<- rbinom(1, 1, prob = Ex_Xit[i, t])
          Exlambda_it <- e_it[i, t] * exp(r.draws[index, t] + s.draws[index, m] + u.draws[index, i] + P_Xit * z_it[i, t] * B.draws[index, 1] + P_Xit * z_it2[i, t] * B.draws[index, 2])
          pred_Y[i, t]<- rpois(1, Exlambda_it)
        }
      }
      Alldata<- Alldata + pred_Y
    }
  }
  Alldata<- round(Alldata/length(thinning))
  return(Alldata)
}

#Input missing data for model evidence computation
Inputdata<- function(y, e_it, inf.object, adjmat, Model){
  if(any(is.na(y))){
    ye<- list(y, e_it)
    inf.y<- InferredData(ye,inf.object,adjmat,Model)
    y<- ifelse(is.na(y), inf.y, y)
  }else{
    y<- y
  }
  return(y)
}

#Linear interpolation for application
LinearInterp<- function(popdata){
  ndept<- nrow(popdata)
  time<- ncol(popdata)*12
  e_it<- matrix(NA, nrow = ndept, ncol = time)
  availabledat<- seq(from = 1, to = time, by = 12)
  e_it[ , availabledat]<- popdata[ , 1:ncol(popdata)]

  for(i in 1:ndept){
    for(t in 1:(ncol(popdata)-1)){
      ind<- availabledat[t]
      ind2<- availabledat[t+1]
      df<- data.frame(x = c(ind,ind2), y = e_it[i, c(ind, ind2)])
      interv<- seq(from = ind+1, to = ind2-1, by = 1)
      e_it[i, (ind+1):(ind2-1)]<- approx(df$x, df$y, xout = interv)$y
    }
  }
  e_it<- e_it[ , -(max(availabledat):time)]
  return(round(e_it))
}


#Datasets for application (Meningococcal)
#Susceptibles<- read.csv("popn.csv")
#countries<- Susceptibles[ , 1]
#eit<- Susceptibles
#names(eit)<- NULL
#eit<- eit[ , -1]
#eit<- as.matrix(eit)
#eit<- LinearInterp(eit)
#dim(eit)
#Infected<- read.csv("Cases.csv")
#y<- matrix(Infected[ ,3], nrow = 28, ncol = 252, byrow = TRUE)
#alldata<- list(y, eit, countries)
#sim.plot(alldata)

#AdjacencyMatrix for application
#poly <- cshapes::cshp(date=as.Date("2019-12-31"), useGW=TRUE)
#Allcountriesnamescodes<- data.frame(countryname=poly$country_name, countrycode=poly$gwcode)
#requiredcountriesnames<-c("Austria","Belgium","Cyprus","Czech Republic",
#                         "Denmark","Estonia","Finland","France","German Federal Republic",
#                         "Greece","Hungary","Iceland","Ireland","Italy/Sardinia","Latvia",
#                         "Lithuania","Luxembourg","Malta","Netherlands","Norway","Poland",
#                         "Portugal","Rumania","Slovakia","Slovenia","Spain","Sweden","United Kingdom")

#1000km for 30%, 820km for 20%
#dmat <- cshapes::distmatrix(as.Date("2019-12-31"), type="capdist")
#colnames(dmat)<- Allcountriesnamescodes$countryname
#rownames(dmat)<- Allcountriesnamescodes$countryname
#dmat<- dmat[requiredcountriesnames, requiredcountriesnames]
#AdjacencyMatrix<- ifelse(dmat>820,0,1)
#diag(AdjacencyMatrix)<- 0

ModComp<- function(allmods, alldata, adjmat){
  y<- alldata[[1]]
  e_it<- alldata[[2]]
  modevid<- matrix(NA, nrow = 1, ncol=length(allmods))
  postmodprob<- matrix(NA, nrow = 1, ncol=length(allmods))
  for(i in 1:length(allmods)){
    Model<- allmods[[i]]
    modevid[,i]<- DetectOutbreaks::ModelEvidence(y=y, e_it = e_it, adjmat = adjmat, Model = i-1, inf.object = Model)
  }
  for(i in 1:length(allmods)){
    postmodprob[,i]<- exp(modevid[,i]- matrixStats::logSumExp(c(modevid)))
  }
  colnames(modevid)<- paste("Model", 0:(length(allmods)-1))
  colnames(postmodprob)<- paste("Model", 0:(length(allmods)-1))
  print(modevid)
  print(postmodprob)
}

secondorderRW<- function(rt, kappa_r, forecastlength){
  rt_forecast<- numeric(forecastlength)
  rt_forecast[1] <- 2 * rt[length(rt)-1] - rt[length(rt)-2] + rnorm(1, mean=0, sd=sqrt(1/kappa_r))
  rt_forecast[2] <- 2 * rt_forecast[1] - rt[length(rt)-1] + rnorm(1, mean=0, sd=sqrt(1/kappa_r))
  for(t in 3:forecastlength){
    rt_forecast[t] <- 2 * rt_forecast[t-1] - rt_forecast[t-2] + rnorm(1, mean=0, sd=sqrt(1/kappa_r))
  }
  return(rt_forecast)
}

MChain<- function(forecastlength, G12, G21, outP){
  transition_matrix <- matrix(c(1-G12, G12, G21, 1-G21), nrow = 2, byrow = TRUE)
  initial_state <- sample(0:1, size = 1, prob = outP)
  states <- numeric(forecastlength)
  states[1] <- initial_state
  for(i in 2:forecastlength) {
    current_state <- states[i - 1]
    next_state <- sample(0:1, size = 1, prob = transition_matrix[current_state + 1, ])
    states[i] <- next_state
  }
  return(states)
}


LinearExtrp<- function(popdata, forecastlength){
  ndept<- nrow(popdata)
  time<- ncol(popdata)
  e_it<- matrix(NA, nrow = ndept, ncol = forecastlength)

  for(i in 1:ndept){
    e_it[i, 1]<- popdata[i, time] + ((popdata[i, time]-popdata[i, time-1])/1) * 1
    e_it[i, 2]<- e_it[i, 1] + ((e_it[i, 1]-popdata[i, time])/1) * 1
    for(t in 3:forecastlength){
      e_it[i, t]<- e_it[i, t-1] + ((e_it[i, t-1]-e_it[i, t-2])/1) * 1
    }
  }
  return(round(e_it))
}


#FORWARD FILTER for new (cyclic) model
forwardfilter2<- function(y, e_it, r, s, u, Gamma, B, Model, adjmat) {

  design_matrix_func <- get(paste0("DesignMatrixModel", Model))
  z_it <- design_matrix_func(y, adjmat)[[1]]
  z_it2 <- design_matrix_func(y, adjmat)[[2]]

  ndept <- nrow(y)
  nstate <- ncol(Gamma)
  time <- ncol(y)
  gamma_11 <- log(Gamma[1, 1])
  gamma_12 <- log(Gamma[1, 2])
  gamma_21 <- log(Gamma[2, 1])
  gamma_22 <- log(Gamma[2, 2])
  init_density<- state_dist_cpp(Gamma[1, 2], Gamma[2, 1])
  init_density<- log(init_density)

  y<- ifelse(is.na(y), -1, y)

  if(Model %in% c(0,1,2,4,5,7)){
    AllForwardprobs<- vector("list", ndept)

    for (i in 1:ndept) {
      Forwardprob<- matrix(NA, nrow = time, ncol = nstate)
      if(y[i, 1] == -1){
        alpha.1 <- init_density[1]
        alpha.2 <- init_density[2]
      }else{
        alpha.1 <- init_density[1] + dpois(y[i, 1], lambda = e_it[i, 1] * exp(r[1] + s[1] + u[i]), log = TRUE)
        alpha.2 <- init_density[2] + dpois(y[i, 1], lambda = e_it[i, 1] * exp(r[1] + s[1] + u[i] + B[1] * z_it[i, 1]), log = TRUE)
      }
      Forwardprob[1, ] <- c(alpha.1, alpha.2)

      for (t in 2:time) {
        if(y[i, t] == -1){
          Alphas<- c(alpha.1, alpha.2)
          alpha.1 <- logSumExp_cpp(c(Alphas[1] + gamma_11))
          alpha.2 <- logSumExp_cpp(c(Alphas[1] + gamma_12))
          Forwardprob[t, ] <- c(alpha.1, alpha.2)
        }else{
          month_index<- (t - 1) %% 12 + 1
          Alphas<- c(alpha.1, alpha.2)
          alpha.1 <- logSumExp_cpp(c(Alphas[1] + gamma_11 + dpois(y[i, t], lambda = e_it[i, t] * exp(r[t] + s[month_index] + u[i]), log = TRUE), Alphas[2] + gamma_21 + dpois(y[i, t], lambda = e_it[i, t] * exp(r[t] + s[month_index] + u[i]), log = TRUE)))
          alpha.2 <- logSumExp_cpp(c(Alphas[1] + gamma_12 + dpois(y[i, t], lambda = e_it[i, t] * exp(r[t] + s[month_index] + u[i] + B[1] * z_it[i,t]), log = TRUE), Alphas[2] + gamma_22 + dpois(y[i, t], lambda = e_it[i, t] * exp(r[t] + s[month_index] + u[i] + B[1] * z_it[i,t]), log = TRUE)))
          Forwardprob[t, ] <- c(alpha.1, alpha.2)
        }
      }
      AllForwardprobs[[i]]<- Forwardprob
    }

    return(AllForwardprobs)
  }

  else if(Model %in% c(3,6)){
    AllForwardprobs<- vector("list", ndept)

    for (i in 1:ndept) {
      Forwardprob<- matrix(NA, nrow = time, ncol = nstate)
      if(y[i, 1] == -1){
        alpha.1 <- init_density[1]
        alpha.2 <- init_density[2]
      }else{
        alpha.1 <- init_density[1] + dpois(y[i, 1], lambda = e_it[i, 1] * exp(r[1] + s[1] + u[i]), log = TRUE)
        alpha.2 <- init_density[2] + dpois(y[i, 1], lambda = e_it[i, 1] * exp(r[1] + s[1] + u[i] + B[1] * z_it[i,1] + B[2] * z_it2[i,1]), log = TRUE)
      }

      Forwardprob[1, ] <- c(alpha.1, alpha.2)

      for (t in 2:time) {
        if(y[i, t] == -1){
          Alphas<- c(alpha.1, alpha.2)
          alpha.1 <- logSumExp_cpp(c(Alphas[1] + gamma_11))
          alpha.2 <- logSumExp_cpp(c(Alphas[1] + gamma_12))
          Forwardprob[t, ] <- c(alpha.1, alpha.2)
        }else{
          month_index<- (t - 1) %% 12 + 1
          Alphas<- c(alpha.1, alpha.2)
          alpha.1 <- logSumExp_cpp(c(Alphas[1] + gamma_11 + dpois(y[i, t], lambda = e_it[i, t] * exp(r[t] + s[month_index] + u[i]), log = TRUE), Alphas[2] + gamma_21 + dpois(y[i, t], lambda = e_it[i, t] * exp(r[t] + s[month_index] + u[i]), log = TRUE)))
          alpha.2 <- logSumExp_cpp(c(Alphas[1] + gamma_12 + dpois(y[i, t], lambda = e_it[i, t] * exp(r[t] + s[month_index] + u[i] + B[1] * z_it[i,t] + B[2] * z_it2[i,t]), log = TRUE), Alphas[2] + gamma_22 + dpois(y[i, t], lambda = e_it[i, t] * exp(r[t] + s[month_index] + u[i] + B[1] * z_it[i,t] + B[2] * z_it2[i,t]), log = TRUE)))
          Forwardprob[t, ] <- c(alpha.1, alpha.2)
        }
      }
      AllForwardprobs[[i]]<- Forwardprob
    }

    return(AllForwardprobs)
  }
}

state_prediction<- function(forecastlength, y, e_it, r, s, u, Gamma, B, Model, adjmat){
  ndept<- nrow(y)
  time<- ncol(y)
  unmlfp<- forwardfilter2(y, e_it, r, s, u, Gamma, B, Model, adjmat)
  Allstatepreds<- matrix(NA, nrow = ndept, ncol = forecastlength)
  for(i in 1:ndept){
    la<- t(unmlfp[[i]])
    Maxla <- max(la[,time])
    llk <- Maxla + log(sum(exp(la[,time] - Maxla)))
    statepreds <- numeric(forecastlength)
    nmlfp<- exp(la[,time] - llk)
    for(f in 1:forecastlength){
      nmlfp <- nmlfp %*% Gamma
      statepreds[f] <- nmlfp[2]
    }
    Allstatepreds[i, ]<- statepreds
  }
  return(Allstatepreds)
}

logI<- function(theta, Model, speed, y, e_it, z_it, z_it2, SMat, R, rankdef, mu, varcov){
  time<- ncol(y)
  ndept<- nrow(y)

  if(Model==0){
    if(theta[1] <= 0 || theta[2] <= 0 || theta[3]<= 0){
      estimate<- NaN
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
      log_lik<- GeneralLoglikelihood_cpp2(y=y, r=r, s=s, u=u, Gamma=G_mat, e_it=e_it, B = ARcoeff, model = Model, z_it=z_it, z_it2=z_it2)
      log_prior<- sum(dgamma(c(kappaR, kappaS), shape = c(1, 1), rate = c(0.0001, 0.001), log = TRUE)) +
        dgamma(kappaU, shape = 1, rate = 0.01, log = TRUE) +
        randomwalk2(r, kappaR) +
        seasonalComp2(s, kappaS, SMat) +
        logIGMRF1(u, kappaU, R, rankdef)
      if(speed==1){
        log_proposal<- mvnfast::dmvt(theta, mu = mu, sigma = varcov, df = 3, log = TRUE)
      }else{
        log_proposal<- mvtnorm::dmvt(theta, delta = mu, sigma = varcov, df = 3, log = TRUE)
      }
      #print(paste("log_lik is ", log_lik))
      #print(paste("log_prior is ", log_prior))
      #print(paste("log_proposal is ", log_proposal))
      estimate<- log_lik + log_prior - log_proposal
    }
  }else if(Model %in% c(1,2,3,4,5,6,7)){
    if(theta[1] <= 0 || theta[1] >= 1 ||theta[2] <= 0 || theta[2] >= 1 || theta[3] <= 0 || theta[4] <= 0 || theta[5]<= 0){
      estimate<- NaN
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
      log_lik<- GeneralLoglikelihood_cpp2(y=y, r=r, s=s, u=u, Gamma=G_mat, e_it=e_it, B = ARcoeff, model = Model, z_it=z_it, z_it2=z_it2)
      log_prior<- (sum(dbeta(Gammas, shape1 = c(2, 2), shape2 = c(2, 2), log = TRUE)) +
        sum(dgamma(c(kappaR, kappaS), shape = c(1, 1), rate = c(0.0001, 0.001), log = TRUE)) +
        dgamma(kappaU, shape = 1, rate = 0.01, log = TRUE) +
        randomwalk2(r, kappaR) +
        seasonalComp2(s, kappaS, SMat) +
        logIGMRF1(u, kappaU, R, rankdef) +
        sum(dgamma(ARcoeff, shape = rep(2, length(ARcoeff)), rate = rep(2, length(ARcoeff)), log = TRUE)))
      if(speed==1){
        log_proposal<- mvnfast::dmvt(theta, mu = mu, sigma = varcov, df = 3, log = TRUE)
      }else{
        log_proposal<- mvtnorm::dmvt(theta, delta = mu, sigma = varcov, df = 3, log = TRUE)
      }
      estimate<- log_lik + log_prior - log_proposal
    }
  }
  return(estimate)
}

logI2<- function(theta, Model, speed, y, e_it, z_it, z_it2, SMat, R, rankdef, mu, varcov){
  time<- ncol(y)
  ndept<- nrow(y)

  if(Model==0){
    if(theta[1] <= 0 || theta[2] <= 0 || theta[3]<= 0){
      estimate<- -Inf
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
      log_lik<- GeneralLoglikelihood_cpp2(y=y, r=r, s=s, u=u, Gamma=G_mat, e_it=e_it, B = ARcoeff, model = Model, z_it=z_it, z_it2=z_it2)
      log_prior<- sum(dgamma(c(kappaR, kappaS), shape = c(1, 1), rate = c(0.0001, 0.001), log = TRUE)) +
        dgamma(kappaU, shape = 1, rate = 0.01, log = TRUE) +
        randomwalk2(r, kappaR) +
        seasonalComp2(s, kappaS, SMat) +
        logIGMRF1(u, kappaU, R, rankdef)
      if(speed==1){
        log_proposal<- mvnfast::dmvn(theta, mu = mu, sigma = varcov, log = T)
      }else{
        log_proposal<- mvtnorm::dmvn(theta, mean = mu, sigma = varcov, log = T)
      }
      print(paste("log_lik is ", log_lik))
      print(paste("log_prior is ", log_prior))
      print(paste("log_proposal is ", log_proposal))
      estimate<- log_lik + log_prior - log_proposal
    }
  }else if(Model %in% c(1,2,3,4,5,6,7)){
    if(theta[1] <= 0 || theta[1] >= 1 ||theta[2] <= 0 || theta[2] >= 1 || theta[3] <= 0 || theta[4] <= 0 || theta[5]<= 0){
      estimate<- -Inf
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
      log_lik<- GeneralLoglikelihood_cpp2(y=y, r=r, s=s, u=u, Gamma=G_mat, e_it=e_it, B = ARcoeff, model = Model, z_it=z_it, z_it2=z_it2)
      log_prior<- (sum(dbeta(Gammas, shape1 = c(2, 2), shape2 = c(2, 2), log = TRUE)) +
                     sum(dgamma(c(kappaR, kappaS), shape = c(1, 1), rate = c(0.0001, 0.001), log = TRUE)) +
                     dgamma(kappaU, shape = 1, rate = 0.01, log = TRUE) +
                     randomwalk2(r, kappaR) +
                     seasonalComp2(s, kappaS, SMat) +
                     logIGMRF1(u, kappaU, R, rankdef) +
                     sum(dgamma(ARcoeff, shape = rep(2, length(ARcoeff)), rate = rep(2, length(ARcoeff)), log = TRUE)))
      if(speed==1){
        log_proposal<- mvnfast::dmvn(theta, mu = mu, sigma = varcov, log = TRUE)
      }else{
        log_proposal<- mvtnorm::dmvn(theta, mean = mu, sigma = varcov, log = TRUE)
      }
      estimate<- log_lik + log_prior - log_proposal
    }
  }
  return(estimate)
}

ESS<- function(log_weights) {
  normalized_log_weights<- log_weights - logSumExp_cpp(log_weights)
  ess<- exp(-logSumExp_cpp(2 * normalized_log_weights))
  return(ess)
}


LSEXPs<- function(logdensities, Paper.Method=TRUE){
  logDensProposal<- logdensities[[1]]
  logDensPosterior<- logdensities[[2]]

  numerator<- numeric(length(logDensProposal))
  denominator<- numeric(length(logDensPosterior))

  s1<- length(logDensPosterior)/(length(logDensProposal)+length(logDensPosterior))
  s2<- length(logDensProposal)/(length(logDensProposal)+length(logDensPosterior))

  if(Paper.Method){
    rt<- 0.00001
  medL<- median(logDensPosterior)
  delta<- 0.3

  while(delta > 0.01){
    for(i in 1:length(logDensProposal)){
      numerator[i]<- exp(logDensProposal[i] - medL)/(s1 * exp(logDensProposal[i] - medL) + s2 * rt)
    }
    for(j in 1:length(logDensPosterior)){
      denominator[j]<- 1/(s1 * exp(logDensPosterior[j] - medL) + s2 * rt)
    }
    numerator<- numerator[numerator != -Inf]
    current.rt<- (sum(numerator)/length(numerator))/(sum(denominator)/length(logDensPosterior))
    delta<- abs(rt-current.rt)
    print(current.rt)
    rt<- current.rt
  }
  current.rt<- log(current.rt) + medL
  return(current.rt)
  }else{
    logPy<- -3305.074
    delta<- 0.03

    while(delta > 0.01){
      for(i in 1:length(logDensProposal)){
        numerator[i]<- logDensProposal[i] - logSumExp_cpp(c(log(s1)+logDensProposal[i], log(s2) + logPy))
      }
      for(j in 1:length(logDensPosterior)){
        denominator[j]<- log(1) - logSumExp_cpp(c(log(s1)+logDensPosterior[j], log(s2) + logPy))
      }
      numerator<- numerator[numerator != -Inf]
      denominator<- denominator[denominator != -Inf]
      currentlogPy<- (-log(length(numerator))+logSumExp_cpp(numerator))-(-log(length(logDensPosterior))+logSumExp_cpp(denominator))
      delta<- abs(logPy-currentlogPy)
      print(currentlogPy)
      logPy<- currentlogPy
    }
    return(currentlogPy)
  }
}

unique_rows<- function(mat, precision=10) {
  mat_rounded<- round(mat, precision)
  row_strings<- apply(mat_rounded, 1, paste, collapse = ",")
  unique_count<- length(unique(row_strings))
  return(unique_count)
}

#Assuming the same TPM for all strains
JointTransitionMatrix<- function(gamma, K) {
  S<- 2^K
  Gamma<- matrix(0, nrow = S, ncol = S)
  for(a in 0:(S - 1)) {
    for(b in 0:(S - 1)) {
      prob<- 1
      for(k in 1:K) {
        from_k<- (a %/% 2^(k - 1)) %% 2
        to_k<- (b %/% 2^(k - 1)) %% 2
        prob<- prob * gamma[from_k + 1, to_k + 1]
      }
      Gamma[a + 1, b + 1]<- prob
    }
  }
  return(Gamma)
}

#Assuming each strain has its unique TPM
JointTransitionMatrix_per_strain<- function(gamma_list){
  K<- length(gamma_list)
  S<- 2^K
  Gamma<- matrix(0, nrow = S, ncol = S)
  for(a in 0:(S - 1)) {
    for(b in 0:(S - 1)) {
      prob<- 1
      for(k in 1:K) {
        from_k<- (a %/% 2^(k - 1)) %% 2
        to_k<- (b %/% 2^(k - 1)) %% 2
        gamma_k<- gamma_list[[k]]
        prob<- prob * gamma_k[from_k + 1, to_k + 1]
      }
      Gamma[a + 1, b + 1]<- prob
    }
  }
  return(Gamma)
}

gaussian_copula_cdf <- function(u,  corrMat) {
  eps <- 1e-12
  q <- qnorm(pmin(pmax(u, eps), 1 - eps))
  mvtnorm::pmvnorm(lower = rep(-Inf, length(u)), upper = q, sigma =  corrMat)
}

.onLoad <- function(lib, pkg){
  ns <- asNamespace(pkg)
  Memoisegaussian_copula_cdf <- memoise::memoise(gaussian_copula_cdf_cpp)
  assign("Memoisegaussian_copula_cdf", Memoisegaussian_copula_cdf, envir = ns)
}

#Memoisegaussian_copula_cdf<- memoise::memoise(gaussian_copula_cdf_cpp)
#assign("Memoisegaussian_copula_cdf", Memoisegaussian_copula_cdf, envir = globalenv())

#Dependence modelling with copula, assuming the same TPM for all strains
JointTransitionMatrix_copula<- function(gamma, K, copulaParams){

  #cop<- copula::normalCopula(param = copulaParam, dim = K, dispstr = "un")
  #cop<- copula::plackettCopula(param = copulaParam)

  gaussian_copula_cdf_cached <- memoise::memoise(gaussian_copula_cdf)

  S <- 2^K
  corrMat <- diag(K)
  uppertriang <- copulaParams
  gdata::upperTriangle(corrMat, byrow=TRUE) <- uppertriang
  gdata::lowerTriangle(corrMat, byrow=FALSE) <- uppertriang

  gamma[1, ] <- gamma[1, ][2:1]

  Gamma <- matrix(0, S, S)

  for(a in 0:(S-1)) {
    for(b in 0:(S-1)) {
      Indices <- integer(0)
      IndicesComplement <- integer(0)
      prob <- numeric(K)

      for(k in 1:K){
        from_k <- (a %/% 2^(k-1)) %% 2
        to_k   <- (b %/% 2^(k-1)) %% 2

        if(from_k == 1) Indices <- c(Indices,k)
        else IndicesComplement <- c(IndicesComplement,k)

        prob[k] <- gamma[from_k+1, to_k+1]
      }

      subsets <- sets::set_power(sets::as.set(IndicesComplement))
      subsets <- lapply(subsets, function(g) unlist(as.vector(g)))

      total <- 0
      for(Tset in subsets){
        sign <- (-1)^length(Tset)
        idx <- c(Indices, Tset)
        u <- rep(1, K)
        if(length(idx)>0) u[idx] <- prob[idx]
        #total <- total + sign * copula::pCopula(u, cop)
        total <- total + sign * gaussian_copula_cdf_cached(u, corrMat)
      }
      Gamma[a+1, b+1] <- total
    }
  }
  Gamma <- Gamma / rowSums(Gamma)
  Gamma
}

#Dependence modelling with copula, assuming the same TPM for all strains
VectorizedJointTransitionMatrix_copula<- function(gamma, K, copulaParams){
  gaussian_copula_cdf_cached <- memoise::memoise(gaussian_copula_cdf_cpp)
  S <- 2^K

  states <- matrix(0, nrow=S, ncol=K)
  for(k in 1:K)
    states[,k] <- (0:(S-1) %/% 2^(k-1)) %% 2

  corrMat <- diag(K)
  uppertriang <- copulaParams
  gdata::upperTriangle(corrMat, byrow=TRUE) <- uppertriang
  gdata::lowerTriangle(corrMat, byrow=FALSE) <- uppertriang

  gamma[1, ] <- gamma[1, ][2:1]
  Gamma_rows <- matrix(0, nrow=S, ncol=S)
  for (a in 0:(S - 1)) {
    from_vec <- states[a+1, ]
    row_vals <- numeric(S)
    for (b in 0:(S - 1)) {
      to_vec <- states[b+1, ]
      prob   <- gamma[cbind(from_vec+1, to_vec+1)]

      Indices<- which(from_vec == 1)
      IndicesComplement<- which(from_vec == 0)
      m <- length(IndicesComplement)
      M <- 2^m

      total <- 0
      # Pre-allocate U matrix
      U <- matrix(1, nrow=M, ncol=K)
      U[, Indices] <- rep(prob[Indices], each=M)
      for (i in 0:(M-1)) {
        bits <- as.logical(intToBits(i)[1:m])
        idx  <- IndicesComplement[bits]
        if (length(idx) > 0)
          U[i+1, idx] <- prob[idx]
        sign <- (-1) ^ sum(bits)
        total <- total + sign * gaussian_copula_cdf_cached(U[i+1, ], corrMat)
      }
      row_vals[b+1] <- total
    }
    Gamma_rows[a+1, ] <- row_vals
  }
  Gamma <- Gamma_rows
  Gamma <- Gamma / rowSums(Gamma)
  return(Gamma)
}

#Dependence modelling with copula, assuming different TPM for all strains
VectorizedJointTransitionMatrix_copula_per_strain<- function(gamma_list, K, copulaParams){
  gaussian_copula_cdf_cached <- memoise::memoise(gaussian_copula_cdf_cpp)
  S <- 2^K

  states <- matrix(0, nrow = S, ncol = K)
  for (k in 1:K)
    states[, k] <- (0:(S - 1) %/% 2^(k - 1)) %% 2

  corrMat <- diag(K)
  gdata::upperTriangle(corrMat, byrow = TRUE)  <- copulaParams
  gdata::lowerTriangle(corrMat, byrow = FALSE) <- copulaParams

  Gamma_rows <- matrix(0, nrow = S, ncol = S)
  for (a in 0:(S - 1)) {
    from_vec <- states[a + 1, ]
    row_vals <- numeric(S)
    for (b in 0:(S - 1)) {
      to_vec <- states[b + 1, ]
      prob <- numeric(K)
      for (k in 1:K) {
        gamma_k<- gamma_list[[k]]
        gamma_k[1, ] <- gamma_k[1, ][2:1]
        prob[k] <- gamma_k[from_vec[k] + 1, to_vec[k] + 1]
      }
      Indices<- which(from_vec == 1)
      IndicesComplement<- which(from_vec == 0)
      m <- length(IndicesComplement)
      M <- 2^m

      total <- 0
      U <- matrix(1, nrow = M, ncol = K)
      if (length(Indices) > 0)
        U[, Indices] <- rep(prob[Indices], each = M)
      for (i in 0:(M - 1)) {
        bits <- as.logical(intToBits(i)[1:m])
        idx  <- IndicesComplement[bits]
        if (length(idx) > 0)
          U[i + 1, idx] <- prob[idx]
        sign  <- (-1)^sum(bits)
        total <- total + sign * gaussian_copula_cdf_cached(U[i + 1, ], corrMat)
      }
      row_vals[b + 1] <- total
    }
    Gamma_rows[a + 1, ] <- row_vals
  }
  Gamma <- Gamma_rows
  Gamma <- Gamma / rowSums(Gamma)
  return(Gamma)
}

Multipurpose_JointTransitionMatrix<- function(gammas, K, Lambdas, Modeltype, gh){
  nstate<- 2^K
  if(Modeltype==1){
    JointTPM<- JointTransitionMatrix_arma_cpp2(G(gammas[1], gammas[2]), K)
  }else if(Modeltype==2){
    Glist<- BuildGamma_list_cpp(gammas)
    JointTPM<- JointTransitionMatrix_per_strain_cpp2(Glist, K)
  }else if(Modeltype==3){
    JointTPM<- ParallelGaussianJointTransitionMatrix_copula_cpp(G(gammas[1], gammas[2]), K, Lambdas, gh$nodes, gh$weights)
  }else if(Modeltype==4){
    Glist<- BuildGamma_list_cpp(gammas)
    JointTPM<- ParallelGaussianJointTransitionMatrix_copula_per_strain_cpp(Glist, K, Lambdas, gh$nodes, gh$weights)
  }
  return(JointTPM)
}


one_factor_copula_cdf <- function(u, lambda) {
  eps <- 1e-12
  q <- qnorm(pmin(pmax(u, eps), 1 - eps))

  integrand <- function(z) {
    out <- numeric(length(z))
    for (i in seq_along(z)) {
      out[i] <- prod(
        pnorm((q - lambda * z[i]) / sqrt(1 - lambda^2))
      )
    }
    out * dnorm(z)
  }

  integrate(integrand, -Inf, Inf, subdivisions = 200L)$value
}

one_factor_copula_cdf2<- function(u, lambda) {
  gh <- statmod::gauss.quad(30, kind = "hermite")
  one_factor_copula_cdf_rcpp(u, lambda, gh$nodes, gh$weights)
}

#Dependence modelling with copula, assuming the same TPM for all strains
FactorJointTransitionMatrix_copula<- function(gamma, K, copulaParams){

  S <- 2^K
  gamma[1, ] <- gamma[1, ][2:1]

  Gamma <- matrix(0, S, S)

  for(a in 0:(S-1)) {
    for(b in 0:(S-1)) {
      Indices <- integer(0)
      IndicesComplement <- integer(0)
      prob <- numeric(K)

      for(k in 1:K){
        from_k <- (a %/% 2^(k-1)) %% 2
        to_k   <- (b %/% 2^(k-1)) %% 2

        if(from_k == 1) Indices <- c(Indices,k)
        else IndicesComplement <- c(IndicesComplement,k)

        prob[k] <- gamma[from_k+1, to_k+1]
      }

      subsets <- sets::set_power(sets::as.set(IndicesComplement))
      subsets <- lapply(subsets, function(g) unlist(as.vector(g)))

      total <- 0
      for(Tset in subsets){
        sign <- (-1)^length(Tset)
        idx <- c(Indices, Tset)
        u <- rep(1, K)
        if(length(idx)>0) u[idx] <- prob[idx]
        total <- total + sign * one_factor_copula_cdf2(u, copulaParams)
      }
      Gamma[a+1, b+1] <- total
    }
  }
  Gamma <- Gamma / rowSums(Gamma)
  Gamma
}


frank_cdf <- function(u, theta) {
  if (theta == 0) {
    cdf<- prod(u)
  }else{
    d <- length(u)
    num <- prod(exp(-theta * u) - 1)
    den <- (exp(-theta) - 1)^(d - 1)
    cdf <- -(1/theta) * log(1 + num / den)
  }
  return(cdf)
}

#Dependence modelling with Frank copula, assuming the same TPM for all strains
FrankJointTransitionMatrix_copula<- function(gamma, K, copulaParams){
  frankcopula_cached <- memoise::memoise(frank_cdf)

  S <- 2^K

  gamma[1, ] <- gamma[1, ][2:1]
  Gamma <- matrix(0, S, S)

  for(a in 0:(S-1)) {
    for(b in 0:(S-1)) {
      Indices <- integer(0)
      IndicesComplement <- integer(0)
      prob <- numeric(K)

      for(k in 1:K){
        from_k <- (a %/% 2^(k-1)) %% 2
        to_k<- (b %/% 2^(k-1)) %% 2

        if(from_k == 1) Indices <- c(Indices,k)
        else IndicesComplement <- c(IndicesComplement,k)

        prob[k] <- gamma[from_k+1, to_k+1]
      }

      subsets <- sets::set_power(sets::as.set(IndicesComplement))
      subsets <- lapply(subsets, function(g) unlist(as.vector(g)))

      total <- 0
      for(Tset in subsets){
        sign <- (-1)^length(Tset)
        idx <- c(Indices, Tset)
        u <- rep(1, K)
        if(length(idx)>0) u[idx] <- prob[idx]
        total <- total + sign * frank_cdf(u, copulaParams)
      }
      Gamma[a+1, b+1] <- total
    }
  }
  Gamma <- Gamma / rowSums(Gamma)
  Gamma
}

#Dependence modelling with Frank copula, assuming different TPM for all strains
FrankJointTransitionMatrix_copula_per_strain<- function(gammalist, K, copulaParams){
  cop<- copula::frankCopula(param = copulaParams, dim = K, use.indepC = "TRUE")

  S <- 2^K
  Gamma <- matrix(0, S, S)

  for(a in 0:(S-1)) {
    for(b in 0:(S-1)) {
      Indices <- integer(0)
      IndicesComplement <- integer(0)
      prob <- numeric(K)

      for(k in 1:K){
        from_k <- (a %/% 2^(k-1)) %% 2
        to_k<- (b %/% 2^(k-1)) %% 2

        if(from_k == 1) Indices <- c(Indices,k)
        else IndicesComplement <- c(IndicesComplement,k)
        gamma_k<- gammalist[[k]]
        gamma_k[1, ] <- gamma_k[1, ][2:1]
        prob[k] <- gamma_k[from_k+1, to_k+1]
      }

      subsets <- sets::set_power(sets::as.set(IndicesComplement))
      subsets <- lapply(subsets, function(g) unlist(as.vector(g)))

      total <- 0
      for(Tset in subsets){
        sign <- (-1)^length(Tset)
        idx <- c(Indices, Tset)
        u <- rep(1, K)
        if(length(idx)>0) u[idx] <- prob[idx]
        total <- total + sign * copula::pCopula(u, cop)
      }
      Gamma[a+1, b+1] <- total
    }
  }
  Gamma <- Gamma / rowSums(Gamma)
  Gamma
}

Multipurpose_JointTransitionMatrix2<- function(gammas, K, copParams, Modeltype){
  if(Modeltype==1){
    JointTPM<- JointTransitionMatrix_arma_cpp2(G(gammas[1], gammas[2]), K)
  }else if(Modeltype==2){
    Glist<- BuildGamma_list_cpp(gammas)
    JointTPM<- JointTransitionMatrix_per_strain_cpp2(Glist, K)
  }else if(Modeltype == 3){
    JointTPM<- JointTransitionMatrix_copula_cpp(G(gammas[1], gammas[2]), K, copParams)
  }else if(Modeltype == 4){
    Glist<- BuildGamma_list_cpp(gammas)
    JointTPM<- JointTransitionMatrix_copula_perstrain_cpp(Glist, K, copParams)
  }else if(Modeltype==5){
    JointTPM<- FrankJointTransitionMatrix_copula(G(gammas[1], gammas[2]), K, copParams)
  }else if(Modeltype==6){
    Glist<- BuildGamma_list_cpp(gammas)
    JointTPM<- FrankJointTransitionMatrix_copula_per_strain(Glist, K, copParams)
  }
  return(JointTPM)
}


#Dependence modelling with copula, assuming the same TPM for all strains
ParallelJointTransitionMatrix_copula<- function(gamma, K, copulaParams){
  S<- 2^K
  corrMat<- diag(K)
  uppertriang <- copulaParams
  gdata::upperTriangle(corrMat, byrow=TRUE) <- uppertriang
  gdata::lowerTriangle(corrMat, byrow=FALSE) <- uppertriang

  gamma[1, ]<- gamma[1, ][2:1]

  ncores <- parallel::detectCores() - 1
  cl <- parallel::makeCluster(ncores)
  doParallel::registerDoParallel(cl)

  Gamma <- matrix(0, nrow = S, ncol = S)
  library(foreach)
  Gamma_rows <- foreach::foreach(a = 0:(S - 1), .combine = rbind,
                        .export = "gaussian_copula_cdf",
                        .packages = c("sets")) %dopar% {

                          row_vals <- numeric(S)
                          for (b in 0:(S - 1)) {
                            Indices <- c()
                            IndicesComplement <- c()
                            prob <- numeric(K)
                            for (k in 1:K) {
                              from_k <- (a %/% 2^(k - 1)) %% 2
                              to_k   <- (b %/% 2^(k - 1)) %% 2
                              if (from_k == 1) {
                                Indices <- c(Indices, k)
                              } else {
                                IndicesComplement <- c(IndicesComplement, k)
                              }
                              prob[k] <- gamma[from_k + 1, to_k + 1]
                            }
                            subsets <- sets::set_power(sets::as.set(IndicesComplement))
                            subsets <- lapply(subsets, function(g) unlist(as.vector(g)))

                            total <- 0
                            for (Tset in subsets) {
                              sign <- (-1)^length(Tset)
                              idx <- c(Indices, Tset)
                              u <- rep(1, K)
                              if (length(idx) > 0) u[idx] <- prob[idx]
                              total <- total + sign * gaussian_copula_cdf(u, corrMat)
                            }
                            row_vals[b + 1] <- total
                          }
                          row_vals
                        }
  Gamma <- Gamma_rows
  parallel::stopCluster(cl)
  Gamma <- Gamma / rowSums(Gamma)
  return(Gamma)
}

#Dependence modelling with copula, assuming the same TPM for all strains - Vine
GaussianJointTransitionMatrix_copula_vine <- function(gamma, K, copulaParams){
  S <- 2^K
  corrMat <- diag(K)
  uppertriang <- copulaParams
  gdata::upperTriangle(corrMat, byrow = TRUE) <- uppertriang
  gdata::lowerTriangle(corrMat, byrow = FALSE) <- uppertriang

  struct <- rvinecopulib::cvine_structure(K)
  pclist <- vector("list", K-1)
  for(t in 1:(K-1)) {
    n_edges <- K - t
    pclist[[t]] <- lapply(1:n_edges, function(i) {
      rho <- corrMat[t, t+i]
      rvinecopulib::bicop_dist('gaussian', parameters = rho)
    })
  }
  vcop <- rvinecopulib::vinecop_dist(pclist, struct)

  gamma[1,] <- gamma[1,][2:1]

  Gamma <- matrix(0, nrow = S, ncol = S)
  for(a in 0:(S - 1)) {
    for(b in 0:(S - 1)){

      Indices <- c()
      IndicesComplement <- c()
      prob <- numeric(K)
      for(k in 1:K){
        from_k <- (a %/% 2^(k - 1)) %% 2
        to_k   <- (b %/% 2^(k - 1)) %% 2
        if(from_k == 1) {
          Indices <- c(Indices, k)
        } else {
          IndicesComplement <- c(IndicesComplement, k)
        }
        prob[k] <- gamma[from_k + 1, to_k + 1]
      }
      # inclusion-exclusion subsets
      subsets <- sets::set_power(sets::as.set(IndicesComplement))
      subsets <- lapply(subsets, function(g) unlist(as.vector(g)))

      total <- 0
      for(Tset in subsets) {
        sign <- (-1)^length(Tset)
        idx  <- c(Indices, Tset)
        u <- rep(1, K)
        if(length(idx) > 0) u[idx] <- prob[idx]
        total <- total + sign * rvinecopulib::pvinecop(u, vcop, n_mc = 10^4, cores = 7)
      }
      Gamma[a + 1, b + 1] <- total
    }
  }
  Gamma <- Gamma / rowSums(Gamma)
  return(Gamma)
}

#Dependence modelling with copula, assuming the same TPM for all strains - Vine
FrankJointTransitionMatrix_copula_vine <- function(gamma, K, copulaParams){
  S <- 2^K
  corrMat <- diag(K)
  uppertriang <- copulaParams
  gdata::upperTriangle(corrMat, byrow = TRUE) <- uppertriang
  gdata::lowerTriangle(corrMat, byrow = FALSE) <- uppertriang

  struct <- rvinecopulib::cvine_structure(K)
  pclist <- vector("list", K-1)
  for(t in 1:(K-1)) {
    n_edges <- K - t
    pclist[[t]] <- lapply(1:n_edges, function(i){
      rho <- corrMat[t, t+i]
      rvinecopulib::bicop_dist('frank', parameters = rho)
    })
  }
  vcop <- rvinecopulib::vinecop_dist(pclist, struct)

  gamma[1,] <- gamma[1,][2:1]

  Gamma <- matrix(0, nrow = S, ncol = S)
  for(a in 0:(S - 1)) {
    for(b in 0:(S - 1)){

      Indices <- c()
      IndicesComplement <- c()
      prob <- numeric(K)
      for(k in 1:K){
        from_k <- (a %/% 2^(k - 1)) %% 2
        to_k   <- (b %/% 2^(k - 1)) %% 2
        if(from_k == 1) {
          Indices <- c(Indices, k)
        } else {
          IndicesComplement <- c(IndicesComplement, k)
        }
        prob[k] <- gamma[from_k + 1, to_k + 1]
      }
      # inclusion-exclusion subsets
      subsets <- sets::set_power(sets::as.set(IndicesComplement))
      subsets <- lapply(subsets, function(g) unlist(as.vector(g)))

      total <- 0
      for(Tset in subsets) {
        sign <- (-1)^length(Tset)
        idx  <- c(Indices, Tset)
        u <- rep(1, K)
        if(length(idx) > 0) u[idx] <- prob[idx]
        total <- total + sign * rvinecopulib::pvinecop(u, vcop, n_mc = 10^4, cores = 7)
      }
      Gamma[a + 1, b + 1] <- total
    }
  }
  Gamma <- Gamma / rowSums(Gamma)
  return(Gamma)
}

#Dependence modelling with copula, assuming each strain has its unique TPM
JointTransitionMatrix_copula_per_strain<- function(gamma_list, copulaParam){
  K<- length(gamma_list)
  S<- 2^K
#  cop<- copula::normalCopula(param = copulaParam, dim = K)
  cop<- copula::plackettCopula(param = copulaParam)

  Gamma<- matrix(0, nrow = S, ncol = S)
  for(a in 0:(S - 1)) {
    for(b in 0:(S - 1)){
      Indices<- c()
      IndicesComplement<- c()
      prob<- numeric(K)
      for(k in 1:K){
        from_k<- (a %/% 2^(k - 1)) %% 2
        to_k<- (b %/% 2^(k - 1)) %% 2
        if(from_k == 1){
          Indices<- c(Indices, k)
        }else{
          IndicesComplement<- c(IndicesComplement, k)
        }
        gamma_k<- gamma_list[[k]]
        gamma_k[1, ]<- gamma_k[1, ][2:1]
        prob[k]<- gamma_k[from_k + 1, to_k + 1]
      }
      subsets<- sets::set_power(sets::as.set(IndicesComplement))
      subsets<- lapply(subsets, function(g) unlist(as.vector(g)))
      total <- 0
      for(Tset in subsets){
        sign <- (-1)^length(Tset)
        idx <- c(Indices, Tset)
        u <- rep(1, K)
        if(length(idx) > 0) u[idx] <- prob[idx]
        total <- total + sign * copula::pCopula(u, cop)
      }
      Gamma[a + 1, b + 1]<- total
    }
  }
  Gamma <- Gamma / rowSums(Gamma)
  return(Gamma)
}


#Build list of per-strain transition matrices
BuildGamma_list<- function(Gs){
  n_mat<-length(Gs)/2
  Gmats<- list()
  index<- 0
  for(i in 1:n_mat){
    Gmats[[i]] = G(Gs[i+index], Gs[i+index+1])
    index<- index + 1
  }
  return(Gmats)
}

encodeBits<- function(K){
S <- 2^K
bits<- matrix(0, nrow = S, ncol = K)
for (i in 0:(S - 1)) {
  for (k in 1:K) {
    bits[i + 1, k] <- (i %/% 2^(k - 1)) %% 2
  }
}
bits <- bits[, K:1]
return(bits)
}

logspace_vecmatmult<- function(log_v, log_M) {
  K <- length(log_v)
  log_result <- numeric(K)

  for (j in 1:K) {
    log_terms <- log_v + log_M[, j]
    log_result[j] <- matrixStats::logSumExp(c(log_terms))
  }
  return(log_result)
}

simulateMarkovChain<- function(nstep, TPM){
  state.space<- 0:(nrow(TPM)-1)
  initial_state <- 0
  states <- numeric(nstep)
  states[1] <- initial_state
  for(i in 2:nstep) {
    current_state <- states[i - 1]
    next_state <- sample(state.space, size = 1, prob = TPM[current_state + 1, ])
    states[i] <- next_state
  }
  return(states)
}

sim.RW2mean0<- function(time, sd=0.0001, init.r1 = 0, init.r2 = 0){
  r <- numeric(time)

  r[1] <- init.r1
  r[2] <- init.r2

  for(t in 3:time){
    epsilon <- rnorm(1, mean = 0, sd = sd)
    r[t] <- 2*(r[t - 1]) - r[t - 2] + epsilon
  }
  r <- r - mean(r)
  return(r)
}

Multstrain.simulate<- function(Model, time, nstrain=2, adj.matrix, Modeltype=1, copulaParam=c(-0.8,0.8,-0.6),
                               e_it=matrix(c(rep(c(rpois(time, 500000), rpois(time, 1000000)), 4), rpois(time, 500000)),
                                           byrow = T, ncol = time),
                               B = runif(nstrain), T.prob = matrix(c(0.9, 0.1, 0.2, 0.8), nrow = 2, byrow = T),
                               r = sim.RW2mean0(time, sd=0.009), s = DetectOutbreaks:::sim.Seasonals2(Amplitude = 1.4),
                               u = DetectOutbreaks:::sim.Spatials(adj.matrix)){
  ndept<- nrow(adj.matrix)
  y_itk<- array(NA, dim=c(ndept, time, nstrain))
  EpidemicIndicator<- matrix(NA, ndept, time)
  Jointstates<- 2^nstrain
  Bits<- encodeBits(K=nstrain)
  aVec<- numeric(nstrain)
  for(k in 1:nstrain){
    a_k<- runif(1, min = -14, max = -12)
    aVec[k]<- a_k
  }
  #Due to lazyloading
  r<- r;  s<- s;  u<- u;  e_it<- e_it;  B<- B


  if(Modeltype == 1){
    copulaParam<- 0
    T.prob<- c(T.prob[1,2],T.prob[2,1])
    JointTPM<- Multipurpose_JointTransitionMatrix2(T.prob, nstrain, copulaParam, Modeltype)
  }else if(Modeltype == 2){
    copulaParam<- 0
    T.prob<- runif(2*nstrain, min = 0.1, max = 0.2)
    JointTPM<- Multipurpose_JointTransitionMatrix2(T.prob, nstrain, copulaParam, Modeltype)
  }else if(Modeltype == 3){
    T.prob<- c(T.prob[1,2],T.prob[2,1])
    JointTPM<- Multipurpose_JointTransitionMatrix2(T.prob, nstrain, copulaParam, Modeltype)
  }else if(Modeltype == 4){
    T.prob<- runif(2*nstrain, min = 0.1, max = 0.2)
    JointTPM<- Multipurpose_JointTransitionMatrix2(T.prob, nstrain, copulaParam, Modeltype)
  }else if(Modeltype == 5){
    T.prob<- c(T.prob[1,2],T.prob[2,1])
    JointTPM<- Multipurpose_JointTransitionMatrix2(T.prob, nstrain, copulaParam, Modeltype)
    #JointTPM<- ParallelJointTransitionMatrix_copula(G(0.1,0.2), K=5, c(0.8,-0.85,0.9,-0.8,-0.86,0.87,-0.85,-0.8,0.8,-0.87))
    JointTPM<- ifelse(JointTPM<=0,1e-6,JointTPM)
    JointTPM<- ifelse(JointTPM>=1,1-1e-6,JointTPM)
  }else if(Modeltype == 6){
    T.prob<- runif(2*nstrain, min = 0.1, max = 0.2)
    JointTPM<- Multipurpose_JointTransitionMatrix2(T.prob, nstrain, copulaParam, Modeltype)
    JointTPM<- ifelse(JointTPM<=0,1e-6,JointTPM)
    JointTPM<- ifelse(JointTPM>=1,1-1e-6,JointTPM)
  }else if(Modeltype == 7){
    JointTPM<- matrix(NA, nrow = Jointstates, ncol = Jointstates)
    T.prob<- 0
    JointTPM<- gtools::rdirichlet(Jointstates, sample(2:7, size = Jointstates, replace = T))
  }

  if(Model == 0){
    for(i in 1:ndept){
      for(t in 1:time){
        m<- (t - 1) %% 12 + 1
        for(k in 1:nstrain){
          lograte <- aVec[k] + r[t] + s[m] + u[i]
          y_itk[i, t, k]<- rpois(1, lambda = e_it[i, t] * exp(lograte))
        }
      }
    }
    return(list("y" =y_itk, "e_it"=e_it, "r"=r, "s"=s, "u"=u, "states"=EpidemicIndicator, "a_k"=aVec))
  }
  else{
    for(i in 1:ndept){
      for(k in 1:nstrain){
        lograte<- aVec[k] + r[1] + s[1] + u[i]
        y_itk[i, 1, k]<- rpois(1, lambda = e_it[i, 1] * exp(lograte))
      }
      EpidemicIndicator[i, ]<- simulateMarkovChain(nstep=time, JointTPM)
    }

    for(t in 2:time){
      m<- (t - 1) %% 12 + 1
      for(i in 1:ndept){
        for(k in 1:nstrain){
          newB<- rep(0, nstrain)
          newB[k]<- B[k]
          lograte<- aVec[k] + r[t] + s[m] + u[i] + (newB %*% Bits[EpidemicIndicator[i, t]+1, ])
          y_itk[i, t, k]<- rpois(1, lambda = e_it[i, t] * exp(lograte))
        }
      }
    }
    return(list("y" =y_itk, "e_it"=e_it, "r"=r, "s"=s, "u"=u, "states"=EpidemicIndicator, "B"=B, "a_k"=aVec, "T.prob"=T.prob, "JointTPM"=JointTPM, "copulaParam"=copulaParam))
  }
}


perstrainOutbreaks <- function(EpidemicIndicator, nstrain) {
  Bits <- encodeBits(K = nstrain)
  ndept <- nrow(EpidemicIndicator)
  time <- ncol(EpidemicIndicator)

  Outbreaks <- array(NA, dim = c(ndept, time, nstrain))

  for (i in 1:ndept) {
    for (t in 1:time) {
      s <- EpidemicIndicator[i, t]   #joint state index
      Outbreaks[i, t, ] <- Bits[s + 1, ]  #per-strain indicators
    }
  }
  return(Outbreaks)
}


stationarydist <- function(Gamma) {
  mT <- t(Gamma)
  eig <- eigen(mT)
  E_values <- eig$values
  E_vectors <- eig$vectors

  # Find the eigenvalue closest to 1 and get its index
  distances <- abs(Re(E_values) - 1)
  index <- which.min(distances)

  stationary_distribution <- Re(E_vectors[, index])
  stationary_distribution <- stationary_distribution / sum(stationary_distribution)
  return(stationary_distribution)
}

multstrain.forwardfilter<- function(y, e_it, nstrain, r, s, u, JointTPM, B, Bits, a_k){
  ndept<- length(u)
  time <- length(r)
  nstate<- 2^nstrain
  AllForwardprobs<- vector("list", ndept)
  loginit.density<- log(stationarydist(JointTPM))
  for(i in 1:ndept){
    forwardprobs<- matrix(NA, nrow = time, ncol = nstate)
    logprodEmission<- rep(0, nstate)
  for(n in 1:nstate){
    for(k in 1:nstrain){
      logprodEmission[n]<- logprodEmission[n] + dpois(y[i, 1, k], lambda = e_it[i, 1] * exp(a_k[k] + r[1] + s[1] + u[i] + B%*%Bits[n, ]), log = TRUE)
    }
  }
  forwardprobs[1, ]<- logprodEmission + loginit.density
  for(t in 2:time){
    month_index<- (t-1) %% 12 + 1
    logprodEmission<- rep(0, nstate)
    for(n in 1:nstate){
      for(k in 1:nstrain){
        newB<- rep(0, nstrain)
        newB[k]<- B[k]
        logprodEmission[n]<- logprodEmission[n] + dpois(y[i, t, k], lambda = e_it[i, t] * exp(a_k[k] + r[t] + s[month_index] + u[i] + newB%*%Bits[n, ]), log = TRUE)
      }
    }
    forwardprobs[t, ]<- logspace_vecmatmult(forwardprobs[t-1, ], log(JointTPM)) + logprodEmission
  }
  AllForwardprobs[[i]]<- forwardprobs
  }
  return(AllForwardprobs)
}

multstrain.backwardsweep<- function(y, e_it, nstrain, r, s, u, JointTPM, B, Bits, a_k){
  ndept<- length(u)
  time <- length(r)
  nstate<- 2^nstrain
  Allbackwardprob<- vector("list", ndept)

    for (i in 1:ndept) {
      backwardprobs<- matrix(NA, nrow = time, ncol = nstate)
      backwardprobs[time, ] <- rep(0, nstate)

      for (t in (time-1):1) {
          month_index<- t %% 12 + 1
          logprodEmission<- rep(0, nstate)
          for(n in 1:nstate){
            for(k in 1:nstrain){
              newB<- rep(0, nstrain)
              newB[k]<- B[k]
              logprodEmission[n]<- logprodEmission[n] + dpois(y[i, t+1, k], lambda = e_it[i, t+1] * exp(a_k[k] + r[t+1] + s[month_index] + u[i] + newB%*%Bits[n, ]), log = TRUE)
            }
          }
          backwardprobs[t, ]<- logspace_vecmatmult(logprodEmission + backwardprobs[t+1, ], log(t(JointTPM)))
      }
      Allbackwardprob[[i]] <- backwardprobs
    }
    return(Allbackwardprob)
}

multstrain.Decoding <- function(y, e_it, nstrain, r, s, u, Gamma, copParam, B, Bits, a_k, state, Modeltype) {
  ndept<- length(u)
  time <- length(r)
  nstate<- 2^nstrain
  if(Modeltype == 1){
    JointTPM<- Multipurpose_JointTransitionMatrix2(Gamma, nstrain, copParam, Modeltype)
  }else if(Modeltype == 2){
    JointTPM<- Multipurpose_JointTransitionMatrix2(Gamma, nstrain, copParam, Modeltype)
  }else if(Modeltype == 3){
    JointTPM<- Multipurpose_JointTransitionMatrix2(Gamma, nstrain, copParam, Modeltype)
    JointTPM<- ifelse(JointTPM<=0,1e-6,JointTPM)
    JointTPM<- ifelse(JointTPM>=1,1-1e-6,JointTPM)
  }else if(Modeltype == 4){
    JointTPM<- Multipurpose_JointTransitionMatrix2(Gamma, nstrain, copParam, Modeltype)
    JointTPM<- ifelse(JointTPM<=0,1e-6,JointTPM)
    JointTPM<- ifelse(JointTPM>=1,1-1e-6,JointTPM)
  }else if(Modeltype == 5){
    JointTPM<- Multipurpose_JointTransitionMatrix2(Gamma, nstrain, copParam, Modeltype)
    JointTPM<- ifelse(JointTPM<=0,1e-6,JointTPM)
    JointTPM<- ifelse(JointTPM>=1,1-1e-6,JointTPM)
  }else if(Modeltype == 6){
    JointTPM<- Multipurpose_JointTransitionMatrix2(Gamma, nstrain, copParam, Modeltype)
    JointTPM<- ifelse(JointTPM<=0,1e-6,JointTPM)
    JointTPM<- ifelse(JointTPM>=1,1-1e-6,JointTPM)
  }else if(Modeltype == 7){
    JointTPM<- matrix(Gamma, nrow = nstate, byrow = TRUE)
  }
    Allforwardprobs<- multstrain.forwardfilter(y, e_it, nstrain, r, s, u, JointTPM, B, Bits, a_k)
    Allbackwardprobs<- multstrain.backwardsweep(y, e_it, nstrain, r, s, u, JointTPM, B, Bits, a_k)
    Res<- matrix(NA, ndept, time)
    for(i in 1:ndept){
      for(j in 1:time){
        post_log <- Allforwardprobs[[i]][j, ] + Allbackwardprobs[[i]][j, ]
        denom <- matrixStats::logSumExp(post_log)
        posterior <- exp(post_log - denom)
        Res[i, j] <- posterior[state]
      }
    }
    return(Res)
}

R_dpois<- function(y, lambda){
  res <- (y * log(lambda)) -lambda -log(factorial(y))
  return(res)
}

Posteriormultstrain.Decoding<- function(y, e_it, inf.object, Modeltype, thinningL=1000, burn.in=1000){
  ndept<- dim(y)[1]
  time<- dim(y)[2]
  nstrain<- dim(y)[3]
  nstate<- 2^nstrain
  Bits<- encodeBits(K=nstrain)

  if(!is.data.frame(inf.object)){
    fullG12.draws<- stack(as.data.frame(inf.object$draws(variables = "G12")[,1,]))[,1]
    fullG21.draws<- stack(as.data.frame(inf.object$draws(variables = "G21")[,1,]))[,1]
    fullr.draws<- as.data.frame(inf.object$draws(variables = "r")[,1,])
    fulls.draws<- as.data.frame(inf.object$draws(variables = "s")[,1,])
    fullu.draws<- as.data.frame(inf.object$draws(variables = "uconstrained")[,1,])
    fullB.draws<- as.data.frame(inf.object$draws(variables = "B")[,1,])
    fulla_k.draws<- as.data.frame(inf.object$draws(variables = "a_k")[,1,])
  }else{
    fullG.draws <- inf.object[-(1:burn.in), startsWith(colnames(inf.object), "G")]
    fullr.draws<- inf.object[-(1:burn.in), startsWith(colnames(inf.object), "r")]
    fulls.draws<- inf.object[-(1:burn.in), startsWith(colnames(inf.object), "s")]
    fullu.draws<- inf.object[-(1:burn.in), startsWith(colnames(inf.object), "u")]
    fullB.draws<- inf.object[-(1:burn.in), startsWith(colnames(inf.object), "B")]
    fulla_k.draws<- inf.object[-(1:burn.in), startsWith(colnames(inf.object), "a")]
    if(Modeltype %in% c(3,4,5,6)) fullcop.draws<- inf.object[-(1:burn.in), startsWith(colnames(inf.object), "c")]
  }

  thinning<- numeric(floor(nrow(fullr.draws)/thinningL))
  thinning[1]<- thinningL
  for(i in 2:length(thinning)){
    thinning[i]<- thinning[i-1] + thinningL
  }

  G.draws<- fullG.draws[thinning, ]
  r.draws<- fullr.draws[thinning, ]
  s.draws<- fulls.draws[thinning, ]
  u.draws<- fullu.draws[thinning, ]
  B.draws<- fullB.draws[thinning, ]
  a_k.draws<- fulla_k.draws[thinning, ]
  if(Modeltype %in% c(3,4)){
    cop.draws<- fullcop.draws[thinning, ]
  }else if(Modeltype %in% c(5,6)){
    cop.draws<- fullcop.draws[thinning]
  }

  decodedOutbreakMatrix<- list()
  for(n in 1:nstate){
    sum_Xit<- matrix(0, nrow = ndept, ncol = time)

    for(index in 1:length(thinning)){
      Gs<- as.numeric(G.draws[index,])
      r<- as.numeric(r.draws[index,])
      s<- as.numeric(s.draws[index,])
      u<- as.numeric(u.draws[index,])
      B<- as.numeric(B.draws[index,])
      a_k<- as.numeric(a_k.draws[index,])
      if(Modeltype %in% c(3,4)){
        cop<- as.numeric(cop.draws[index, ])
      }else if(Modeltype %in% c(5,6)){
        cop<- cop.draws[index]
      }else{
        cop<- 0
      }
      Ex_Xit <- multstrain.Decoding(y = y, e_it = e_it, nstrain=nstrain, r = r, s = s, u = u, Gamma = Gs, copParam=cop, Bits = Bits, B = B, a_k = a_k, state = n, Modeltype = Modeltype)
      sum_Xit<- sum_Xit + Ex_Xit
    }
    decodedOutbreakMatrix[[n]]<- sum_Xit/length(thinning)
  }

  perStrainProbs <- array(0, dim = c(ndept, time, nstrain))
  for(n in 1:nstate){
    for(k in 1:nstrain){
      perStrainProbs[,,k] <- perStrainProbs[,,k] + Bits[n, k] * decodedOutbreakMatrix[[n]]
    }
  }
  return(perStrainProbs)
}

perstrainPosteriorOutbreaks <- function(decodedOutbreakMatrix, nstrain){
  time<- ncol(decodedOutbreakMatrix[[1]])
  ndept<- nrow(decodedOutbreakMatrix[[1]])
  perStrainProbs <- array(0, dim = c(ndept, time, nstrain))
  nstate<- 2^nstrain
  Bits<- encodeBits(K=nstrain)
  for(n in 1:nstate){
    for(k in 1:nstrain){
      perStrainProbs[,,k] <- perStrainProbs[,,k] + Bits[n, k] * decodedOutbreakMatrix[[n]]
    }
  }
  return(perStrainProbs)
}


multstrainLoglikelihood<- function(y, e_it, nstrain, r, s, u, Gamma, B, Bits, a_k, Model, independentChains=FALSE){
  ndept<- nrow(y[,,1])
  time <- ncol(y[,,1])
  if(Model == 0){
    loglike<- 0
    for(i in 1:ndept){
      for(t in 1:time){
        month_index<- (t-1) %% 12 + 1
        for(k in 1:nstrain){
          loglike<- loglike + dpois(y[i, t, k], lambda = e_it[i, t] * exp(a_k[k] + r[t] + s[month_index] + u[i]), log=T)
        }
      }
    }
    return(loglike)
  }else{
    if(independentChains){
      nstate<- ncol(Gamma)
      AllForwardprobs<- numeric(ndept)
      loginit.density<- log(stationarydist(Gamma))
      for(i in 1:ndept){
        forwardprobs<- matrix(NA, nrow = time, ncol = nstate)
        logprodEmission<- rep(0, nstate)
        for(n in 1:nstate){
          for(k in 1:nstrain){
            logprodEmission[n]<- logprodEmission[n] + dpois(y[i, 1, k], lambda = e_it[i, 1] * exp(a_k[k] + r[1] + s[1] + u[i] + B[k] * (n-1)), log=T)
          }
        }
        forwardprobs[1, ]<- logprodEmission + loginit.density
        for(t in 2:time){
          month_index<- (t-1) %% 12 + 1
          logprodEmission<- rep(0, nstate)
          for(n in 1:nstate){
            for(k in 1:nstrain){
              logprodEmission[n]<- logprodEmission[n] + dpois(y[i, t, k], lambda = e_it[i, t] * exp(a_k[k] + r[t] + s[month_index] + u[i] + B[k] * (n-1)), log=T)
            }
          }
          forwardprobs[t, ]<- logspace_vecmatmult(forwardprobs[t-1, ], log(Gamma)) + logprodEmission
        }
        AllForwardprobs[i]<- logSumExp_cpp(forwardprobs[time, ])
        #print(AllForwardprobs[[i]])
      }
      return(sum(AllForwardprobs))
  }else{
  nstate<- 2^nstrain
  AllForwardprobs<- numeric(ndept)
  JointTPM<- JointTransitionMatrix(gamma = Gamma, K = nstrain)
  loginit.density<- log(stationarydist(JointTPM))
  for(i in 1:ndept){
    forwardprobs<- matrix(NA, nrow = time, ncol = nstate)
    logprodEmission<- rep(0, nstate)
    for(n in 1:nstate){
      for(k in 1:nstrain){
        newB<- rep(0, nstrain)
        newB[k]<- B[k]
        logprodEmission[n]<- logprodEmission[n] + dpois(y[i, 1, k], lambda = e_it[i, 1] * exp(a_k[k] + r[1] + s[1] + u[i] + as.numeric(newB%*%Bits[n, ])), log=T)
      }
    }
    forwardprobs[1, ]<- logprodEmission + loginit.density
    for(t in 2:time){
      month_index<- (t-1) %% 12 + 1
      logprodEmission<- rep(0, nstate)
      for(n in 1:nstate){
        for(k in 1:nstrain){
          newB<- rep(0, nstrain)
          newB[k]<- B[k]
          logprodEmission[n]<- logprodEmission[n] + dpois(y[i, t, k], lambda = e_it[i, t] * exp(a_k[k] + r[t] + s[month_index] + u[i] + as.numeric(newB%*%Bits[n, ])), log=T)
        }
      }
      forwardprobs[t, ]<- logspace_vecmatmult(forwardprobs[t-1, ], log(JointTPM)) + logprodEmission
    }
    AllForwardprobs[i]<- logSumExp_cpp(forwardprobs[time, ])
    #print(AllForwardprobs[[i]])
  }
  return(sum(AllForwardprobs))
   }
  }
}

mcmc.plot2<- function(inf.object){

  if(is.data.frame(inf.object)){

    # par(mfrow=c(3, 3))
    #  for (i in 1:ncol(inf.object)) {
    #    hist(inf.object[-(1:2000), i], main = colnames(inf.object)[i], xlab ="", col = "white", border = "black")
    # }

    par(mfrow=c(3, 3))
    for (i in 1:ncol(inf.object)) {
      plot(inf.object[, i], type = "l", main = colnames(inf.object)[i], xlab ="MCMC iterations", ylab = "", col = "purple")
      grid()
    }
  }else{
    cond<- checkB(inf.object)
    if(cond){
      inf.object<- inf.object$draws(variables = c("G12", "G21", "kappa_r", "kappa_s", "kappa_u", "r", "s", "uconstrained", "B", "a_k", "state1_stationary_dist", "stationaryDistribution", "log_lik"))
    }else{
      inf.object<- inf.object$draws(variables = c("G12", "G21", "kappa_r", "kappa_s", "kappa_u", "r", "s", "uconstrained", "a_k", "state1_stationary_dist", "log_lik"))
    }
    inf.object<- inf.object[,1,]
    inf.object<- as.data.frame(inf.object)
    colnames(inf.object) <- gsub("^1.", "", colnames(inf.object))

    par(mfrow=c(3, 3))
    for (i in 1:ncol(inf.object)) {
      hist(inf.object[, i], main = colnames(inf.object)[i], xlab ="", col = "white", border = "black")
    }

    par(mfrow=c(3, 3))
    for (i in 1:ncol(inf.object)) {
      plot(inf.object[, i], type = "l", main = colnames(inf.object)[i], xlab ="HMC iterations", ylab = "", col = "red")
      grid()
    }
  }
}

scalingmultstrainLoglikelihood<- function(y, e_it, nstrain, r, s, u, Gamma, B, Bits, a_k, Model, independentChains=FALSE){
  ndept<- nrow(y[,,1])
  time <- ncol(y[,,1])
  if(Model == 0){
    loglike<- 0
    for(i in 1:ndept){
      for(t in 1:time){
        month_index<- (t-1) %% 12 + 1
        for(k in 1:nstrain){
          loglike<- loglike + dpois(y[i, t, k], lambda = e_it[i, t] * exp(a_k[k] + r[t] + s[month_index] + u[i]), log=T)
        }
      }
    }
    return(loglike)
  }else{
    if(independentChains){
      nstate<- ncol(Gamma)
      AllForwardprobs<- numeric(ndept)
      loginit.density<- log(stationarydist(Gamma))
      for(i in 1:ndept){
        forwardprobs<- matrix(NA, nrow = time, ncol = nstate)
        logprodEmission<- rep(0, nstate)
        for(n in 1:nstate){
          for(k in 1:nstrain){
            logprodEmission[n]<- logprodEmission[n] + dpois(y[i, 1, k], lambda = e_it[i, 1] * exp(a_k[k] + r[1] + s[1] + u[i] + B[k] * (n-1)),log=T)
          }
        }
        forwardprobs[1, ]<- logprodEmission + loginit.density
        for(t in 2:time){
          month_index<- (t-1) %% 12 + 1
          logprodEmission<- rep(0, nstate)
          for(n in 1:nstate){
            for(k in 1:nstrain){
              logprodEmission[n]<- logprodEmission[n] + dpois(y[i, t, k], lambda = e_it[i, t] * exp(a_k[k] + r[t] + s[month_index] + u[i] + B[k] * (n-1)), log = T)
            }
          }
          forwardprobs[t, ]<- logspace_vecmatmult(forwardprobs[t-1, ], log(Gamma)) + logprodEmission
        }
        AllForwardprobs[i]<- logSumExp_cpp(forwardprobs[time, ])
        #print(AllForwardprobs[[i]])
      }
      return(sum(AllForwardprobs))
    }else{
      nstate<- 2^nstrain
      AllForwardprobs<- numeric(ndept)
      JointTPM<- JointTransitionMatrix(gamma = Gamma, K = nstrain)
      init.density<- stationarydist(JointTPM)
      for(i in 1:ndept){
        forwardprobs<- matrix(NA, nrow = time, ncol = nstate)
        prodEmission<- matrix(rep(1, time*nstate), nrow = time, ncol = nstate)
        init.density<- init.density/sum(init.density)
        for(n in 1:nstate){
          for(k in 1:nstrain){
            newB<- rep(0, nstrain)
            newB[k]<- B[k]
            prodEmission[1, n]<- prodEmission[1, n] * dpois(y[i, 1, k], lambda = e_it[i, 1] * exp(a_k[k] + r[1] + s[1] + u[i] + as.numeric(newB%*%Bits[n, ])))
          }
        }
        prodEmission[1, ]<- init.density * prodEmission[1, ]
        logsumforwardprobs<- log(sum(prodEmission[1, ]))
        prodEmission[1, ]<- prodEmission[1, ]/sum(prodEmission[1, ])
        forwardprobs[1, ]<- prodEmission[1, ]
        for(t in 2:time){
          month_index<- (t-1) %% 12 + 1
          for(n in 1:nstate){
            for(k in 1:nstrain){
              newB<- rep(0, nstrain)
              newB[k]<- B[k]
              prodEmission[t, n]<- prodEmission[t, n] * dpois(y[i, t, k], lambda = e_it[i, t] * exp(a_k[k] + r[t] + s[month_index] + u[i] + as.numeric(newB%*%Bits[n, ])))
            }
          }
          forwardprobs[t, ]<- (forwardprobs[t-1, ] %*% (JointTPM) * prodEmission[t, ])/sum(forwardprobs[t-1, ] %*% (JointTPM) * prodEmission[t, ])
          logsumforwardprobs<- logsumforwardprobs + log(sum(forwardprobs[t-1, ] %*% (JointTPM) * prodEmission[t, ]))
        }
        AllForwardprobs[i]<- logsumforwardprobs
        #print(AllForwardprobs[[i]])
      }
      return(sum(AllForwardprobs))
    }
  }
}

newscalingmultstrainLoglikelihood<- function(y, e_it, nstrain, r, s, u, Gamma, B, Bits, a_k, Model, independentChains=FALSE){
  ndept<- nrow(y[,,1])
  time <- ncol(y[,,1])
  if(Model == 0){
    loglike<- 0
    for(i in 1:ndept){
      for(t in 1:time){
        month_index<- (t-1) %% 12 + 1
        for(k in 1:nstrain){
          loglike<- loglike + dpois(y[i, t, k], lambda = e_it[i, t] * exp(a_k[k] + r[t] + s[month_index] + u[i]), log=T)
        }
      }
    }
    return(loglike)
  }else{
      nstate<- 2^nstrain
      loglikelihood<- numeric(ndept)
      gradR<- matrix(NA, nrow = ndept, ncol = time)
      grad_r<- matrix(NA, nrow=ndept, ncol=time)
      grad_s<- matrix(NA, nrow=ndept, ncol=12)
      grad_u<- matrix(NA, nrow=ndept, ncol=ndept)
      JointTPM<- JointTransitionMatrix(gamma = Gamma, K = nstrain)
      init.density<- stationarydist(JointTPM)

      for(i in 1:ndept){
        scalingFactor<- numeric(time)
        gradR_scalingFactor<- numeric(time)
        gradS_scalingFactor<- numeric(12)
        grad_forward<- array(rep(1, time*time*nstate), dim=c(time, time, nstate))
        forwardprobs<- matrix(NA, nrow = time, ncol = nstate)
        gradR_forward<- matrix(rep(1, time*nstate), nrow = time, ncol = nstate)
        gradS_forward<- matrix(NA, nrow = time, ncol = nstate)
        gradU_scalingFactor<- matrix(NA, nrow = time, ncol = nstate)
        prodEmission<- matrix(rep(1, time*nstate), nrow = time, ncol = nstate)
        init.density<- init.density
        for(n in 1:nstate){
          for(k in 1:nstrain){
            newB<- rep(0, nstrain)
            newB[k]<- B[k]
            Lambda<- exp(a_k[k] + r[1] + s[1] + u[i] + as.numeric(newB%*%Bits[n, ]))
            prodEmission[1, n]<- prodEmission[1, n] * dpois(y[i, 1, k], lambda = e_it[i, 1] * Lambda)
            gradR_forward[1, n]<- gradR_forward[1, n] * (e_it[i, 1]/factorial(y[i, 1, k])) * exp(e_it[i, 1] * Lambda) * Lambda * (1 - e_it[i, 1] * Lambda)
          }
        }
        for(gradt in 1:time){
          month_index<- (gradt-1) %% 12 + 1
          for(n in 1:nstate){
            for(k in 1:nstrain){
              newB<- rep(0, nstrain)
              newB[k]<- B[k]
              Lambda<- exp(a_k[k] + r[gradt] + s[month_index] + u[i] + as.numeric(newB%*%Bits[n, ]))
              grad_forward[1, gradt, n]<- grad_forward[1, gradt, n] * (e_it[i, gradt]/factorial(y[i, gradt, k])) * exp(e_it[i, gradt] * Lambda) * Lambda * (1 - e_it[i, gradt] * Lambda)
            }
          }
        }
        prodEmission[1, ]<- init.density * prodEmission[1, ]
        forwardprobs[1, ]<- prodEmission[1, ]
        scalingFactor[1]<- sum(forwardprobs[1,])
        grad_forward[1,1,]<- grad_forward[1,1,] * init.density
        for(t in 2:time){
          month_index<- (t-1) %% 12 + 1
          for(n in 1:nstate){
            for(k in 1:nstrain){
              newB<- rep(0, nstrain)
              newB[k]<- B[k]
              Lambda<- exp(a_k[k] + r[t] + s[month_index] + u[i] + as.numeric(newB%*%Bits[n, ]))
              prodEmission[t, n]<- prodEmission[t, n] * dpois(y[i, t, k], lambda = e_it[i, t] * Lambda)
            }
          }
          forwardprobs[t, ]<- ((forwardprobs[t-1, ] %*% JointTPM) * prodEmission[t, ])/sum(forwardprobs[t-1, ])
          scalingFactor[t]<- sum(forwardprobs[t, ])
          for(gradt in 1:time){
            month_index<- (gradt-1) %% 12 + 1
            derivRt<- rep(1, nstate)
            for(n in 1:nstate){
              for(k in 1:nstrain){
                newB<- rep(0, nstrain)
                newB[k]<- B[k]
                Lambda<- exp(a_k[k] + r[gradt] + s[month_index] + u[i] + as.numeric(newB%*%Bits[n, ]))
                derivRt[n]<- derivRt[n] * (e_it[i, t]/factorial(y[i, t, k])) * exp(e_it[i, t] * Lambda) * Lambda * (1 - e_it[i, t] * Lambda)
              }
            }
            grad_forward[t, gradt, ]<- (((grad_forward[t-1, gradt, ] %*% JointTPM) * prodEmission[t, ]) + (forwardprobs[t-1, ] %*% JointTPM) * derivRt)/sum(forwardprobs[t-1, ])
          }
        }
        loglikelihood[i]<- sum(log(scalingFactor))
        gradR[i, ]<- rowSums(grad_forward[time,,])/scalingFactor[time]
      }
      return(list("loglike"=sum(loglikelihood), "grad_r"= sum(gradR)))
  }
}

gradmultstrainLoglikelihood<- function(y, e_it, nstrain, r, s, u, Gamma, B, Bits, a_k, Model, Q_r, Q_s, Q_u){
  ndept <- nrow(e_it)
  time  <- ncol(e_it)
  nstrain<- dim(y)[3]
  if(Model==0){
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
    nstate <- 2^nstrain
    JointTPM <- JointTransitionMatrix(gamma = Gamma, K = nstrain)
    loginit.density <- log(stationarydist(JointTPM))
    AllForwardprobs <- numeric(ndept)

    grad_r <- numeric(time)
    grad_s <- numeric(12)
    grad_u <- numeric(ndept)

    for(i in 1:ndept){
      f <- matrix(-Inf, nrow=time, ncol=nstate)
      g_r <- array(0, dim=c(time, nstate, time))
      g_s <- array(0, dim=c(time, nstate, 12))
      g_u <- matrix(0, nrow=time, ncol=nstate)

      # initialization
      for(n in 1:nstate){
        lambda <- sapply(1:nstrain, function(k)
          e_it[i,1] * exp(a_k[k] + r[1] + s[1] + u[i] + as.numeric(B %*% Bits[n, ]))
        )
        logEm <- sum(dpois(y[i,1,], lambda=lambda, log=TRUE))
        f[1,n] <- loginit.density[n] + logEm
        emission_grad <- sum(y[i,1,] - lambda)
        g_r[1,n,1] <- emission_grad
        g_s[1,n,1] <- emission_grad
        g_u[1,n]   <- emission_grad
      }

      # recursion
      for(t in 2:time){
        month_index <- (t-1) %% 12 + 1
        for(n in 1:nstate){
          lambda <- sapply(1:nstrain, function(k)
            e_it[i,t] * exp(a_k[k] + r[t] + s[month_index] + u[i] + as.numeric(B %*% Bits[n, ]))
          )
          logEm <- sum(dpois(y[i,t,], lambda=lambda, log=TRUE))
          logtrans <- f[t-1,] + log(JointTPM[,n])
          f[t,n] <- logSumExp_cpp(logtrans) + logEm

          weights <- exp(logtrans - logSumExp_cpp(logtrans))

          # propagate derivatives
          for(tau in 1:(t-1)){
            g_r[t,n,tau] <- sum(weights * g_r[t-1,,tau])
          }
          for(m in 1:12){
            g_s[t,n,m] <- sum(weights * g_s[t-1,,m])
          }
          g_u[t,n] <- sum(weights * g_u[t-1,])

          # add emission gradient
          emission_grad <- sum(y[i,t,] - lambda)
          g_r[t,n,t] <- g_r[t,n,t] + emission_grad
          g_s[t,n,month_index] <- g_s[t,n,month_index] + emission_grad
          g_u[t,n] <- g_u[t,n] + emission_grad
        }
      }

      AllForwardprobs[i] <- logSumExp_cpp(f[time,])
      weightsT <- exp(f[time,] - AllForwardprobs[i])

      for(tau in 1:time){
        grad_r[tau] <- grad_r[tau] + sum(weightsT * g_r[time,,tau])
      }
      for(m in 1:12){
        grad_s[m] <- grad_s[m] + sum(weightsT * g_s[time,,m])
      }
      grad_u[i] <- grad_u[i] + sum(weightsT * g_u[time,])
    }
    grad_r<- grad_r - as.numeric(Q_r %*% r)
    grad_s<- grad_s - as.numeric(Q_s %*% s)
    grad_u<- grad_u - as.numeric(Q_u %*% u)
    loglike <- sum(AllForwardprobs)
    return(list(loglike=loglike, grad_r=grad_r, grad_s=grad_s, grad_u=grad_u))
  }
}

#Gradients, inverse-fisher, loglikelihood
FIRSTgradmultstrainLoglikelihood2<- function(y, e_it, nstrain, r, s, u, Gamma, B, Bits, a_k, Model, Q_r, Q_s, Q_u){
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
    grad_r <- numeric(time)
    grad_s <- numeric(12)
    grad_u <- numeric(ndept)

    Fisher_r <- matrix(0, nrow = time, ncol = time)
    Fisher_s <- matrix(0, nrow = 12,   ncol = 12)
    Fisher_u <- matrix(0, nrow = ndept, ncol = ndept)

    loglike_total <- 0

    JointTPM <- JointTransitionMatrix(gamma = Gamma, K = nstrain)
    logJointTPM <- log(JointTPM)

    for(i in 1:ndept){
      logEmissions <- matrix(NA, nrow = time, ncol = nstate)
      lambda_array  <- array(0, dim = c(time, nstate, nstrain))

      for(t in 1:time){
        month_index <- (t-1) %% 12 + 1
        for(n in 1:nstate){
          for(k in 1:nstrain){
            newB<- rep(0, nstrain)
            newB[k]<- B[k]
            lambda_array[t,n,k] <- e_it[i,t] * exp(a_k[k] + r[t] + s[month_index] + u[i] + as.numeric(newB %*% Bits[n, ]))
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

      loggamma <- matrix(NA, nrow = time, ncol = nstate)
      gamma    <- matrix(NA, nrow = time, ncol = nstate)
      for(t in 1:time){
        loggamma[t, ] <- logalpha[t, ] + logbeta[t, ] - loglik_i
        gamma[t, ] <- exp(loggamma[t, ])
        gamma[t, ] <- gamma[t, ] / sum(gamma[t, ])
      }

      E_lambda_tk <- matrix(0, nrow = time, ncol = nstrain)
      for(t in 1:time){
        for(k in 1:nstrain){
          E_lambda_tk[t,k] <- sum(gamma[t, ] * lambda_array[t, , k])
        }
      }

      for(t in 1:time){
        month_index <- (t-1) %% 12 + 1
        for(k in 1:nstrain){
          Elambda <- E_lambda_tk[t,k]
          grad_r[t] <- grad_r[t] + (y[i,t,k] - Elambda)
          grad_s[month_index] <- grad_s[month_index] + (y[i,t,k] - Elambda)
          grad_u[i] <- grad_u[i] + (y[i,t,k] - Elambda)
          Fisher_r[t,t] <- Fisher_r[t,t] + Elambda
          Fisher_s[month_index, month_index] <- Fisher_s[month_index, month_index] + Elambda
          Fisher_u[i,i] <- Fisher_u[i,i] + Elambda
        }
      }
    }
    grad_r<- grad_r - as.numeric(Q_r %*% r)
    grad_s<- grad_s - as.numeric(Q_s %*% s)
    grad_u<- grad_u - as.numeric(Q_u %*% u)
    cov_r<- solve(Fisher_r + Q_r + diag(1e-8, time))
    cov_s<- solve(Fisher_s + Q_s)

    return(list(loglike = loglike_total, grad_r = grad_r, grad_s = grad_s, grad_u = grad_u, cov_r = cov_r,cov_s = cov_s))
  }
}
