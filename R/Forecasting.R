forecast<- function(Model, inf.object, forecastlength, y, popdata, adj.matrix, thinningL=10, burn.in=2000){
  ndept<- nrow(popdata)
  time<- ncol(popdata)
  OutP<- DetectOutbreaks::OutbreakProbability(y=y,e_it=popdata,inf.object = inf.object,adjmat = adj.matrix,Model = Model)
  OutP<- OutP[,time]

  if(!is.data.frame(inf.object)){
    fullG12.draws<- stack(as.data.frame(inf.object$draws(variables = "G12")[,1,]))[,1]
    fullG21.draws<- stack(as.data.frame(inf.object$draws(variables = "G21")[,1,]))[,1]
    fullkappaR.draws<- stack(as.data.frame(inf.object$draws(variables = "kappa_r")[,1,]))[,1]
    fullr.draws<- as.data.frame(inf.object$draws(variables = "r")[,1,])
    fulls.draws<- as.data.frame(inf.object$draws(variables = "s")[,1,])
    fullu.draws<- as.data.frame(inf.object$draws(variables = "uconstrained")[,1,])
  }else{
    fullG12.draws<- as.numeric(inf.object[-(1:burn.in), 1])
    fullG21.draws<- as.numeric(inf.object[-(1:burn.in), 2])
    fullkappaR.draws<- as.numeric(inf.object[-(1:burn.in), 3])
    fullr.draws<- inf.object[-(1:burn.in), 5+(1:time)]
    fulls.draws<- inf.object[-(1:burn.in), 5+time+(1:12)]
    fullu.draws<- inf.object[-(1:burn.in), 5+time+12+(1:ndept)]
  }

  thinning<- numeric(floor(nrow(fullr.draws)/thinningL))
  thinning[1]<- thinningL
  for(i in 2:length(thinning)){
    thinning[i]<- thinning[i-1] + thinningL
  }

  G12.draws<- fullG12.draws[thinning]
  G21.draws<- fullG21.draws[thinning]
  kappaR.draws<- fullkappaR.draws[thinning]
  r.draws<- fullr.draws[thinning, ]
  s.draws<- fulls.draws[thinning, ]
  u.draws<- fullu.draws[thinning, ]

  Alldata<- matrix(0, nrow = ndept, ncol = forecastlength)
  Alldata.list<- vector("list", length(thinning))
  Outbreaks<- matrix(0, ndept, forecastlength)
  allRtforecasts<- matrix(NA, nrow = length(thinning), ncol = forecastlength)

  for(index in 1:length(thinning)){
  rt<- as.numeric(r.draws[index, ])
  s<- as.numeric(s.draws[index, ])
  u<- as.numeric(u.draws[index, ])
  kappa_r<- kappaR.draws[index]
  G12<- G12.draws[index]
  G21<- G21.draws[index]

  rt_forecast<- numeric(forecastlength)
  rt_forecast[1] <- 2 * rt[time-1] - rt[time-2] + rnorm(1, mean=0, sd=sqrt(1/kappa_r))
  rt_forecast[2] <- 2 * rt_forecast[1] - rt[time-1] + rnorm(1, mean=0, sd=sqrt(1/kappa_r))
  for(t in 3:forecastlength){
    rt_forecast[t] <- 2 * rt_forecast[t-1] - rt_forecast[t-2] + rnorm(1, mean=0, sd=sqrt(1/kappa_r))
  }

  e_it<- matrix(NA, nrow = ndept, ncol = forecastlength)
  for(i in 1:ndept){
    e_it[i, 1]<- popdata[i, time] + ((popdata[i, time]-popdata[i, time-1])/1) * 1
    e_it[i, 2]<- e_it[i, 1] + ((e_it[i, 1]-popdata[i, time])/1) * 1
    for(t in 3:forecastlength){
      e_it[i, t]<- e_it[i, t-1] + ((e_it[i, t-1]-e_it[i, t-2])/1) * 1
    }
  }

  y_it<- matrix(NA, ndept, forecastlength+1)
  y_it[,1]<- y[,time]

  EpidemicIndicator<- matrix(NA, ndept, forecastlength)
  transition_matrix <- matrix(c(1-G12, G12, G21, 1-G21), nrow = 2, byrow = TRUE)

  if(Model == 0){
    for(i in 1:ndept){
      for(t in 1:forecastlength){
        m<- (t - 1) %% 12 + 1
        lograte <- rt_forecast[t] + s[m] + u[i]
        y_it[i, t+1]<- rpois(1, lambda = e_it[i, t] * exp(lograte))
      }
    }
    Alldata<- Alldata + y_it[,-1]
    Alldata.list[[index]]<- y_it[,-1]
    allRtforecasts[index,]<- rt_forecast
    EpidemicIndicator<- matrix(0, ndept, forecastlength)
  }

  else if(Model == 1){
    if(!is.data.frame(inf.object)){
      B.draws<- stack(as.data.frame(inf.object$draws(variables = "B")[,1,]))[,1]
    }else{
      B.draws<- as.numeric(inf.object[-(1:burn.in), 5+time+12+ndept+1])
    }
    B.draws<- B.draws[thinning]
    B<- B.draws[index]

    for(i in 1:ndept) {
      states<- numeric(forecastlength)
      initial_state <- sample(0:1, 1, prob = c(1 - OutP[i], OutP[i]))
      states[1] <- initial_state

      for(t in 2:forecastlength){
        current_state <- states[t - 1]
        next_state <- sample(0:1, size = 1, prob = transition_matrix[current_state + 1, ])
        states[t] <- next_state
      }
      EpidemicIndicator[i,] <- states
    }

  for(t in 1:forecastlength){
      m<- (t - 1) %% 12 + 1
      for(i in 1:ndept){
        z_it<- ifelse(y_it[i, t]>0, 1, 0)
        lograte<- rt_forecast[t] + s[m] + u[i] + z_it * EpidemicIndicator[i, t] * B[1]
        y_it[i, t+1]<- rpois(1, lambda = e_it[i, t] * exp(lograte))
      }
    }
    Alldata<- Alldata + y_it[, -1]
    Alldata.list[[index]]<- y_it[,-1]
    allRtforecasts[index,]<- rt_forecast
    Outbreaks<- Outbreaks + EpidemicIndicator
  }

  else if(Model == 2){
    if(!is.data.frame(inf.object)){
      B.draws<- stack(as.data.frame(inf.object$draws(variables = "B")[,1,]))[,1]
    }else{
      B.draws<- as.numeric(inf.object[-(1:burn.in), 5+time+12+ndept+1])
    }
    B.draws<- B.draws[thinning]
    B<- B.draws[index]
    for(i in 1:ndept) {
      states<- numeric(forecastlength)
      initial_state <- sample(0:1, 1, prob = c(1 - OutP[i], OutP[i]))
      states[1] <- initial_state

      for(t in 2:forecastlength){
        current_state <- states[t - 1]
        next_state <- sample(0:1, size = 1, prob = transition_matrix[current_state + 1, ])
        states[t] <- next_state
      }
      EpidemicIndicator[i,] <- states
    }

    for(t in 1:forecastlength){
      m<- (t - 1) %% 12 + 1
      for(i in 1:ndept){
        indexes <- which(adj.matrix[i, ] > 0 & 1:ndept != i)
        z_it<-  ifelse(y_it[i, t]>0 || any(y_it[indexes, t]>0), 1, 0)
        lograte<- rt_forecast[t] + s[m] + u[i] + z_it * EpidemicIndicator[i, t] * B[1]
        y_it[i, t+1]<- rpois(1, lambda = e_it[i, t] * exp(lograte))
      }
    }
    Alldata<- Alldata + y_it[, -1]
    Alldata.list[[index]]<- y_it[,-1]
    allRtforecasts[index,]<- rt_forecast
    Outbreaks<- Outbreaks + EpidemicIndicator
  }

  else if(Model == 3){
    if(!is.data.frame(inf.object)){
      B.draws<- as.data.frame(inf.object$draws(variables = "B")[,1,])
    }else{
      B.draws<- inf.object[-(1:burn.in), 5+time+12+ndept+(1:2)]
    }
    B.draws<- B.draws[thinning, ]
    B<- as.numeric(B.draws[index,])

    for(i in 1:ndept) {
      states<- numeric(forecastlength)
      initial_state <- sample(0:1, 1, prob = c(1 - OutP[i], OutP[i]))
      states[1] <- initial_state

      for(t in 2:forecastlength){
        current_state <- states[t - 1]
        next_state <- sample(0:1, size = 1, prob = transition_matrix[current_state + 1, ])
        states[t] <- next_state
      }
      EpidemicIndicator[i,] <- states
    }

    for(t in 1:forecastlength){
      m<- (t - 1) %% 12 + 1
      for(i in 1:ndept){
        indexes <- which(adj.matrix[i, ] > 0 & 1:ndept != i)
        z_it1<- ifelse(y_it[i, t]>0, 1, 0)
        z_it2<- ifelse(any(y_it[indexes, t-1] > 0), 1, 0)
        lograte <- rt_forecast[t] + s[m] + u[i] + (z_it1 * EpidemicIndicator[i, t] * B[1]) + (z_it2 * EpidemicIndicator[i, t] * B[2])
        y_it[i, t+1]<- rpois(1, lambda = e_it[i, t] * exp(lograte))
      }
    }
    Alldata<- Alldata + y_it[, -1]
    Alldata.list[[index]]<- y_it[,-1]
    allRtforecasts[index,]<- rt_forecast
    Outbreaks<- Outbreaks + EpidemicIndicator
  }

  else if(Model == 4){
    if(!is.data.frame(inf.object)){
      B.draws<- stack(as.data.frame(inf.object$draws(variables = "B")[,1,]))[,1]
    }else{
      B.draws<- as.numeric(inf.object[-(1:burn.in), 5+time+12+ndept+1])
    }
    B.draws<- B.draws[thinning]
    B<- B.draws[index]
    for(i in 1:ndept) {
      states<- numeric(forecastlength)
      initial_state <- sample(0:1, 1, prob = c(1 - OutP[i], OutP[i]))
      states[1] <- initial_state

      for(t in 2:forecastlength){
        current_state <- states[t - 1]
        next_state <- sample(0:1, size = 1, prob = transition_matrix[current_state + 1, ])
        states[t] <- next_state
      }
      EpidemicIndicator[i,] <- states
    }

    for(t in 1:forecastlength){
      m<- (t - 1) %% 12 + 1
      for(i in 1:ndept){
        z_it<- log(y_it[i, t] + 1)
        lograte<- rt_forecast[t] + s[m] + u[i] + z_it * EpidemicIndicator[i, t] * B[1]
        y_it[i, t+1]<- rpois(1, lambda = e_it[i, t] * exp(lograte))
      }
    }
    Alldata<- Alldata + y_it[, -1]
    Alldata.list[[index]]<- y_it[,-1]
    allRtforecasts[index,]<- rt_forecast
    Outbreaks<- Outbreaks + EpidemicIndicator
  }

  else if(Model == 5){
    if(!is.data.frame(inf.object)){
      B.draws<- stack(as.data.frame(inf.object$draws(variables = "B")[,1,]))[,1]
    }else{
      B.draws<- as.numeric(inf.object[-(1:burn.in), 5+time+12+ndept+1])
    }
    B.draws<- B.draws[thinning]
    B<- B.draws[index]
    for(i in 1:ndept) {
      states<- numeric(forecastlength)
      initial_state <- sample(0:1, 1, prob = c(1 - OutP[i], OutP[i]))
      states[1] <- initial_state

      for(t in 2:forecastlength){
        current_state <- states[t - 1]
        next_state <- sample(0:1, size = 1, prob = transition_matrix[current_state + 1, ])
        states[t] <- next_state
      }
      EpidemicIndicator[i,] <- states
    }

    for(t in 1:forecastlength){
      m<- (t - 1) %% 12 + 1
      for(i in 1:ndept){
        indexes <- which(adj.matrix[i, ] > 0 & 1:ndept != i)
        z_it<- log(y_it[i, t] + sum(y_it[indexes, t]) + 1)
        lograte<- rt_forecast[t] + s[m] + u[i] + z_it * EpidemicIndicator[i, t] * B[1]
        y_it[i, t+1]<- rpois(1, lambda = e_it[i, t] * exp(lograte))
      }
    }
    Alldata<- Alldata + y_it[, -1]
    Alldata.list[[index]]<- y_it[,-1]
    allRtforecasts[index,]<- rt_forecast
    Outbreaks<- Outbreaks + EpidemicIndicator
  }

  else if(Model == 6){
    if(!is.data.frame(inf.object)){
      B.draws<- as.data.frame(inf.object$draws(variables = "B")[,1,])
    }else{
      B.draws<- inf.object[-(1:burn.in), 5+time+12+ndept+(1:2)]
    }
    B.draws<- B.draws[thinning, ]
    B<- as.numeric(B.draws[index,])

    for(i in 1:ndept) {
      states<- numeric(forecastlength)
      initial_state <- sample(0:1, 1, prob = c(1 - OutP[i], OutP[i]))
      states[1] <- initial_state

      for(t in 2:forecastlength){
        current_state <- states[t - 1]
        next_state <- sample(0:1, size = 1, prob = transition_matrix[current_state + 1, ])
        states[t] <- next_state
      }
      EpidemicIndicator[i,] <- states
    }

    for(t in 1:forecastlength){
      m<- (t - 1) %% 12 + 1
      for(i in 1:ndept){
        indexes <- which(adj.matrix[i, ] > 0 & 1:ndept != i)
        z_it1<- log(y_it[i, t] + 1)
        z_it2<- log(sum(y_it[indexes, t-1]) + 1)
        lograte <- rt_forecast[t] + s[m] + u[i] + (z_it1 * EpidemicIndicator[i, t] * B[1]) + (z_it2 * EpidemicIndicator[i, t] * B[2])
        y_it[i, t+1]<- rpois(1, lambda = e_it[i, t] * exp(lograte))
      }
    }
    Alldata<- Alldata + y_it[, -1]
    Alldata.list[[index]]<- y_it[,-1]
    allRtforecasts[index,]<- rt_forecast
    Outbreaks<- Outbreaks + EpidemicIndicator
  }

  else if(Model == 7){
    if(!is.data.frame(inf.object)){
      B.draws<- stack(as.data.frame(inf.object$draws(variables = "B")[,1,]))[,1]
    }else{
      B.draws<- as.numeric(inf.object[-(1:burn.in), 5+time+12+ndept+1])
    }
    B.draws<- B.draws[thinning]
    B<- B.draws[index]
    for(i in 1:ndept) {
      states<- numeric(forecastlength)
      initial_state <- sample(0:1, 1, prob = c(1 - OutP[i], OutP[i]))
      states[1] <- initial_state

      for(t in 2:forecastlength){
        current_state <- states[t - 1]
        next_state <- sample(0:1, size = 1, prob = transition_matrix[current_state + 1, ])
        states[t] <- next_state
      }
      EpidemicIndicator[i,] <- states
    }

    for(t in 1:forecastlength){
      m<- (t - 1) %% 12 + 1
      for(i in 1:ndept){
        z_it<- 1
        lograte<- rt_forecast[t] + s[m] + u[i] + z_it * EpidemicIndicator[i, t] * B[1]
        y_it[i, t+1]<- rpois(1, lambda = e_it[i, t] * exp(lograte))
      }
    }
    Alldata<- Alldata + y_it[, -1]
    Alldata.list[[index]]<- y_it[,-1]
    allRtforecasts[index,]<- rt_forecast
    Outbreaks<- Outbreaks + EpidemicIndicator
    }
  }
  return(list("Counts" = Alldata/length(thinning), "Population" = e_it, "TrendComp" = allRtforecasts, "SeasonalComps" = s.draws, "SpatialComps" = colMeans(fullu.draws), "Outbreaks" = Outbreaks/length(thinning), "FullCounts" = Alldata.list))
}

forecastplot<- function(forecast.object, countries=NULL){
  Alldata.list<- forecast.object$FullCounts
  Trendcomps<- forecast.object$TrendComp
  Seasonalcomps<- forecast.object$SeasonalComps
  Outbreaks<- forecast.object$Outbreaks
  samples<- length(Alldata.list)
  ndept<- nrow(Alldata.list[[1]])
  time<- ncol(Alldata.list[[1]])

  forecastsamples<- matrix(NA, nrow = samples, ncol = time)
  sortedforecasts<- vector("list", length(ndept))
  for(i in 1:ndept){
    for(t in 1:samples){
      forecastsamples[t, ]<- Alldata.list[[t]][i,]
    }
    sortedforecasts[[i]]<- forecastsamples
  }
  #ind<- seq(1,ndept-2,by=2)
  par(mfrow=c(3,3))
  for(i in 1:ndept){
    #Counts
    mean.Y<- colMeans(sortedforecasts[[i]])
    uCI.Y<- posterior_interval_custom(as.matrix.data.frame(sortedforecasts[[i]]))[,2]
    lCI.Y<- posterior_interval_custom(as.matrix.data.frame(sortedforecasts[[i]]))[,1]
    plot(0, type = "n", xlim = c(1,time), ylim = c(min(sortedforecasts[[i]], lCI.Y), max(sortedforecasts[[i]], uCI.Y)), ylab = "case counts", xlab = "Time [month]", cex.axis = 1.8, cex.lab=1.8)
    polygon(c(1:time, rev(1:time)), c(lCI.Y, rev(uCI.Y)),
            col = "pink", border = NA)
    lines(1:time, mean.Y, col = "red", lty=1)

    #mean.Y<- colMeans(sortedforecasts[[i+1]])
    #uCI.Y<- posterior_interval_custom(as.matrix.data.frame(sortedforecasts[[i+1]]))[,2]
    #lCI.Y<- posterior_interval_custom(as.matrix.data.frame(sortedforecasts[[i+1]]))[,1]
    #polygon(c(1:time, rev(1:time)), c(lCI.Y, rev(uCI.Y)),
    #        col = "pink", border = NA)
    #lines(1:time, mean.Y, col = "red", lty=1)

    #mean.Y<- colMeans(sortedforecasts[[i+2]])
    #uCI.Y<- posterior_interval_custom(as.matrix.data.frame(sortedforecasts[[i+2]]))[,2]
    #lCI.Y<- posterior_interval_custom(as.matrix.data.frame(sortedforecasts[[i+2]]))[,1]
    #polygon(c(1:time, rev(1:time)), c(lCI.Y, rev(uCI.Y)),
    #        col = "pink", border = NA)
    #lines(1:time, mean.Y, col = "red", lty=1)
    #points(1:ncol(y), colSums(y, na.rm = T), pch = 19)
    grid()
  }
  par(mfrow=c(1,3))
  #Trend
  mean.r<- colMeans(Trendcomps)
  uCI.r<- posterior_interval_custom(as.matrix.data.frame(Trendcomps))[,2]
  lCI.r<- posterior_interval_custom(as.matrix.data.frame(Trendcomps))[,1]
  #par(mar = c(4, 7.5, 4, 1))
  plot(0, type = "n", xlim = c(1,time), ylim = c(min(Trendcomps, lCI.r), max(Trendcomps, uCI.r)), ylab = "trend component", xlab = "Time [month]", cex.axis = 1.8, cex.lab=1.8)
  polygon(c(1:time, rev(1:time)), c(lCI.r, rev(uCI.r)),
          col = "pink", border = NA)
  lines(1:time, mean.r, col = "red", lty=1)
  #points(1:ncol(y), colSums(y, na.rm = T), pch = 19)
  grid()

  #Seasonal
  mean.s<- colMeans(Seasonalcomps)
  uCI.s<- posterior_interval_custom(as.matrix.data.frame(Seasonalcomps))[,2]
  lCI.s<- posterior_interval_custom(as.matrix.data.frame(Seasonalcomps))[,1]
  #par(mar = c(4, 7.5, 4, 1))
  plot(0, type = "n", xlim = c(1,12), ylim = c(min(Seasonalcomps, lCI.s), max(Seasonalcomps, uCI.s)), ylab = "seasonal component", xlab = "Time [month]", cex.axis = 1.8, cex.lab=1.8)
  polygon(c(1:12, rev(1:12)), c(lCI.s, rev(uCI.s)),
          col = "pink", border = NA)
  lines(1:12, mean.s, col = "red", lty=1)
  #points(1:ncol(y), colSums(y, na.rm = T), pch = 19)
  grid()

  #Outbreak
  par(mar = c(4, 7.5, 4, 1))
  Outbreaks<- Outbreaks[ndept:1, ]
  image(x=1:time, y=1:ndept, t(Outbreaks), axes=F, col = hcl.colors(50, "YlOrRd", rev = TRUE), main ="", ylab="", xlab="Time [month/year]", cex.lab=1.80)
  #custom Y-axis
  axis(2, at=seq(1, length(countries), length.out=length(countries)), labels=rev(countries), lwd.ticks = 1, las = 1, lwd=0, cex.axis = 1.40)
  #custom X-axis
  years<- 2020:2021
  axis(1, at = seq(1, 24, by = 12), labels = years, cex.axis = 1.8)
}

forecastplot2<- function(forecast.object){
  Alldata.list<- forecast.object$FullCounts
  Trendcomps<- forecast.object$TrendComp
  Seasonalcomps<- forecast.object$SeasonalComps
  Outbreaks<- forecast.object$Outbreaks
  samples<- length(Alldata.list)
  ndept<- nrow(Alldata.list[[1]])
  time<- ncol(Alldata.list[[1]])
  forecastsamples<- matrix(NA, nrow = samples, ncol = time)
  sortedforecasts<- vector("list", length(ndept))
  for(i in 1:ndept){
    for(t in 1:samples){
      forecastsamples[t, ]<- Alldata.list[[t]][i,]
    }
    sortedforecasts[[i]]<- forecastsamples
  }
    uCI.Y<- c()
    lCI.Y<- c()
  for(i in 1:ndept){
    uCI.Y<- c(uCI.Y, posterior_interval_custom(as.matrix.data.frame(sortedforecasts[[i]]))[,2])
    lCI.Y<- c(lCI.Y, posterior_interval_custom(as.matrix.data.frame(sortedforecasts[[i]]))[,1])
  }
    ci_data<- data.frame(Time=1:time, variable=rep(paste("u", 1:ndept, sep = ""), each=time), lCI.Y=lCI.Y, uCI.Y=uCI.Y)

  if (is.list(forecast.object)) {
    spatdata <- forecast.object[[1]]
  } else {
    spatdata <- forecast.object
  }

  ts_spatdata <- as.data.frame(t(spatdata))
  ts_spatdata$Time <- 1:ncol(spatdata)
  colnames(ts_spatdata) <- c(paste("u", 1:(ncol(ts_spatdata) - 1), sep = ""), "Time")
  long_data <- reshape2::melt(ts_spatdata, id.vars = "Time")

  library(ggplot2)

  p <- ggplot(data = long_data, aes(x = Time, y = value, color = variable, fill = variable)) +
    geom_line() +
    labs(x = "Time [month/year]", y = "Case counts", color = "Location", fill = "Location") +
    #guides(color = guide_legend("Location"), linetype = guide_legend("Location")) +
    theme(axis.title.y = element_text(size = 18),
          axis.title.x = element_text(size = 18),
          axis.text.x = element_text(size = 16),
          axis.text.y = element_text(size = 16),
          legend.title = element_text(size = 18),
          legend.text = element_text(size = 16))

    p <- p + geom_ribbon(
      data = ci_data,
      aes(x = Time, ymin = lCI.Y, ymax = uCI.Y, fill = variable, group = variable),
      alpha = 0.2,
      inherit.aes = FALSE
    )
  return(p)
}
