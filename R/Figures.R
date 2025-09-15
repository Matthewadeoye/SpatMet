Posteriorpredictive<- function(all.infobjects,  realdata, adjmat){
  time<- ncol(realdata[[1]])
  ndept<- nrow(realdata[[1]])

  pdf("PosteriorpredictiveOutbreak.pdf", paper="special", width=18,height=9, pointsize=12)
  par(mfrow=c(2,4))

  for(i in 1:8){
    inf.object<- all.infobjects[[i]]
    Model<- i-1

    y<- realdata[[1]]
    e_it<- realdata[[2]]

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

    sum_Y<- matrix(NA, nrow = length(thinning), ncol = time)
    sum_Xit<- matrix(0, nrow = ndept, ncol = time)

    if(Model == 0){
      for(index in 1:length(thinning)){
        for(i in 1:ndept){
          for(t in 1:time){
            m<- (t - 1) %% 12 + 1
            Exlambda_it <- e_it[i, t] * exp(r.draws[index, t] + s.draws[index, m] + u.draws[index, i])
            pred_Y[i, t]<- rpois(1, Exlambda_it)
          }
        }
        sum_Y[index, ]<- colSums(pred_Y)
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
        sum_Xit<- sum_Xit + Ex_Xit
        for(i in 1:ndept){
          for(t in 1:time){
            m<- (t - 1) %% 12 + 1
            P_Xit<- rbinom(1, 1, prob = Ex_Xit[i, t])
            Exlambda_it <- e_it[i, t] * exp(r.draws[index, t] + s.draws[index, m] + u.draws[index, i] + P_Xit * z_it[i, t] * B.draws[index])
            pred_Y[i, t]<- rpois(1, Exlambda_it)
          }
        }
        sum_Y[index, ]<- colSums(pred_Y)
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
        sum_Y[index, ]<- colSums(pred_Y)
      }
    }
    inf.r<- colMeans(fullr.draws)
    uCI.r<- posterior_interval_custom(as.matrix.data.frame(fullr.draws))[,2]
    lCI.r<- posterior_interval_custom(as.matrix.data.frame(fullr.draws))[,1]
    inf.s<- colMeans(fulls.draws)
    uCI.s<- posterior_interval_custom(as.matrix.data.frame(fulls.draws))[,2]
    lCI.s<- posterior_interval_custom(as.matrix.data.frame(fulls.draws))[,1]
    inf.u<- colMeans(fullu.draws)
    uCI.u<- posterior_interval_custom(as.matrix.data.frame(fullu.draws))[,2]
    lCI.u<- posterior_interval_custom(as.matrix.data.frame(fullu.draws))[,1]
    inf.Y<- colMeans(sum_Y)
    uCI.Y<- posterior_interval_custom(as.matrix.data.frame(sum_Y))[,2]
    lCI.Y<- posterior_interval_custom(as.matrix.data.frame(sum_Y))[,1]

    plot(0, type = "n", xlim = c(1,ncol(y)), ylim = c(min(lCI.Y, y, na.rm = T), max(uCI.Y, y, na.rm = T)), ylab = "overall case counts", xlab = "Time [month]", cex.axis = 1.7, cex.lab=1.8)
    polygon(c(1:length(inf.Y), rev(1:length(inf.Y))), c(lCI.Y, rev(uCI.Y)),
            col = "pink", border = NA)
    lines(1:length(inf.Y), inf.Y, col = "red", lty=1)
    points(1:ncol(y), colSums(y, na.rm = T), pch = 19)
    grid()
    legendary::labelFig(LETTERS[Model+1], adj = c(-0.15, 0.10), font=2, cex=1.8)
  }
  add_legend("topright", legend=c("Truth", "Posterior means"), lty=c(NA, 1),
             pch=c(19, NA), col=c("black", "red"),
             horiz=TRUE, bty='n', cex=1.8)
  dev.off()
}

Temporalcomps<- function(all.infobjects){
  rmat<- matrix(NA, nrow = 8, ncol = 60)
  smat<- matrix(NA, nrow = 8, ncol = 12)
  for(i in 1:8){
    inf.object<- all.infobjects[[i]]
    Model<- i

    if(!is.data.frame(inf.object)){
      r.draws<- colMeans(as.data.frame(inf.object$draws(variables = "r")[,1,]))
      s.draws<- colMeans(as.data.frame(inf.object$draws(variables = "s")[,1,]))
    }else{
      r.draws<- colMeans(inf.object[-(1:burn.in), 5+(1:time)])
      s.draws<- colMeans(inf.object[-(1:burn.in), 5+time+(1:12)])
    }
    rmat[i, ]<- r.draws
    smat[i, ]<- s.draws
  }
  r_names<- c("0", "I", "II", "III", "IV", "V", "VI", "VII")
  ts_rdata <- as.data.frame(t(rmat))
  ts_sdata <- as.data.frame(t(smat))
  ts_rdata$Time <- 1:nrow(ts_rdata)
  ts_sdata$Time <- 1:12

  colnames(ts_rdata) <- c(r_names, "Time")
  colnames(ts_sdata) <- c(r_names, "Time")

  long_data <- reshape2::melt(ts_rdata, id.vars = "Time")
  library(ggplot2)
  a<- (ggplot2::ggplot(data = long_data, mapping = aes(x = Time, y = value, color = variable)) +
         geom_line() +
         labs(x = "Time [month]", y = "Trend component", color = "Model") +
         theme(axis.title.y = element_text(size=18),
               axis.title.x = element_text(size=18),
               axis.text.x = element_text(size=16),
               axis.text.y = element_text(size=16),
               legend.title = element_text(size = 18),
               legend.text = element_text(size = 16), legend.position = "none"))

  long_data2 <- reshape2::melt(ts_sdata, id.vars = "Time")
  library(ggplot2)
  b<- (ggplot2::ggplot(data = long_data2, mapping = aes(x = Time, y = value, color = variable)) +
         geom_line() +
         labs(x = "Time [month]", y = "Seasonal component", color = "Model") +
         theme(axis.title.y = element_text(size=18),
               axis.title.x = element_text(size=18),
               axis.text.x = element_text(size=16),
               axis.text.y = element_text(size=16),
               legend.title = element_text(size = 18),
               legend.text = element_text(size = 16)))
  plotlists<- list(a, b)
  print(cowplot::plot_grid(plotlist = plotlists, ncol = 2, labels = c("A", "B"), label_size = 17, rel_widths = c(1.08, 1.15)))
}



AllPosteriorpredictive<- function(all.infobjects,  realdata, adjmat){
  time<- ncol(realdata[[1]])
  ndept<- nrow(realdata[[1]])

  pdf("AllPosteriorpredictive.pdf", paper="special", width=28,height=24, pointsize=12)
  par(mfrow=c(8,9))

  for(m in 1:8){
    inf.object<- all.infobjects[[m]]
    Model<- m-1

    y<- realdata[[1]]
    e_it<- realdata[[2]]
    for(k in 1:ndept){

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

    sum_Y<- matrix(NA, nrow = length(thinning), ncol = time)
    sum_Xit<- matrix(0, nrow = ndept, ncol = time)

    if(Model == 0){
      for(index in 1:length(thinning)){
        for(i in 1:ndept){
          for(t in 1:time){
            m<- (t - 1) %% 12 + 1
            Exlambda_it <- e_it[i, t] * exp(r.draws[index, t] + s.draws[index, m] + u.draws[index, i])
            pred_Y[i, t]<- rpois(1, Exlambda_it)
          }
        }
        sum_Y[index, ]<- pred_Y[k, ]
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
        sum_Xit<- sum_Xit + Ex_Xit
        for(i in 1:ndept){
          for(t in 1:time){
            m<- (t - 1) %% 12 + 1
            P_Xit<- rbinom(1, 1, prob = Ex_Xit[i, t])
            Exlambda_it <- e_it[i, t] * exp(r.draws[index, t] + s.draws[index, m] + u.draws[index, i] + P_Xit * z_it[i, t] * B.draws[index])
            pred_Y[i, t]<- rpois(1, Exlambda_it)
          }
        }
        sum_Y[index, ]<- pred_Y[k,]
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
        sum_Y[index, ]<- pred_Y[k, ]
      }
    }

    inf.Y<- colMeans(sum_Y)
    uCI.Y<- posterior_interval_custom(as.matrix.data.frame(sum_Y))[,2]
    lCI.Y<- posterior_interval_custom(as.matrix.data.frame(sum_Y))[,1]

    plot(0, type = "n", xlim = c(1,ncol(y)), ylim = c(min(lCI.Y, y, na.rm = T), max(uCI.Y, y, na.rm = T)), main = paste0("Model = ", Model, ",", " Location = ", k), ylab = "overall case counts", xlab = "Time [month]", cex.axis = 1.7, cex.lab=1.8, cex.main=2.0)
    polygon(c(1:length(inf.Y), rev(1:length(inf.Y))), c(lCI.Y, rev(uCI.Y)),
            col = "pink", border = NA)
    lines(1:length(inf.Y), inf.Y, col = "red", lty=1)
    points(1:ncol(y), y[k, ], pch = 19)
    grid()
    }
  }
#  add_legend("topright", legend=c("Truth", "Posterior means"), lty=c(NA, 1),
#             pch=c(19, NA), col=c("black", "red"),
#             horiz=TRUE, bty='n', cex=1.8)
  dev.off()
}

Heatmaps<- function(all.infobjects, all.data, adjmat){
  y<- all.data[[1]]
  e_it<- all.data[[2]]
  time<- ncol(y)
  ndept<- nrow(y)

  #pdf("HeatmapsOutbreak.pdf", paper="special", width=21,height=12, pointsize=14)
  par(mfrow=c(2,4))

  for(i in 1:7){
    inf.object<- all.infobjects[[i+1]]
    Model<- i

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

    sum_Xit<- matrix(0, nrow = ndept, ncol = time)

    if(Model %in% c(1,2,4,5,7)){
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
        sum_Xit<- sum_Xit + Ex_Xit
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
      }
    }

    if(Model != 0){
      mean_Xit<- sum_Xit/length(thinning)
      smallmeanxit<- mean_Xit[c(1,3,5,7,9), ]
      bigmeanxit<- mean_Xit[c(2,4,6,8), ]
      bigsmallmeanxit<- mean_Xit[c(2,4,6,8,1,3,5,7,9), ]
      image(x=1:time, y=1:ndept, t(bigsmallmeanxit), main = "", axes=F, ylab = "spatial location", xlab = "Time [month]", cex.lab=1.80)
      abline(h=4.5, col="black", lty=2)
      #custom Y-axis
      axis(2, at=seq(1, 4, length.out=4), labels=c("u2", "u4", "u6", "u8"), col = "red", col.axis="red", lwd.ticks = 1, las = 1, lwd=0, cex.axis = 1.8, cex.lab=1.8)
      axis(2, at=seq(5, 9, length.out=5), labels=c("u1", "u3", "u5", "u7", "u9"), col = "blue", col.axis="blue", lwd.ticks = 1, las = 1, lwd=0, cex.axis = 1.8, cex.lab=1.8)
      #custom X-axis
      axis(1, cex.axis = 1.8)
      legendary::labelFig(LETTERS[Model], adj = c(-0.15, 0.10), font=2, cex=1.8)
    }
  }
  #dev.off()
}


multstraindatafig2<- function(y, maxAll, Modeltype = ""){
  nstrain<- dim(y)[3]
  maxY<- maxAll
  plotlists<- list()

  for(i in 1:nstrain){
    Strain<- i
    sim.object<- y[,,i]

    spatdata<- sim.object

    ts_spatdata <- as.data.frame(t(spatdata))
    ts_spatdata$Time <- 1:ncol(spatdata)
    naming<- c(paste("u", 1:(ncol(ts_spatdata)-1), sep=""), "Time")
    colnames(ts_spatdata)<- naming

    Colors <- rep(c("blue", "red"), length.out = nrow(spatdata))
    Linetypes <- rep(c("dotted", "dashed", "dotdash", "longdash", "twodash"), length.out = nrow(spatdata))

    long_data <- reshape2::melt(ts_spatdata, id.vars = "Time")
    library(ggplot2)
    if(Strain==nstrain){
      rfigs<- ggplot2::ggplot(data = long_data, mapping = aes(x = Time, y = value, color = variable, linetype = variable)) +
        geom_line() +
        scale_color_manual(values = Colors) +
        scale_linetype_manual(values = Linetypes) +
        ylim(0, maxY) +
        labs(x = "Time [month]", y = "Case counts", color = "Location") +
        guides(color = guide_legend("Location"), linetype = guide_legend("Location")) +
        theme(axis.title.y = element_text(size=18),
              axis.title.x = element_text(size=18),
              axis.text.x = element_text(size=16),
              axis.text.y = element_text(size=16),
              legend.title = element_text(size = 18),
              legend.text = element_text(size = 16))
    }else{
      rfigs<- ggplot2::ggplot(data = long_data, mapping = aes(x = Time, y = value, color = variable, linetype = variable)) +
        geom_line() +
        scale_color_manual(values = Colors) +
        scale_linetype_manual(values = Linetypes) +
        ylim(0, maxY) +
        labs(x = "Time [month]", y = "Case counts", color = "Location") +
        guides(color = guide_legend("Location"), linetype = guide_legend("Location")) +
        theme(axis.title.y = element_text(size=18),
              axis.title.x = element_text(size=18),
              axis.text.x = element_text(size=16),
              axis.text.y = element_text(size=16),
              legend.position = "none")
    }
    plotlists[[Strain]]<- rfigs
  }
  row_1<- cowplot::plot_grid(plotlist = plotlists[1:3], ncol = 3, labels = c("A", "B", "C"), label_size = 17)
  row_2<- cowplot::plot_grid(plotlist = plotlists[4:5], ncol = 3, labels = c("D", "E"), label_size = 17, rel_widths = c(1, 1.35, 0.65))
  print(cowplot::plot_grid(row_1, row_2, nrow = 2))
  add_legend(0.6, -0.3, legend=substitute(paste(bold(Modeltype))),
              col="black",
             horiz=TRUE, bty='n', cex=3.0)
}


RecoverInf.plot<- function(inf.object, true_r, true_s, Modeltype=""){

  if(is.data.frame(inf.object)){
    rPosterior<- inf.object[, startsWith(colnames(inf.object), "r")]
    sPosterior<- inf.object[, startsWith(colnames(inf.object), "s")]
    inf.r<- colMeans(rPosterior)
    inf.s<- colMeans(sPosterior)
    uCI.r<- posterior_interval_custom(as.matrix.data.frame(rPosterior))[,2]
    lCI.r<- posterior_interval_custom(as.matrix.data.frame(rPosterior))[,1]
    uCI.s<- posterior_interval_custom(as.matrix.data.frame(sPosterior))[,2]
    lCI.s<- posterior_interval_custom(as.matrix.data.frame(sPosterior))[,1]
  }else{
    inf.r<- colMeans(as.data.frame(inf.object$draws(variables = "r")[,1,]))
    uCI.r<- posterior_interval_custom(as.matrix.data.frame(inf.object$draws(variables = "r")[,1,]))[,2]
    lCI.r<- posterior_interval_custom(as.matrix.data.frame(inf.object$draws(variables = "r")[,1,]))[,1]
    inf.s<- colMeans(as.data.frame(inf.object$draws(variables = "s")[,1,]))
    uCI.s<- posterior_interval_custom(as.matrix.data.frame(inf.object$draws(variables = "s")[,1,]))[,2]
    lCI.s<- posterior_interval_custom(as.matrix.data.frame(inf.object$draws(variables = "s")[,1,]))[,1]
    inf.u<- colMeans(as.data.frame(inf.object$draws(variables = "uconstrained")[,1,]))
    uCI.u<- posterior_interval_custom(as.matrix.data.frame(inf.object$draws(variables = "uconstrained")[,1,]))[,2]
    lCI.u<- posterior_interval_custom(as.matrix.data.frame(inf.object$draws(variables = "uconstrained")[,1,]))[,1]
  }

  par(mfrow=c(1,2))
  plot(0, type = "n", xlim = c(1,length(inf.r)), ylim = c(min(lCI.r, inf.r), max(uCI.r, inf.r)), ylab = "Trend component", xlab = "Time")
  polygon(c(1:length(inf.r), rev(1:length(inf.r))), c(lCI.r, rev(uCI.r)),
          col = "pink", border = NA)
  lines(1:length(inf.r), inf.r, col="red")
  points(1:length(inf.r), true_r, pch = 19)
  grid()

  plot(0, type = "n", xlim = c(1,length(inf.s)), ylim = c(min(lCI.s, inf.s), max(uCI.s, inf.s)), ylab = "Seasonal component", xlab = "Season")
  polygon(c(1:length(inf.s), rev(1:length(inf.s))), c(lCI.s, rev(uCI.s)),
          col = "pink", border = NA)
  lines(1:length(inf.s), inf.s, col="red")
  points(1:length(inf.s), true_s, pch = 19)
  grid()
  add_legend("topright", legend=c("Truth", "Posterior means"), lty=c(NA, 1),
             pch=c(19, NA), col=c("black", "red"),
             horiz=TRUE, bty='n', cex=1.4)
  add_legend("topleft", legend=substitute(paste(bold(Modeltype))),
             col=c("black", "red"),
             horiz=TRUE, bty='n', cex=1.5)

}

RecoverInfU.plot<- function(inf.object, true_u, burn.in=100){

  if(!is.data.frame(inf.object)){
    fullu.draws<- as.data.frame(inf.object$draws(variables = "uconstrained")[,1,])
  }else{
    fullu.draws<- inf.object[-(1:burn.in), 5+time+12+(1:ndept)]
  }

  thinning<- numeric(floor(nrow(fullu.draws)/10))
  thinning[1]<- 10
  for(i in 2:length(thinning)){
    thinning[i]<- thinning[i-1] + 10
  }
  u.draws<- fullu.draws[thinning, ]
  sim.u<- true_u

  # Violin plot for spatial components and intercepts
  spatcomp<- data.frame(value = c(u.draws[, 1], u.draws[, 2], u.draws[, 3], u.draws[, 4],
                                  u.draws[, 5], u.draws[, 6], u.draws[, 7], u.draws[, 8],
                                  u.draws[, 9]), group = factor(rep(c("u1", "u2", "u3", "u4", "u5", "u6", "u7", "u8", "u9"), each = nrow(u.draws))))
  spatcomp$group <- factor(spatcomp$group, levels = unique(spatcomp$group))
  library(ggplot2)
  library(RColorBrewer)
  mycolors <- c(rep(c("blue", "red"), 4), "blue")
  rfigs<- ggplot(spatcomp, aes(x = group, y = value, fill = group)) +
    geom_violin(trim = FALSE, alpha = 0.7) +
    geom_point(x = 1, y = sim.u[1], size = 2, shape = 19) +
    geom_point(x = 2, y = sim.u[2], size = 2, shape = 19) +
    geom_point(x = 3, y = sim.u[3], size = 2, shape = 19) +
    geom_point(x = 4, y = sim.u[4], size = 2, shape = 19) +
    geom_point(x = 5, y = sim.u[5], size = 2, shape = 19) +
    geom_point(x = 6, y = sim.u[6], size = 2, shape = 19) +
    geom_point(x = 7, y = sim.u[7], size = 2, shape = 19) +
    geom_point(x = 8, y = sim.u[8], size = 2, shape = 19) +
    geom_point(x = 9, y = sim.u[9], size = 2, shape = 19) +
    #ylim(-0.90, 0.70) +
    labs(title = "", x = "Location", y = "Value", fill = "") +
    theme_minimal() +
    scale_fill_manual(values = mycolors) +
    theme(axis.title.y = element_text(size=18),
          axis.title.x = element_text(size=18),
          axis.text.x = element_text(size=16),
          axis.text.y = element_text(size=16),
          legend.title = element_text(size = 18),
          legend.text = element_text(size = 16),legend.position = "none")

  print(rfigs)
  add_legend(0.85, 1.15, legend="Truth",
             pch=19, col="black",
             horiz=TRUE, bty='n', cex=1.8)
  #add_legend("topright", legend="Truth",
  #           pch=19, col="black",
  #           horiz=TRUE, bty='n', cex=1.1)
}

RecoverInfU_ak.plot<- function(inf.object, true_u, true_a_k, Modeltype="", burn.in=100){

  nstrain<- length(true_a_k)
    if(!is.data.frame(inf.object)){
      fullu.draws<- as.data.frame(inf.object$draws(variables = "uconstrained")[,1,])
      fullak.draws<- as.data.frame(inf.object$draws(variables = "a_k")[,1,])
    }else{
      fullu.draws<- inf.object[-(1:burn.in), 5+time+12+(1:ndept)]
      fullak.draws<- inf.object[-(1:burn.in), 5+time+12+ndept+nstrain+(1:nstrain)]
    }

    thinning<- numeric(floor(nrow(fullu.draws)/10))
    thinning[1]<- 10
    for(i in 2:length(thinning)){
      thinning[i]<- thinning[i-1] + 10
    }
    u.draws<- fullu.draws[thinning, ]
    ak.draws<- fullak.draws[thinning, ]
    sim.u<- true_u
    sim.ak<- true_a_k

    # Violin plot for spatial components and intercepts
    spatcomp<- data.frame(value = c(u.draws[, 1], u.draws[, 2], u.draws[, 3], u.draws[, 4],
                                    u.draws[, 5], u.draws[, 6], u.draws[, 7], u.draws[, 8],
                                    u.draws[, 9]), group = factor(rep(c("u1", "u2", "u3", "u4", "u5", "u6", "u7", "u8", "u9"), each = nrow(u.draws))))
    spatcomp$group <- factor(spatcomp$group, levels = unique(spatcomp$group))
    library(ggplot2)
    library(RColorBrewer)
    mycolors <- c(rep(c("blue", "red"), 4), "blue")
    rfigs<- ggplot(spatcomp, aes(x = group, y = value, fill = group)) +
      geom_violin(trim = FALSE, alpha = 0.7) +
      geom_point(x = 1, y = sim.u[1], size = 2, shape = 19) +
      geom_point(x = 2, y = sim.u[2], size = 2, shape = 19) +
      geom_point(x = 3, y = sim.u[3], size = 2, shape = 19) +
      geom_point(x = 4, y = sim.u[4], size = 2, shape = 19) +
      geom_point(x = 5, y = sim.u[5], size = 2, shape = 19) +
      geom_point(x = 6, y = sim.u[6], size = 2, shape = 19) +
      geom_point(x = 7, y = sim.u[7], size = 2, shape = 19) +
      geom_point(x = 8, y = sim.u[8], size = 2, shape = 19) +
      geom_point(x = 9, y = sim.u[9], size = 2, shape = 19) +
      ylim(-0.90, 0.70) +
      labs(title = "", x = "Location", y = "Value", fill = "") +
      theme_minimal() +
      scale_fill_manual(values = mycolors) +
      theme(axis.title.y = element_text(size=18),
            axis.title.x = element_text(size=18),
            axis.text.x = element_text(size=16),
            axis.text.y = element_text(size=16),
            legend.title = element_text(size = 18),
            legend.text = element_text(size = 16),legend.position = "none")

    spatcomp<- data.frame(value = c(ak.draws[, 1], ak.draws[, 2], ak.draws[, 3], ak.draws[, 4],
                                    ak.draws[, 5]), group = factor(rep(c("a1", "a2", "a3", "a4", "a5"), each = nrow(ak.draws))))
    spatcomp$group <- factor(spatcomp$group, levels = unique(spatcomp$group))
    mycolors <- c(rep("purple", 5))
    rfigs2<- ggplot(spatcomp, aes(x = group, y = value, fill = group)) +
      geom_violin(trim = FALSE, alpha = 0.7) +
      geom_point(x = 1, y = sim.ak[1], size = 2, shape = 19) +
      geom_point(x = 2, y = sim.ak[2], size = 2, shape = 19) +
      geom_point(x = 3, y = sim.ak[3], size = 2, shape = 19) +
      geom_point(x = 4, y = sim.ak[4], size = 2, shape = 19) +
      geom_point(x = 5, y = sim.ak[5], size = 2, shape = 19) +
      #ylim(-0.90, 0.70) +
      labs(title = "", x = "Location", y = "Value", fill = "") +
      theme_minimal() +
      scale_fill_manual(values = mycolors) +
      theme(axis.title.y = element_text(size=18),
            axis.title.x = element_text(size=18),
            axis.text.x = element_text(size=16),
            axis.text.y = element_text(size=16),
            legend.title = element_text(size = 18),
            legend.text = element_text(size = 16),legend.position = "none")
    plotlists<- list(rfigs,rfigs2)
    print(cowplot::plot_grid(plotlist = plotlists, ncol = 2, labels = c("A", "B"), label_size = 17))
  add_legend(0.85, 1.15, legend="Truth",
             pch=19, col="black",
             horiz=TRUE, bty='n', cex=1.8)
  add_legend("topleft", legend=substitute(paste(bold(Modeltype))),
             col=c("black", "red"),
             horiz=TRUE, bty='n', cex=1.5)
  #add_legend("topright", legend="Truth",
  #           pch=19, col="black",
  #           horiz=TRUE, bty='n', cex=1.1)
}


RecoverInfUABG.plot<- function(inf.object, true_u, true_a_k, true_B, true_G, Modeltype="", time=60, ndept=9, burn.in=100){

  nstrain<- length(true_a_k)
  if(!is.data.frame(inf.object)){
    fullu.draws<- as.data.frame(inf.object$draws(variables = "uconstrained")[,1,])
    fullak.draws<- as.data.frame(inf.object$draws(variables = "a_k")[,1,])
    fullB.draws<- as.data.frame(inf.object$draws(variables = "B")[,1,])
    fullG.draws<- as.data.frame(inf.object$draws(variables = c("G12","G21"))[,1,])
  }else{
    fullu.draws<- inf.object[-(1:burn.in), 5+time+12+(1:ndept)]
    fullak.draws<- inf.object[-(1:burn.in), 5+time+12+ndept+nstrain+(1:nstrain)]
    fullB.draws<- inf.object[-(1:burn.in), 5+time+12+ndept+(1:nstrain)]
    fullG.draws<- inf.object[-(1:burn.in), 1:2]
  }

  thinning<- numeric(floor(nrow(fullu.draws)/10))
  thinning[1]<- 10
  for(i in 2:length(thinning)){
    thinning[i]<- thinning[i-1] + 10
  }
  u.draws<- fullu.draws[thinning, ]
  ak.draws<- fullak.draws[thinning, ]
  B.draws<- fullB.draws[thinning, ]
  G.draws<- fullG.draws[thinning, ]
  sim.u<- true_u
  sim.ak<- true_a_k
  sim.B<- true_B
  sim.G<- true_G

  # Violin plot for spatial components and intercepts
  spatcomp<- data.frame(value = c(u.draws[, 1], u.draws[, 2], u.draws[, 3], u.draws[, 4],
                                  u.draws[, 5], u.draws[, 6], u.draws[, 7], u.draws[, 8],
                                  u.draws[, 9]), group = factor(rep(c("u1", "u2", "u3", "u4", "u5", "u6", "u7", "u8", "u9"), each = nrow(u.draws))))
  spatcomp$group <- factor(spatcomp$group, levels = unique(spatcomp$group))
  library(ggplot2)
  library(RColorBrewer)
  mycolors <- c(rep(c("blue", "red"), 4), "blue")
  rfigs<- ggplot(spatcomp, aes(x = group, y = value, fill = group)) +
    geom_violin(trim = FALSE, alpha = 0.7) +
    geom_point(x = 1, y = sim.u[1], size = 2, shape = 19) +
    geom_point(x = 2, y = sim.u[2], size = 2, shape = 19) +
    geom_point(x = 3, y = sim.u[3], size = 2, shape = 19) +
    geom_point(x = 4, y = sim.u[4], size = 2, shape = 19) +
    geom_point(x = 5, y = sim.u[5], size = 2, shape = 19) +
    geom_point(x = 6, y = sim.u[6], size = 2, shape = 19) +
    geom_point(x = 7, y = sim.u[7], size = 2, shape = 19) +
    geom_point(x = 8, y = sim.u[8], size = 2, shape = 19) +
    geom_point(x = 9, y = sim.u[9], size = 2, shape = 19) +
    ylim(-0.90, 0.70) +
    labs(title = "", x = "Location", y = "Value", fill = "") +
    theme_minimal() +
    scale_fill_manual(values = mycolors) +
    theme(axis.title.y = element_text(size=18),
          axis.title.x = element_text(size=18),
          axis.text.x = element_text(size=16),
          axis.text.y = element_text(size=16),
          legend.title = element_text(size = 18),
          legend.text = element_text(size = 16),legend.position = "none")

  Acomp<- data.frame(value = c(ak.draws[, 1], ak.draws[, 2], ak.draws[, 3], ak.draws[, 4],
                                  ak.draws[, 5]), group = factor(rep(c("a1", "a2", "a3", "a4", "a5"), each = nrow(ak.draws))))
  Acomp$group <- factor(Acomp$group, levels = unique(Acomp$group))
  mycolors <- c(rep("purple", 5))
  rfigs2<- ggplot(Acomp, aes(x = group, y = value, fill = group)) +
    geom_violin(trim = FALSE, alpha = 0.7) +
    geom_point(x = 1, y = sim.ak[1], size = 2, shape = 19) +
    geom_point(x = 2, y = sim.ak[2], size = 2, shape = 19) +
    geom_point(x = 3, y = sim.ak[3], size = 2, shape = 19) +
    geom_point(x = 4, y = sim.ak[4], size = 2, shape = 19) +
    geom_point(x = 5, y = sim.ak[5], size = 2, shape = 19) +
    #ylim(-0.90, 0.70) +
    labs(title = "", x = "Location", y = "Value", fill = "") +
    theme_minimal() +
    scale_x_discrete(labels = c("a1" = expression(a[1]),
                                "a2" = expression(a[2]),
                                "a3" = expression(a[3]),
                                "a4" = expression(a[4]),
                                "a5" = expression(a[5]))) +
    scale_fill_manual(values = c(rep("purple", 5)),
                      labels = c(expression(a[1]),
                                 expression(a[2]),
                                 expression(a[3]),
                                 expression(a[4]),
                                 expression(a[5])),
                      name = "True value") +
    theme(axis.title.y = element_text(size=18),
          axis.title.x = element_text(size=18),
          axis.text.x = element_text(size=16),
          axis.text.y = element_text(size=16),
          legend.title = element_text(size = 18),
          legend.text = element_text(size = 16),legend.position = "none")

  Bcomp<- data.frame(value = c(G.draws[, 1], G.draws[, 2], B.draws[, 1], B.draws[, 2], B.draws[, 3], B.draws[, 4],
                               B.draws[, 5]), group = factor(rep(paste0("B",1:(nstrain+2)), each = nrow(B.draws))))
  Bcomp$group <- factor(Bcomp$group, levels = unique(Bcomp$group))
  rfigs3<- ggplot(Bcomp, aes(x = group, y = value, fill = group)) +
    geom_violin(trim = FALSE, alpha = 0.7) +
    geom_point(x = 1, y = sim.G[1], size = 2, shape = 19) +
    geom_point(x = 2, y = sim.G[2], size = 2, shape = 19) +
    geom_point(x = 3, y = sim.B[1], size = 2, shape = 19) +
    geom_point(x = 4, y = sim.B[2], size = 2, shape = 19) +
    geom_point(x = 5, y = sim.B[3], size = 2, shape = 19) +
    geom_point(x = 6, y = sim.B[4], size = 2, shape = 19) +
    geom_point(x = 7, y = sim.B[5], size = 2, shape = 19) +
    #ylim(-0.90, 0.70) +
    labs(title = "", x = "Location", y = "Value", fill = "") +
    theme_minimal() +
    scale_fill_brewer(palette = "Set2") +
    scale_x_discrete(labels = c("B1" = expression(gamma[0][1]),
                                "B2" = expression(gamma[1][0]),
                                "B3" = expression(beta[1]),
                                "B4" = expression(beta[2]),
                                "B5" = expression(beta[3]),
                                "B6" = expression(beta[4]),
                                "B7" = expression(beta[5]))) +
    scale_fill_manual(values = c(rep("purple", 7)),
                      labels = c(expression(gamma[0][1]),
                                 expression(gamma[1][0]),
                                 expression(beta[1]),
                                 expression(beta[2]),
                                 expression(beta[3]),
                                 expression(beta[4]),
                                 expression(beta[5])),
                      name = "True value") +
    theme(axis.title.y = element_text(size=18),
          axis.title.x = element_text(size=18),
          axis.text.x = element_text(size=16),
          axis.text.y = element_text(size=16),
          legend.title = element_text(size = 18),
          legend.text = element_text(size = 16),legend.position = "none")
  plotlists<- list(rfigs,rfigs2,rfigs3)
  print(cowplot::plot_grid(plotlist = plotlists, ncol = 3, labels = c("A", "B", "C"), rel_widths = c(1,0.8,1.2), label_size = 17))
  add_legend(0.85, 1.15, legend="Truth",
             pch=19, col="black",
             horiz=TRUE, bty='n', cex=1.8)
  add_legend("topleft", legend=substitute(paste(bold(Modeltype))),
             horiz=TRUE, bty='n', cex=1.5)
  #add_legend("topright", legend="Truth",
  #           pch=19, col="black",
  #           horiz=TRUE, bty='n', cex=1.1)
}

Outbreakfigures<- function(matrix_list, BitsMatrix, labelLetter=""){
  time<- ncol(matrix_list[[1]])
  ndept<- nrow(matrix_list[[1]])
  pdf(paste0("Alloutbreaks",labelLetter,".pdf"), paper="special", width=24,height=24, pointsize=12)
  par(mfrow=c(4,8))
  for(i in 1:length(matrix_list)){
    X_it<- matrix_list[[i]]
    smallxit<- X_it[c(1,3,5,7,9), ]
    bigxit<- X_it[c(2,4,6,8), ]
    bigsmallxit<- X_it[c(2,4,6,8,1,3,5,7,9), ]
    image(x=1:time, y=1:ndept, t(bigsmallxit), main =paste(Bits[i,], collapse = ","), axes=F, ylab="spatial location", xlab="Time [month]", cex.lab=1.80, cex.main=2.5)
    abline(h=4.5, col="black", lty=2)
    #custom Y-axis
    axis(2, at=seq(1, 4, length.out=4), labels=c("u2", "u4", "u6", "u8"), col = "red", col.axis="red", lwd.ticks = 1, las = 1, lwd=0, cex.axis = 1.8, cex.lab=1.8)
    axis(2, at=seq(5, 9, length.out=5), labels=c("u1", "u3", "u5", "u7", "u9"), col = "blue", col.axis="blue", lwd.ticks = 1, las = 1, lwd=0, cex.axis = 1.8, cex.lab=1.8)
    #custom X-axis
    axis(1, cex.axis = 1.8)
  }
  add_legend("topleft", legend=substitute(paste(bold(labelLetter))),
             horiz=TRUE, bty='n', cex=2.0)

  dev.off()
}
