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
