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

  pdf("AllPosteriorpredictiveData7.pdf", paper="special", width=28,height=24, pointsize=12)
  par(mfrow=c(8,9), mar = c(5.1, 6.8, 4, 0.8))
  #c(bottom, left, top, right) ==> specification order

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

    thinning<- numeric(floor(nrow(fullr.draws)/50))
    thinning[1]<- 50
    for(i in 2:length(thinning)){
      thinning[i]<- thinning[i-1] + 50
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

    if(Model==0){
      plot(0, type = "n", xlim = c(1,ncol(y)), ylim = c(min(lCI.Y, y, na.rm = T), max(uCI.Y, y, na.rm = T)), main = paste("Location ", k), ylab = "Case counts", xlab = "Time [month]", cex.axis = 1.7, cex.lab=1.8, cex.main=2.5)
    }else{
      plot(0, type = "n", xlim = c(1,ncol(y)), ylim = c(min(lCI.Y, y, na.rm = T), max(uCI.Y, y, na.rm = T)), main = "", ylab = "Case counts", xlab = "Time [month]", cex.axis = 1.7, cex.lab=1.8, cex.main=2.5)
    }

    #plot(0, type = "n", xlim = c(1,ncol(y)), ylim = c(min(lCI.Y, y, na.rm = T), max(uCI.Y, y, na.rm = T)), main = paste0("Model = ", Model, ",", " Location = ", k), ylab = "Case counts", xlab = "Time [month]", cex.axis = 1.7, cex.lab=1.8, cex.main=2.0)
    polygon(c(1:length(inf.Y), rev(1:length(inf.Y))), c(lCI.Y, rev(uCI.Y)),
            col = "pink", border = NA)
    lines(1:length(inf.Y), inf.Y, col = "red", lty=1)
    points(1:ncol(y), y[k, ], pch = 19)
    grid()
    }
    legendary::labelFig(LETTERS[Model+1], adj = c(-12.25, -0.50), font=2, cex=2.8)
    #c(left/right, top/bottom) ==> specification order
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

RecoverInfUA.plot <- function(inf.object, true_u, true_a_k,
                               Modeltype = "", burn.in = 100) {
  library(ggplot2)
  library(cowplot)

  # parameter dimensions
  n_u  <- length(true_u)
  n_ak <- length(true_a_k)

  # Extract posterior draws
  if (!is.data.frame(inf.object)) {
    fullu.draws <- as.data.frame(inf.object$draws(variables = "uconstrained")[, 1, ])
    fullak.draws <- as.data.frame(inf.object$draws(variables = "a_k")[, 1, ])
  } else {
    fullu.draws <- inf.object[-(1:burn.in), startsWith(colnames(inf.object), "u")]
    fullak.draws <- inf.object[-(1:burn.in), startsWith(colnames(inf.object), "a")]
  }

  # thinning every 10
  thinning <- seq(10, nrow(fullu.draws), by = 10)
  u.draws  <- fullu.draws[thinning, , drop = FALSE]
  ak.draws <- fullak.draws[thinning, , drop = FALSE]

  # ---- U violin plot ----
  u_labels <- do.call(expression, lapply(1:n_u, function(i) {
    bquote(u[.(i)])
  }))

  spatcomp_u <- data.frame(
    value = as.vector(as.matrix(u.draws)),
    group = factor(rep(paste0("u", 1:n_u), each = nrow(u.draws)))
  )

  rfigs_u <- ggplot(spatcomp_u, aes(x = group, y = value, fill = group)) +
    geom_violin(trim = FALSE, alpha = 0.7) +
    geom_point(data = data.frame(x = factor(paste0("u", 1:n_u), levels = levels(spatcomp_u$group)),
                                 y = true_u),
               aes(x = x, y = y),
               size = 2, shape = 19, inherit.aes = FALSE) +
    ylim(-0.40, 0.40) +
    labs(x = "Location", y = "Value", fill = "") +
    theme_minimal() +
    scale_fill_manual(values = rep(c("blue", "red"), length.out = n_u)) +
    scale_x_discrete(labels = u_labels) +   # <-- add dynamic labels here
    theme(axis.title = element_text(size = 17),
          axis.text = element_text(size = 16),
          legend.position = "none")

  # ---- a_k violin plot ----
  ak_labels <- do.call(expression, lapply(1:n_ak, function(i) {
    bquote(a[.(i)])
  }))

  comp_ak <- data.frame(
    value = as.vector(as.matrix(ak.draws)),
    group = factor(rep(paste0("a", 1:n_ak), each = nrow(ak.draws)))
  )

  rfigs_ak <- ggplot(comp_ak, aes(x = group, y = value, fill = group)) +
    geom_violin(trim = FALSE, alpha = 0.7) +
    geom_point(data = data.frame(x = factor(paste0("a", 1:n_ak), levels = levels(comp_ak$group)),
                                 y = true_a_k),
               aes(x = x, y = y),
               size = 2, shape = 19, inherit.aes = FALSE) +
    ylim(-14.5, -12) +
    labs(x = "Intercepts", y = "Value", fill = "") +
    theme_minimal() +
    scale_fill_manual(values = rep("purple", n_ak)) +
    scale_x_discrete(labels = ak_labels) +   # <-- add dynamic labels here
    theme(axis.title = element_text(size = 17),
          axis.text = element_text(size = 16),
          legend.position = "none")

  # ---- combine plots ----
  plotlists <- list(rfigs_u, rfigs_ak)
  print(cowplot::plot_grid(plotlist = plotlists, ncol = 2,
                           labels = c("A", "B"),
                           rel_widths = c(1.25, 1), label_size = 17))

  # Legends
  add_legend(0.85, 1.15, legend = "Truth",
             pch = 19, col = "black",
             horiz = TRUE, bty = 'n', cex = 1.8)
  add_legend("topleft", legend = substitute(paste(bold(Modeltype))),
             horiz = TRUE, bty = 'n', cex = 1.5)
}


RecoverInfUAB.plot <- function(inf.object, true_u, true_a_k, true_B,
                                 Modeltype = "", burn.in = 100) {
  library(ggplot2)
  library(cowplot)

  # parameter dimensions
  n_u  <- length(true_u)
  n_ak <- length(true_a_k)
  n_B  <- length(true_B)

  # Extract posterior draws
  if (!is.data.frame(inf.object)) {
    fullu.draws <- as.data.frame(inf.object$draws(variables = "uconstrained")[, 1, ])
    fullak.draws <- as.data.frame(inf.object$draws(variables = "a_k")[, 1, ])
    fullB.draws <- as.data.frame(inf.object$draws(variables = "B")[, 1, ])
  } else {
    fullu.draws <- inf.object[-(1:burn.in), startsWith(colnames(inf.object), "u")]
    fullak.draws <- inf.object[-(1:burn.in), startsWith(colnames(inf.object), "a")]
    fullB.draws <- inf.object[-(1:burn.in), startsWith(colnames(inf.object), "B")]
  }

  # thinning every 10
  thinning <- seq(10, nrow(fullu.draws), by = 10)
  u.draws  <- fullu.draws[thinning, , drop = FALSE]
  ak.draws <- fullak.draws[thinning, , drop = FALSE]
  B.draws  <- fullB.draws[thinning, , drop = FALSE]

  # ---- U violin plot ----
  u_labels <- do.call(expression, lapply(1:n_u, function(i) {
    bquote(u[.(i)])
  }))

  spatcomp_u <- data.frame(
    value = as.vector(as.matrix(u.draws)),
    group = factor(rep(paste0("u", 1:n_u), each = nrow(u.draws)))
  )

  rfigs_u <- ggplot(spatcomp_u, aes(x = group, y = value, fill = group)) +
    geom_violin(trim = FALSE, alpha = 0.7) +
    geom_point(data = data.frame(x = factor(paste0("u", 1:n_u), levels = levels(spatcomp_u$group)),
                                 y = true_u),
               aes(x = x, y = y),
               size = 2, shape = 19, inherit.aes = FALSE) +
    ylim(-0.40, 0.40) +
    labs(x = "Location", y = "Value", fill = "") +
    theme_minimal() +
    scale_fill_manual(values = rep(c("blue", "red"), length.out = n_u)) +
    scale_x_discrete(labels = u_labels) +   # <-- add dynamic labels here
    theme(axis.title = element_text(size = 17),
          axis.text = element_text(size = 16),
          legend.position = "none")

  # ---- a_k violin plot ----
  ak_labels <- do.call(expression, lapply(1:n_ak, function(i) {
    bquote(a[.(i)])
  }))

  comp_ak <- data.frame(
    value = as.vector(as.matrix(ak.draws)),
    group = factor(rep(paste0("a", 1:n_ak), each = nrow(ak.draws)))
  )

  rfigs_ak <- ggplot(comp_ak, aes(x = group, y = value, fill = group)) +
    geom_violin(trim = FALSE, alpha = 0.7) +
    geom_point(data = data.frame(x = factor(paste0("a", 1:n_ak), levels = levels(comp_ak$group)),
                                 y = true_a_k),
               aes(x = x, y = y),
               size = 2, shape = 19, inherit.aes = FALSE) +
    ylim(-14.5, -12) +
    labs(x = "Intercepts", y = "Value", fill = "") +
    theme_minimal() +
    scale_fill_manual(values = rep("purple", n_ak)) +
    scale_x_discrete(labels = ak_labels) +   # <-- add dynamic labels here
    theme(axis.title = element_text(size = 17),
          axis.text = element_text(size = 16),
          legend.position = "none")

  # ---- B violin plot ----
  B_labels <- do.call(expression, lapply(1:n_B, function(i) {
    bquote(beta[.(i)])
  }))

  comp_B <- data.frame(
    value = as.vector(as.matrix(B.draws)),
    group = factor(rep(paste0("B", 1:n_B), each = nrow(B.draws)))
  )

  rfigs_B <- ggplot(comp_B, aes(x = group, y = value, fill = group)) +
    geom_violin(trim = FALSE, alpha = 0.7) +
    geom_point(data = data.frame(x = factor(paste0("B", 1:n_B), levels = levels(comp_B$group)),
                                 y = true_B),
               aes(x = x, y = y),
               size = 2, shape = 19, inherit.aes = FALSE) +
    ylim(0.5, 2.5) +
    labs(x = "Regression coeff.", y = "Value", fill = "") +
    theme_minimal() +
    scale_fill_manual(values = rep("purple", n_B)) +
    scale_x_discrete(labels = B_labels) +   # <-- add dynamic labels here
    theme(axis.title = element_text(size = 17),
          axis.text = element_text(size = 16),
          legend.position = "none")

  # ---- combine plots ----
  plotlists <- list(rfigs_u, rfigs_ak, rfigs_B)
  print(cowplot::plot_grid(plotlist = plotlists, ncol = 3,
                           labels = c("A", "B", "C"),
                           rel_widths = c(1.25, 1, 1), label_size = 17))

  # Legends
  add_legend(0.85, 1.15, legend = "Truth",
             pch = 19, col = "black",
             horiz = TRUE, bty = 'n', cex = 1.8)
  add_legend("topleft", legend = substitute(paste(bold(Modeltype))),
             horiz = TRUE, bty = 'n', cex = 1.5)
}


RecoverInfUABG.plot <- function(inf.object, true_u, true_a_k, true_B, true_G,
                                Modeltype = "", burn.in = 100) {
  library(ggplot2)
  library(cowplot)

  # parameter dimensions
  n_u  <- length(true_u)
  n_ak <- length(true_a_k)
  n_B  <- length(true_B)
  n_G  <- length(true_G)

  # Extract posterior draws
  if (!is.data.frame(inf.object)) {
    fullu.draws <- as.data.frame(inf.object$draws(variables = "uconstrained")[, 1, ])
    fullak.draws <- as.data.frame(inf.object$draws(variables = "a_k")[, 1, ])
    fullB.draws <- as.data.frame(inf.object$draws(variables = "B")[, 1, ])
    fullG.draws <- as.data.frame(inf.object$draws(variables = paste0("G", 1:n_G))[, 1, ])
  } else {
    fullu.draws <- inf.object[-(1:burn.in), startsWith(colnames(inf.object), "u")]
    fullak.draws <- inf.object[-(1:burn.in), startsWith(colnames(inf.object), "a")]
    fullB.draws <- inf.object[-(1:burn.in), startsWith(colnames(inf.object), "B")]
    fullG.draws <- inf.object[-(1:burn.in), startsWith(colnames(inf.object), "G")]
  }

  # thinning every 10
  thinning <- seq(10, nrow(fullu.draws), by = 10)
  u.draws  <- fullu.draws[thinning, , drop = FALSE]
  ak.draws <- fullak.draws[thinning, , drop = FALSE]
  B.draws  <- fullB.draws[thinning, , drop = FALSE]
  G.draws  <- fullG.draws[thinning, , drop = FALSE]

  # ---- U violin plot ----
  u_labels <- do.call(expression, lapply(1:n_u, function(i) {
    bquote(u[.(i)])
  }))

  spatcomp_u <- data.frame(
    value = as.vector(as.matrix(u.draws)),
    group = factor(rep(paste0("u", 1:n_u), each = nrow(u.draws)))
  )

  rfigs_u <- ggplot(spatcomp_u, aes(x = group, y = value, fill = group)) +
    geom_violin(trim = FALSE, alpha = 0.7) +
    geom_point(data = data.frame(x = factor(paste0("u", 1:n_u), levels = levels(spatcomp_u$group)),
                                 y = true_u),
               aes(x = x, y = y),
               size = 2, shape = 19, inherit.aes = FALSE) +
    ylim(-0.40, 0.40) +
    labs(x = "Location", y = "Value", fill = "") +
    theme_minimal() +
    scale_fill_manual(values = rep(c("blue", "red"), length.out = n_u)) +
    scale_x_discrete(labels = u_labels) +   # <-- add dynamic labels here
    theme(axis.title = element_text(size = 17),
          axis.text = element_text(size = 16),
          legend.position = "none")

  # ---- a_k violin plot ----
  ak_labels <- do.call(expression, lapply(1:n_ak, function(i) {
    bquote(a[.(i)])
  }))

  comp_ak <- data.frame(
    value = as.vector(as.matrix(ak.draws)),
    group = factor(rep(paste0("a", 1:n_ak), each = nrow(ak.draws)))
  )

  rfigs_ak <- ggplot(comp_ak, aes(x = group, y = value, fill = group)) +
    geom_violin(trim = FALSE, alpha = 0.7) +
    geom_point(data = data.frame(x = factor(paste0("a", 1:n_ak), levels = levels(comp_ak$group)),
                                 y = true_a_k),
               aes(x = x, y = y),
               size = 2, shape = 19, inherit.aes = FALSE) +
    ylim(-14.5, -12) +
    labs(x = "Intercepts", y = "Value", fill = "") +
    theme_minimal() +
    scale_fill_manual(values = rep("purple", n_ak)) +
    scale_x_discrete(labels = ak_labels) +   # <-- add dynamic labels here
    theme(axis.title = element_text(size = 17),
          axis.text = element_text(size = 16),
          legend.position = "none")

  # ---- B violin plot ----
  B_labels <- do.call(expression, lapply(1:n_B, function(i) {
    bquote(beta[.(i)])
  }))

  comp_B <- data.frame(
    value = as.vector(as.matrix(B.draws)),
    group = factor(rep(paste0("B", 1:n_B), each = nrow(B.draws)))
  )

  rfigs_B <- ggplot(comp_B, aes(x = group, y = value, fill = group)) +
    geom_violin(trim = FALSE, alpha = 0.7) +
    geom_point(data = data.frame(x = factor(paste0("B", 1:n_B), levels = levels(comp_B$group)),
                                 y = true_B),
               aes(x = x, y = y),
               size = 2, shape = 19, inherit.aes = FALSE) +
    ylim(0.5, 2.5) +
    labs(x = "Regression coeff.", y = "Value", fill = "") +
    theme_minimal() +
    scale_fill_manual(values = rep("purple", n_B)) +
    scale_x_discrete(labels = B_labels) +   # <-- add dynamic labels here
    theme(axis.title = element_text(size = 17),
          axis.text = element_text(size = 16),
          legend.position = "none")

  # ---- G violin plot ----
  G_index_pairs <- rep(list(c(0,1), c(1,0)), ncol(G.draws)/2)

  G_labels <- do.call(expression, lapply(G_index_pairs, function(idx) {
    bquote(gamma[.(idx[1])][.(idx[2])])
  }))

  comp_G <- data.frame(
    value = as.vector(as.matrix(G.draws)),
    group = factor(rep(paste0("G", 1:n_G), each = nrow(G.draws)))
  )

  rfigs_G <- ggplot(comp_G, aes(x = group, y = value, fill = group)) +
    geom_violin(trim = FALSE, alpha = 0.7) +
    geom_point(data = data.frame(x = factor(paste0("G", 1:n_G),
                                            levels = levels(comp_G$group)),
                                 y = true_G),
               aes(x = x, y = y),
               size = 2, shape = 19, inherit.aes = FALSE) +
    ylim(0, 1) +
    labs(x = "Transition prob.", y = "Value", fill = "") +
    theme_minimal() +
    scale_fill_manual(values = rep("purple", n_G)) +
    scale_x_discrete(labels = G_labels) +   # <-- add dynamic labels here
    theme(axis.title = element_text(size = 17),
          axis.text = element_text(size = 16),
          legend.position = "none")

  # ---- combine plots ----
  plotlists <- list(rfigs_u, rfigs_ak, rfigs_B, rfigs_G)
  print(cowplot::plot_grid(plotlist = plotlists, ncol = 4,
                           labels = c("A", "B", "C", "D"),
                           rel_widths = c(1.25, 1, 1, 1.5), label_size = 17))

  # Legends
  add_legend(0.85, 1.15, legend = "Truth",
             pch = 19, col = "black",
             horiz = TRUE, bty = 'n', cex = 1.8)
  add_legend("topleft", legend = substitute(paste(bold(Modeltype))),
             horiz = TRUE, bty = 'n', cex = 1.5)
}


#Copula parameters plot
RecoverInfcopula.plot <- function(all.infobjects, true_copParams, burn.in = 100){
  library(ggplot2)
  library(cowplot)

  n_c <- length(true_copParams)

  n <- (sqrt(8*n_c + 1) - 1)/2

  c_labels <- list()
  k <- 1
  for (i in 1:n) {
    for (j in i:n) {
      c_labels[[k]] <- bquote(rho[.(i) * .(j+1)])
      k <- k + 1
    }
  }

  AllFigs<- list()

  for(r in 1:2){
    inf.object<- all.infobjects[[r]]
  if (!is.data.frame(inf.object)) {
    fullc.draws <- as.data.frame(inf.object$draws(variables = "copParams")[, 1, ])
  } else {
    fullc.draws <- inf.object[-(1:burn.in), startsWith(colnames(inf.object), "co")]
  }

  # thinning every 10
  thinning <- seq(10, nrow(fullc.draws), by = 10)
  c.draws  <- fullc.draws[thinning, , drop = FALSE]

  copulaP <- data.frame(
    value = as.vector(as.matrix(c.draws)),
    group = factor(rep(paste0("c", 1:n_c), each = nrow(c.draws)))
  )

  rfig_c <- ggplot(copulaP, aes(x = group, y = value, fill = group)) +
    geom_violin(trim = FALSE, alpha = 0.7) +
    geom_point(data = data.frame(x = factor(paste0("c", 1:n_c), levels = levels(copulaP$group)),
                                 y = true_copParams),
               aes(x = x, y = y),
               size = 2, shape = 19, inherit.aes = FALSE) +
    ylim(-1.1, 1.1) +
    labs(x = "Copula parameter", y = "Value", fill = "") +
    theme_minimal() +
    scale_fill_manual(values = rep("purple", n_c)) +
    scale_x_discrete(labels = c_labels) +
    theme(axis.title = element_text(size = 17),
          axis.text = element_text(size = 16),
          legend.position = "none")

  AllFigs[[r]]<- rfig_c
}

print(cowplot::plot_grid(plotlist = AllFigs, ncol = 2,
                         labels = c("A", "B"),
                          label_size = 17))

  # Legends
  add_legend(0.85, 1.15, legend = "Truth",
             pch = 19, col = "black",
             horiz = TRUE, bty = 'n', cex = 1.8)
}


#RecoverInfG2.plot <- function(inf.object, true_G, nstrain,
#                                 Modeltype = "", burn.in = 100) {
#  library(ggplot2)
#  library(cowplot)

  # parameter dimensions
#  nstate<- 2^nstrain
#  true_G<- as.numeric(t(true_G))

  # Extract posterior draws
#  if (!is.data.frame(inf.object)) {
#    fullG.draws <- as.data.frame(inf.object$draws(variables = paste0("G", 1:nstate))[, 1, ])
#  } else {
#    fullG.draws <- inf.object[-(1:burn.in), startsWith(colnames(inf.object), "G")]
#  }

  # thinning every 10
#  thinning <- seq(10, nrow(fullG.draws), by = 10)

#  G_index_pairs <- lapply(0:(nstate-1), function(i) {
#    lapply(0:(nstate-1), function(j) c(i,j))
#  })
#  G_index_pairs <- unlist(G_index_pairs, recursive = FALSE)

#  G_labels <- do.call(expression, lapply(G_index_pairs, function(idx) {
#    bquote(gamma[.(idx[1])][.(idx[2])])
#  }))

#  allFigs<- list()

#  for(i in 1:nstate){
#    index<- nstate * (i-1) + 1
#    G.draws  <- fullG.draws[thinning, index:(i*nstate), drop = FALSE]

#  comp_G <- data.frame(
#    value = as.vector(as.matrix(G.draws)),
#    group = factor(rep(paste0("G", 1:nstate), each = nrow(G.draws)))
#  )

#  rfigs_G <- ggplot(comp_G, aes(x = group, y = value, fill = group)) +
#    geom_violin(trim = FALSE, alpha = 0.7) +
#    geom_point(data = data.frame(x = factor(paste0("G", 1:nstate),
#                                            levels = levels(comp_G$group)),
#                                 y = true_G[index:(i*nstate)]),
#               aes(x = x, y = y),
#               size = 2, shape = 19, inherit.aes = FALSE) +
#    labs(x = NULL, y = NULL, fill = NULL) +
#    theme_minimal() +
#    theme(
#      axis.text = element_blank(),
#      axis.ticks = element_blank(),
#      legend.position = "none"
#    ) +
#    scale_fill_manual(values = rep("purple", ncol(G.draws)))
#    labs(x = "Transition prob.", y = "Value", fill = "") +
#    theme_minimal() +
#    scale_fill_manual(values = rep("purple", nstate)) +
#    scale_x_discrete(labels = G_labels[index:(i*nstate)]) +   # <-- add dynamic labels here
#    theme(axis.title = element_text(size = 17),
#          axis.text = element_text(size = 16),
#          legend.position = "none")

#  allFigs[[i]]<- rfigs_G
#}

  # ---- combine plots ----
#  print(cowplot::plot_grid(plotlist = allFigs, ncol = 1))

  # Legends
#  add_legend(0.85, 1.15, legend = "Truth",
#             pch = 19, col = "black",
#             horiz = TRUE, bty = 'n', cex = 1.8)
#  add_legend("topleft", legend = substitute(paste(bold(Modeltype))),
#             horiz = TRUE, bty = 'n', cex = 1.5)
#}

RecoverInfG.plot <- function(inf.object, true_G, burn.in = 100) {
  # parameter dimensions
  true_G<- as.numeric(t(true_G))

  # Extract posterior draws
  if (!is.data.frame(inf.object)) {
    fullG.draws <- as.data.frame(inf.object$draws(variables = paste0("G", 1:nstate))[, 1, ])
  } else {
    fullG.draws <- inf.object[-(1:burn.in), startsWith(colnames(inf.object), "G")]
  }

  # thinning every 10
  thinning <- seq(10, nrow(fullG.draws), by = 10)
  G.draws  <- fullG.draws[thinning, ]

  par(mfrow=c(3, 3))
  for (i in 1:ncol(G.draws)) {
    hist(G.draws[, i], main = colnames(fullG.draws)[i], xlab ="", col = "white", border = "black")
    abline(v=true_G[i], col="red")
  }
  # Legends
  add_legend(0.85, 1.15, legend = "Truth",
             lty = 1, col = "red",
             horiz = TRUE, bty = 'n', cex = 1.8)
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


# Super-impose true values
perstrainOutbreakfigures<- function(Truth_array, matrix_array, Outbreaktype=""){
  time<- ncol(matrix_array[,,1])
  ndept<- nrow(matrix_array[,,1])
  nstrain<- dim(matrix_array)[3]
  pdf(paste0("Alloutbreaksperstrain",Outbreaktype,".pdf"), paper="special", width=15,height=9, pointsize=12)
  par(mfrow=c(2,3))
  for(i in 1:nstrain){
    Truth<- Truth_array[,,i]
    smallTruth<- Truth[c(1,3,5,7,9), ]
    bigTruth<- Truth[c(2,4,6,8), ]
    bigsmallTruth<- Truth[c(2,4,6,8,1,3,5,7,9), ]

    X_it<- matrix_array[,,i]
    smallxit<- X_it[c(1,3,5,7,9), ]
    bigxit<- X_it[c(2,4,6,8), ]
    bigsmallxit<- X_it[c(2,4,6,8,1,3,5,7,9), ]
    image(x=1:time, y=1:ndept, t(bigsmallxit), zlim = c(0,1),  main =paste("Strain", i), axes=F, ylab="spatial location", xlab="Time [month]", cex.lab=1.80, cex.main=2.5)
    abline(h=4.5, col="black", lty=2)
    #custom Y-axis
    axis(2, at=seq(1, 4, length.out=4), labels=c("u2", "u4", "u6", "u8"), col = "red", col.axis="red", lwd.ticks = 1, las = 1, lwd=0, cex.axis = 1.8, cex.lab=1.8)
    axis(2, at=seq(5, 9, length.out=5), labels=c("u1", "u3", "u5", "u7", "u9"), col = "blue", col.axis="blue", lwd.ticks = 1, las = 1, lwd=0, cex.axis = 1.8, cex.lab=1.8)
    #custom X-axis
    axis(1, cex.axis = 1.8)
    # Build grid of row/col indices
    grid <- expand.grid(x = seq(nrow(bigsmallTruth)),
                        y = seq(ncol(bigsmallTruth)))
    out <- transform(grid, z = bigsmallTruth[as.matrix(grid)])

    # Select true outbreak cells
    outbreakcell <- out$z == 1

    # Overlay truth as colored dots at cell centers
    points(out$y[outbreakcell],   # time axis
           out$x[outbreakcell],   # location axis
           pch = 16, cex = 1.5, col = "magenta")
  }
  add_legend("topright", legend=substitute(paste(bold("Truth"))),
             pch=16, col="magenta",
             horiz=TRUE, bty='n', cex=1.8)

  add_legend(0.50, -0.4, legend=substitute(paste(bold(Outbreaktype))),
             col="black",
             horiz=TRUE, bty='n', cex=3.0)

  dev.off()
}



RecoverInfAllRS.plot<- function(allinf.object, alltrue_r, alltrue_s){

  allobject<- length(allinf.object)
  par(mfrow=c(2,4))
  Modeltype<- c("Model 0", "Model 1", "Model 2", "Model 5")

for(a in 1:allobject){
  inf.object<- allinf.object[[a]]
  true_r<- alltrue_r[[a]]
  true_s<- alltrue_s[[a]]

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
  }

  plot(0, type = "n", xlim = c(1,length(inf.r)), ylim = c(-0.5, 0.5), ylab = "Trend component", xlab = "Time [Month]", main = Modeltype[[a]])
  polygon(c(1:length(inf.r), rev(1:length(inf.r))), c(lCI.r, rev(uCI.r)),
          col = "pink", border = NA)
  lines(1:length(inf.r), inf.r, col="red")
  points(1:length(inf.r), true_r, pch = 19)
  grid()

  plot(0, type = "n", xlim = c(1,length(inf.s)), ylim = c(-1.6, 1.6), ylab = "Seasonal component", xlab = "Season [Month]", main = Modeltype[[a]])
  polygon(c(1:length(inf.s), rev(1:length(inf.s))), c(lCI.s, rev(uCI.s)),
          col = "pink", border = NA)
  lines(1:length(inf.s), inf.s, col="red")
  points(1:length(inf.s), true_s, pch = 19)
  grid()
  }
  add_legend(0.46, 1.15, legend=c("Truth", "Posterior means"), lty=c(NA, 1),
             pch=c(19, NA), col=c("black", "red"),
             horiz=TRUE, bty='n', cex=1.5)
  }


publicationfig4<- function(all.infobjects, all.simobjects, adjmat){
  time<- ncol(all.simobjects[[2]][[1]])
  ndept<- nrow(all.simobjects[[2]][[1]])
  pdf("posterheatmap.pdf", paper="special", width=21,height=12, pointsize=14)
  #plot.new()
  par(mfrow=c(2,4), mar = c(5.1, 6.8, 6.8, 0.8))
  #par(mfrow=c(8,8), mar = c(5.1, 6.8, 4, 0.8))
  #c(bottom, left, top, right) ==> specification order

  X_it<- all.simobjects[[2]][[6]]
  smallxit<- X_it[c(1,3,5,7,9), ]
  bigxit<- X_it[c(2,4,6,8), ]
  bigsmallxit<- X_it[c(2,4,6,8,1,3,5,7,9), ]
  #image(x=1:time, y=1:ndept, t(bigsmallxit), main ="", axes=F, ylab="spatial location", xlab="Time [month]", cex.lab=1.80)

  #abline(h=4.5, col="black", lty=2)
  #custom Y-axis
  #axis(2, at=seq(1, 4, length.out=4), labels=c("u2", "u4", "u6", "u8"), col = "red", col.axis="red", lwd.ticks = 1, las = 1, lwd=0, cex.axis = 1.8, cex.lab=1.8)
  #axis(2, at=seq(5, 9, length.out=5), labels=c("u1", "u3", "u5", "u7", "u9"), col = "blue", col.axis="blue", lwd.ticks = 1, las = 1, lwd=0, cex.axis = 1.8, cex.lab=1.8)
  #custom X-axis
  #axis(1, cex.axis = 1.8)
  #legendary::labelFig("A", adj = c(-0.15, 0.10), font=2, cex=1.8)

  for(i in 1:7){
    sim.object<- all.simobjects[[i+1]]
    inf.object<- all.infobjects[[i+1]]
    Model<- i

    y<- sim.object[[1]]
    e_it<- sim.object[[2]]
    sim.r<- sim.object[[3]]
    sim.s<- sim.object[[4]]
    sim.u<- sim.object[[5]]
    X_it<- sim.object[[6]]

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
      image(x=1:time, y=1:ndept, t(bigsmallmeanxit), main = "", axes=F, ylab = "Spatial location", xlab = "Time [month]", cex.lab=1.80)
      abline(h=4.5, col="black", lty=2)
      #custom Y-axis
      axis(2, at=seq(1, 4, length.out=4), labels=c("u2", "u4", "u6", "u8"), col = "red", col.axis="red", lwd.ticks = 1, las = 1, lwd=0, cex.axis = 1.8, cex.lab=1.8)
      axis(2, at=seq(5, 9, length.out=5), labels=c("u1", "u3", "u5", "u7", "u9"), col = "blue", col.axis="blue", lwd.ticks = 1, las = 1, lwd=0, cex.axis = 1.8, cex.lab=1.8)
      #custom X-axis
      axis(1, cex.axis = 1.8)
      # Build grid of row/col indices
      grid <- expand.grid(x = seq(nrow(bigsmallxit)),
                          y = seq(ncol(bigsmallxit)))
      out <- transform(grid, z = bigsmallxit[as.matrix(grid)])

      # Select true outbreak cells
      outbreakcell <- out$z == 1

      # Overlay truth as colored dots at cell centers
      points(out$y[outbreakcell],   # time axis
             out$x[outbreakcell],   # location axis
             pch = 16, cex = 2.0, col = "magenta")
      legendary::labelFig(LETTERS[Model], adj = c(-0.15, 0.10), font=2, cex=2.0)
    }
  }
  par(mar = c(5, 16, 7, 15))
  #c(bottom, left, top, right) ==> specification order

  zseq <- seq(0, 1, length.out = 101)
  xseq <- c(0, 1)
  zmat <- matrix(zseq[-1], nrow = 1)

  image(x = xseq, y = zseq, z = zmat,
        axes = FALSE, xlab = "", ylab = "", main = "")

  axis(4, at = seq(0, 1, 0.2), labels = seq(0, 1, 0.2), las = 1, cex.axis=2.0)
  mtext("Posterior probability of outbreak", side = 4, line = 5.0, cex = 1.7)
  box()
  add_legend("topright", legend=c("Truth", "Small cities", "Large cities"), lty=c(NA, 1, 1),
             pch=c(16, NA, NA), col=c("magenta", "blue", "red"),
             horiz=TRUE, bty='n', cex=2.0)
  dev.off()
}


publicationfig7<- function(all.infobjects, realdata, adjmat){
  y<- realdata[[1]]
  e_it<- realdata[[2]]
  countries<- realdata[[3]]
  time<- ncol(y)
  ndept<- nrow(y)
  #listofOutProbs<- list()

  pdf("Allconnectedheatmap.pdf", paper="special", width=21,height=12, pointsize=14)
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
      mean_Xit<- mean_Xit[ndept:1, ]
      par(mar = c(4, 7.5, 4, 1))
      image(x=1:time, y=1:ndept, t(mean_Xit), main = "", axes=F, ylab = "", xlab = "Time [month/year]", cex.lab=1.8)
      #custom Y-axis
      axis(2, at=seq(1, length(countries), length.out=length(countries)), labels=rev(countries), lwd.ticks = 1, las = 1, lwd=0, cex.axis = 1.40)
      #custom X-axis
      years<- 2013:2019
      axis(1, at = seq(1, 84, by = 12), labels = years, cex.axis = 1.8)
      legendary::labelFig(LETTERS[Model], adj = c(-0.15, 0.05), font=2, cex=2.0)
    }
    #listofOutProbs[[Model]]<- sum_Xit/length(thinning)
  }

  par(mar = c(4, 15, 4, 14))
  zseq <- seq(0, 1, length.out = 101)
  xseq <- c(0, 1)
  zmat <- matrix(zseq[-1], nrow = 1)

  image(x = xseq, y = zseq, z = zmat,
        axes = FALSE, xlab = "", ylab = "", main = "")

  axis(4, at = seq(0, 1, 0.2), labels = seq(0, 1, 0.2), las = 1, cex.axis=2.0)
  mtext("Posterior probability of outbreak", side = 4, line = 5.0, cex = 2.0)
  box()
  dev.off()
  #return(listofOutProbs)
}




#par(fig = c(0.9, 1, 0.1, 0.9), new = TRUE, mar = c(5, 0, 4, 4))
# Define coordinates properly
#zseq <- seq(0, 1, length.out = 101)
#xseq <- c(0, 1)

#zmat <- matrix(zseq[-1], nrow = 1)
#image(x = xseq, y = zseq, z = zmat
#      , axes = FALSE, xlab = "", ylab = "", yaxt="n", xaxt="n", main="")

# Add axis and label
#axis(4, at = seq(0, 1, 0.3), labels = seq(0, 1, 0.3), las = 1)
#mtext("Posterior probability of outbreak", side = 4, line = 2.5)
#box()


#excessFits<- list()
#for(i in 1:8){
#  for(j in 1:8){
#    index<- 8*(i-1)+j
#    excessFits[[index]] <- readRDS(paste0("C:/Users/Matthew Adeoye/Downloads/simfit",i-1,j-1,".RDS"))
#  }
#}

publicationfigS3<- function(all.infobjects, all.simobjects, adjmat){
  time<- ncol(all.simobjects[[2]][[1]])
  ndept<- nrow(all.simobjects[[2]][[1]])

  pdf("supplmA.pdf", paper="special", width=24,height=24, pointsize=12)
  par(mfrow=c(8,8), mar = c(5.1, 6.8, 4, 0.8))
  #c(bottom, left, top, right) ==> specification order

  for(j in 1:8){
    sim.object<- all.simobjects[[j]]
    y<- sim.object[[1]]
    e_it<- sim.object[[2]]
    ind<- j*8-8

    for(k in 1:8){
      inf.object<- all.infobjects[[ind+k]]
      Model<- k-1

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

      thinning<- numeric(floor(nrow(fullr.draws)/100))
      thinning[1]<- 100
      for(i in 2:length(thinning)){
        thinning[i]<- thinning[i-1] + 100
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
      inf.Y<- colMeans(sum_Y)
      uCI.Y<- posterior_interval_custom(as.matrix.data.frame(sum_Y))[,2]
      lCI.Y<- posterior_interval_custom(as.matrix.data.frame(sum_Y))[,1]

      if(j==1){
      if(Model == 0){
        romanfig<- 0
      }else{
        romanfig<- as.roman(Model)
      }
      plot(0, type = "n", xlim = c(1,ncol(y)), ylim = c(0, 300), ylab = "overall case counts", xlab = "Time [month]", cex.lab=1.80, cex.axis = 1.7, main=paste("Model", romanfig), cex.main=2.5)
      }else{
        plot(0, type = "n", xlim = c(1,ncol(y)), ylim = c(0, 300), ylab = "overall case counts", xlab = "Time [month]", cex.lab=1.80, cex.axis = 1.7, main="", cex.main=2.5)
      }
      polygon(c(1:length(inf.Y), rev(1:length(inf.Y))), c(lCI.Y, rev(uCI.Y)),
              col = "pink", border = NA)
      lines(1:length(inf.Y), inf.Y, col = "red", lty=1)
      points(1:ncol(y), colSums(y), pch = 19)
      grid()
    }
    legendary::labelFig(LETTERS[j], adj = c(-10.97, -0.50), font=2, cex=2.8)
  }
#  add_legend("topright", legend=c("Truth", "Posterior means"), lty=c(NA, 1),
#             pch=c(19, NA), col=c("black", "red"),
#             horiz=TRUE, bty='n', cex=1.8)
  dev.off()
}


PosteriorFOI<- function(all.infobjects,  all.simobjects, adjmat){
  time<- ncol(all.simobjects[[1]][[1]])
  ndept<- nrow(all.simobjects[[1]][[1]])
  month_indexes <- ((1:time - 1) %% 12) + 1

  pdf("PosteriorFOI.pdf", paper="special", width=28,height=24, pointsize=12)
  par(mfrow=c(8,9), mar = c(5.1, 6.8, 4, 0.8))
  #c(bottom, left, top, right) ==> specification order

  for(m in 1:8){
    inf.object<- all.infobjects[[m]]
    sim.object<- all.simobjects[[m]]
    Model<- m-1

    design_matrix_func <- get(paste0("DesignMatrixModel", Model))
    z_it <- design_matrix_func(y, adjmat)[[1]]
    z_it2 <- design_matrix_func(y, adjmat)[[2]]

    y<- sim.object[[1]]
    e_it<- sim.object[[2]]
    true_r<- sim.object[[3]]
    true_s<- sim.object[[4]]
    true_u<- sim.object[[5]]
    if(Model==0){
      r_mat <- matrix(rep(true_r, each = ndept), nrow = ndept)
      s_mat <- matrix(rep(true_s[month_indexes], each = ndept), nrow = ndept)
      u_mat <- matrix(rep(true_u, times = time), nrow = ndept)
      truelambda <- exp(r_mat + s_mat + u_mat)
    }else if(Model %in% c(1,2,7)){
      r_mat <- matrix(rep(true_r, each = ndept), nrow = ndept)
      s_mat <- matrix(rep(true_s[month_indexes], each = ndept), nrow = ndept)
      u_mat <- matrix(rep(true_u, times = time), nrow = ndept)
      trueXit<- sim.object[[6]]
      truelambda <- exp(r_mat + s_mat + u_mat + trueXit * (1.65 * z_it))
    }else if(Model==3){
      r_mat <- matrix(rep(true_r, each = ndept), nrow = ndept)
      s_mat <- matrix(rep(true_s[month_indexes], each = ndept), nrow = ndept)
      u_mat <- matrix(rep(true_u, times = time), nrow = ndept)
      trueXit<- sim.object[[6]]
      truelambda <- exp(r_mat + s_mat + u_mat + trueXit * ((1.25 * z_it) + (0.75 * z_it2)))
    }else if(Model==4){
      r_mat <- matrix(rep(true_r, each = ndept), nrow = ndept)
      s_mat <- matrix(rep(true_s[month_indexes], each = ndept), nrow = ndept)
      u_mat <- matrix(rep(true_u, times = time), nrow = ndept)
      trueXit<- sim.object[[6]]
      truelambda <- exp(r_mat + s_mat + u_mat + trueXit * (0.55 * z_it))
    }else if(Model==5){
      r_mat <- matrix(rep(true_r, each = ndept), nrow = ndept)
      s_mat <- matrix(rep(true_s[month_indexes], each = ndept), nrow = ndept)
      u_mat <- matrix(rep(true_u, times = time), nrow = ndept)
      trueXit<- sim.object[[6]]
      truelambda <- exp(r_mat + s_mat + u_mat + trueXit * (0.49 * z_it))
    }else if(Model==6){
      r_mat <- matrix(rep(r, each = ndept), nrow = ndept)
      s_mat <- matrix(rep(s[month_indexes], each = ndept), nrow = ndept)
      u_mat <- matrix(rep(u, times = time), nrow = ndept)
      trueXit<- sim.object[[6]]
      truelambda <- exp(r_mat + s_mat + u_mat + trueXit * ((0.35 * z_it) + (0.20 * z_it2)))
    }

    for(k in 1:ndept){

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
      pred_lambda<- matrix(NA, nrow = ndept, ncol = time)

      thinning<- numeric(floor(nrow(fullr.draws)/50))
      thinning[1]<- 50
      for(i in 2:length(thinning)){
        thinning[i]<- thinning[i-1] + 50
      }

      G12.draws<- fullG12.draws[thinning]
      G21.draws<- fullG21.draws[thinning]
      r.draws<- fullr.draws[thinning, ]
      s.draws<- fulls.draws[thinning, ]
      u.draws<- fullu.draws[thinning, ]
      d.draws<- fulld.draws[thinning]

      sum_lambda<- matrix(NA, nrow = length(thinning), ncol = time)
      sum_Xit<- matrix(0, nrow = ndept, ncol = time)

      if(Model == 0){
        for(index in 1:length(thinning)){
          for(i in 1:ndept){
            for(t in 1:time){
              m<- (t - 1) %% 12 + 1
              pred_lambda[i, t] <- exp(r.draws[index, t] + s.draws[index, m] + u.draws[index, i])
              #pred_Y[i, t]<- rpois(1, Exlambda_it)
            }
          }
          sum_lambda[index, ]<- pred_lambda[k, ]
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
              pred_lambda[i, t] <- exp(r.draws[index, t] + s.draws[index, m] + u.draws[index, i] + P_Xit * z_it[i, t] * B.draws[index])
              #pred_Y[i, t]<- rpois(1, Exlambda_it)
            }
          }
          sum_lambda[index, ]<- pred_lambda[k,]
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
              pred_lambda[i, t] <- exp(r.draws[index, t] + s.draws[index, m] + u.draws[index, i] + P_Xit * z_it[i, t] * B.draws[index, 1] + P_Xit * z_it2[i, t] * B.draws[index, 2])
              #pred_Y[i, t]<- rpois(1, Exlambda_it)
            }
          }
          sum_lambda[index, ]<- pred_lambda[k, ]
        }
      }

      inf.lambda<- colMeans(sum_lambda)
      uCI.lambda<- posterior_interval_custom(as.matrix.data.frame(sum_lambda))[,2]
      lCI.lambda<- posterior_interval_custom(as.matrix.data.frame(sum_lambda))[,1]

      if(Model==0){
        plot(0, type = "n", xlim = c(1,ncol(truelambda)), yaxt="n", ylim = c(0, 8e-05), main = paste("Location ", k), ylab = expression("Force of infection ("%*% 10^-5 * ")"), xlab = "Time [month]", cex.axis = 1.7, cex.lab=1.8, cex.main=2.5)
        axis(2,
             at = seq(0, 8e-5, by = 1e-5),
             labels = seq(0, 8, by = 1), cex.axis = 1.7)
        #        axis(2,
#             at = seq(0, 6e-6, by = 1e-6),
#             labels = seq(0, 0.6, by = 0.1))
        }else{
        plot(0, type = "n", xlim = c(1,ncol(truelambda)), yaxt="n", ylim = c(0, 8e-05), main = "", ylab = expression("Force of infection ("%*% 10^-5 * ")"), xlab = "Time [month]", cex.axis = 1.7, cex.lab=1.8, cex.main=2.5)
          axis(2,
               at = seq(0, 8e-5, by = 1e-5),
               labels = seq(0, 8, by = 1), cex.axis = 1.7)
      }

      polygon(c(1:length(inf.lambda), rev(1:length(inf.lambda))), c(lCI.lambda, rev(uCI.lambda)),
              col = "pink", border = NA)
      lines(1:length(inf.lambda), inf.lambda, col = "red", lty=1)
      points(1:ncol(truelambda), truelambda[k, ], pch = 19)
      grid()
    }
    legendary::labelFig(LETTERS[Model+1], adj = c(-12.25, -0.50), font=2, cex=2.8)
    #c(left/right, top/bottom) ==> specification order
  }
  #  add_legend("topright", legend=c("Truth", "Posterior means"), lty=c(NA, 1),
  #             pch=c(19, NA), col=c("black", "red"),
  #             horiz=TRUE, bty='n', cex=1.8)
  dev.off()
}


# Section 4 - Application to invasive meningococcal disease
# Monthly reported cases (Figure 6 in paper)
#data(ApplicationCounts, package = "DetectOutbreaks")
#data(ApplicationPopulation, package = "DetectOutbreaks")
#data(ApplicationAdjMat, package = "DetectOutbreaks")
#y<- ApplicationCounts
#e_it<- ApplicationPopulation
#Adjmat<- ApplicationAdjMat
#DetectOutbreaks::sim.plot(y, names=colnames(Adjmat))

#start_time <- Sys.time()
#ApplicationfitsHMC<- list()
#for(i in 1:8){
#  Model<- i-1
#  HMCFit<- DetectOutbreaks::infer(y = y, e_it = e_it, Model = Model, adjmat = Adjmat, Stan = T, iter = 5000, GPU = F, verbose = T)
#  ApplicationfitsHMC[[i]]<- HMCFit
#}
#end_time <- Sys.time()
#time_taken<- end_time - start_time
#print(time_taken)

#start_time <- Sys.time()
#ApplicationfitsMCMC<- list()
#for(i in 1:8){
#  Model<- i-1
#  MCMCFit<- DetectOutbreaks::infer(y = y, e_it = e_it, Model = Model, adjmat = Adjmat, Stan = F, iter = 5000, GPU = F, verbose = T, num_iteration = 5000)
#  ApplicationfitsMCMC[[i]]<- MCMCFit
#}
#end_time <- Sys.time()
#time_taken<- end_time - start_time
#print(time_taken)
#cmdstanr::check_cmdstan_toolchain(fix = TRUE)


#cumulative_standardizedcasesper100k<- rowSums((y/e_it) * 100000, na.rm = T)
publicationfigs8<- function(cumulative_standardizedcasesper100k){

    #get required regions
    poly <- cshapes::cshp(date=as.Date("2019-12-31"), useGW=TRUE)
    #newpoly<- poly[poly$country_name %in% c("Austria","Belgium","Cyprus","Czech Republic",
    #                                      "Denmark","Estonia","Finland","France","German Federal Republic",
    #                                      "Greece","Hungary","Iceland","Ireland","Italy/Sardinia","Latvia",
    #                                      "Lithuania","Luxembourg","Malta","Netherlands","Norway","Poland",
    #                                      "Portugal","Rumania","Slovakia","Slovenia","Spain","Sweden","United Kingdom"), ]


    newpoly<- poly[poly$country_name %in% c("Albania", "Austria", "Belarus (Byelorussia)", "Belgium", "Bosnia-Herzegovina",
                                            "Bulgaria", "Croatia", "Cyprus", "Czech Republic", "Denmark", "Estonia", "Finland",
                                            "France", "German Federal Republic", "Greece", "Hungary", "Iceland",
                                            "Ireland", "Italy/Sardinia", "Kosovo", "Latvia", "Lithuania", "Luxembourg", "Malta",
                                            "Macedonia (FYROM/North Macedonia)", "Moldova", "Montenegro", "Netherlands", "Norway",
                                            "Poland", "Portugal", "Rumania", "Russia (Soviet Union)", "Serbia", "Slovakia", "Slovenia", "Spain", "Sweden",
                                            "Switzerland", "Ukraine", "United Kingdom"), ]

    Europe<- c("Austria","Belgium","Cyprus","Czech Republic",
               "Denmark","Estonia","Finland","France","German Federal Republic",
               "Greece","Hungary","Iceland","Ireland","Italy/Sardinia","Latvia",
               "Lithuania","Luxembourg","Malta","Netherlands","Norway","Poland",
               "Portugal","Rumania","Slovakia","Slovenia","Spain","Sweden","United Kingdom",
               "Albania", "Belarus (Byelorussia)", "Bosnia-Herzegovina", "Bulgaria",
               "Croatia", "Kosovo", "Russia (Soviet Union)",
               "Moldova", "Montenegro", "Macedonia (FYROM/North Macedonia)",
               "Serbia", "Switzerland", "Ukraine")

    sortedpoly <- newpoly[order(newpoly$country_name), ]
    #sortedpoly$country_name<- c("Austria","Belgium","Cyprus","Czechia",
    #                          "Denmark","Estonia","Finland","France","Germany",
    #                          "Greece","Hungary","Iceland","Ireland","Italy","Latvia",
    #                          "Lithuania","Luxembourg","Malta","Netherlands","Norway","Poland",
    #                          "Portugal","Romania","Slovakia","Slovenia","Spain","Sweden","United Kingdom")

    u.draws<- cumulative_standardizedcasesper100k
    #expUi<- c(exp(u.draws), rep(NA, 13))
    expUi<- c(u.draws, rep(NA, 13))
    uidf<- data.frame(expUi = expUi, countryname = Europe)
    sorteduidf<- uidf[order(uidf$countryname), ]
    rrDF<- data.frame(rr=sorteduidf$expUi, gwcode=sortedpoly$gwcode, index=c(1:nrow(uidf)), name=sorteduidf$countryname)

    library(sf)
    sortedpoly<- st_make_valid(newpoly)
    bbox <- st_bbox(c(xmin = -29, ymin = 30, xmax = 40.17875, ymax = 71.15471), crs = st_crs(sortedpoly))
    sortedpoly<- st_crop(sortedpoly, bbox)
    shapefile <- ggplot2::fortify(sortedpoly, region = 'gwcode')

    shp_rrDF<- sp::merge(shapefile, rrDF,
                         by.x="gwcode",
                         by.y="gwcode",
                         all.x=F)

    library(ggplot2)
    library(viridis)

      rfig<- (ggplot(data = shp_rrDF) +
                geom_sf(aes(fill = rr)) +
                #geom_sf_text(aes(label = country_name), size = 3) +
                #scale_fill_viridis(option = "turbo", direction = 1, alpha=1, begin=0.6, end=1, na.value = "lightgrey") +  # Reverse the color scale
                scale_fill_gradient2(low = "lightblue", mid = "blue", high = "red", midpoint = 4, na.value = "lightgrey") +
                coord_sf() + theme_void() +
                ggtitle("") + theme(plot.title = element_text(hjust = 0.55)) +
                labs(fill = "Cumulative cases per 100k") +
                theme(legend.title = element_text(size = 18),
                      legend.text = element_text(size = 16)))
   print(rfig)
}


#postprobsfulldata<- DetectOutbreaks::OutbreakProbability(y=y,e_it=e_it,inf.object=CorrectedGPU_7,adjmat=Adjmat,Model = 7)
#postprobs73Jan2019<- DetectOutbreaks::OutbreakProbability(y=y[,1:73],e_it=e_it[,1:73],inf.object=nowcastfitGPU_11,adjmat=Adjmat,Model = 7)
#postprobs74Feb2019<- DetectOutbreaks::OutbreakProbability(y=y[,1:74],e_it=e_it[,1:74],inf.object=nowcastfitGPU_10,adjmat=Adjmat,Model = 7)
#postprobs75Mar2019<- DetectOutbreaks::OutbreakProbability(y=y[,1:75],e_it=e_it[,1:75],inf.object=nowcastfitGPU_9,adjmat=Adjmat,Model = 7)
#postprobs76Apr2019<- DetectOutbreaks::OutbreakProbability(y=y[,1:76],e_it=e_it[,1:76],inf.object=nowcastfitGPU_8,adjmat=Adjmat,Model = 7)
#postprobs77May2019<- DetectOutbreaks::OutbreakProbability(y=y[,1:77],e_it=e_it[,1:77],inf.object=nowcastfitGPU_7,adjmat=Adjmat,Model = 7)
#postprobs78Jun2019<- DetectOutbreaks::OutbreakProbability(y=y[,1:78],e_it=e_it[,1:78],inf.object=nowcastfitGPU_6,adjmat=Adjmat,Model = 7)
#postprobs79Jul2019<- DetectOutbreaks::OutbreakProbability(y=y[,1:79],e_it=e_it[,1:79],inf.object=nowcastfitGPU_5,adjmat=Adjmat,Model = 7)
#postprobs80Aug2019<- DetectOutbreaks::OutbreakProbability(y=y[,1:80],e_it=e_it[,1:80],inf.object=nowcastfitGPU_4,adjmat=Adjmat,Model = 7)
#postprobs81Sep2019<- DetectOutbreaks::OutbreakProbability(y=y[,1:81],e_it=e_it[,1:81],inf.object=nowcastfitGPU_3,adjmat=Adjmat,Model = 7)
#postprobs82Oct2019<- DetectOutbreaks::OutbreakProbability(y=y[,1:82],e_it=e_it[,1:82],inf.object=nowcastfitGPU_2,adjmat=Adjmat,Model = 7)
#postprobs83Nov2019<- DetectOutbreaks::OutbreakProbability(y=y[,1:83],e_it=e_it[,1:83],inf.object=nowcastfitGPU_1,adjmat=Adjmat,Model = 7)

#par(mfrow=c(3,4))
#plot(postprobsfulldata[,73],postprobs73Jan2019[,73], col="purple", pch=15, main="January 2019", xlab = "Posterior prob. - full data", ylab = "Posterior prob. - shortened data", cex.lab=1.3, cex.main=1.5, cex.axis=1.2, ylim=c(0,1), cex=1.5);grid()
#abline(a = 0, b = 1, col = "grey", lwd = 2.5)
#plot(postprobsfulldata[,74],postprobs74Feb2019[,74], col="purple", pch=15, main="February 2019", xlab = "Posterior prob. - full data", ylab = "Posterior prob. - shortened data", cex.lab=1.3, cex.main=1.5, cex.axis=1.2, ylim=c(0,1), cex=1.5);grid()
#abline(a = 0, b = 1, col = "grey", lwd = 2.5)
#plot(postprobsfulldata[,75],postprobs75Mar2019[,75], col="purple", pch=15, main="March 2019", xlab = "Posterior prob. - full data", ylab = "Posterior prob. - shortened data", cex.lab=1.3, cex.main=1.5, cex.axis=1.2, ylim=c(0,1), cex=1.5);grid()
#abline(a = 0, b = 1, col = "grey", lwd = 2.5)
#plot(postprobsfulldata[,76],postprobs76Apr2019[,76], col="purple", pch=15, main="April 2019", xlab = "Posterior prob. - full data", ylab = "Posterior prob. - shortened data", cex.lab=1.3, cex.main=1.5, cex.axis=1.2, ylim=c(0,1), cex=1.5);grid()
#abline(a = 0, b = 1, col = "grey", lwd = 2.5)
#plot(postprobsfulldata[,77],postprobs77May2019[,77], col="purple", pch=15, main="May 2019", xlab = "Posterior prob. - full data", ylab = "Posterior prob. - shortened data", cex.lab=1.3, cex.main=1.5, cex.axis=1.2, ylim=c(0,1), cex=1.5);grid()
#abline(a = 0, b = 1, col = "grey", lwd = 2.5)
#plot(postprobsfulldata[,78],postprobs78Jun2019[,78], col="purple", pch=15, main="June 2019", xlab = "Posterior prob. - full data", ylab = "Posterior prob. - shortened data", cex.lab=1.3, cex.main=1.5, cex.axis=1.2, ylim=c(0,1), cex=1.5);grid()
#abline(a = 0, b = 1, col = "grey", lwd = 2.5)
#plot(postprobsfulldata[,79],postprobs79Jul2019[,79], col="purple", pch=15, main="July 2019", xlab = "Posterior prob. - full data", ylab = "Posterior prob. - shortened data", cex.lab=1.3, cex.main=1.5, cex.axis=1.2, ylim=c(0,1), cex=1.5);grid()
#abline(a = 0, b = 1, col = "grey", lwd = 2.5)
#plot(postprobsfulldata[,80],postprobs80Aug2019[,80], col="purple", pch=15, main="August 2019", xlab = "Posterior prob. - full data", ylab = "Posterior prob. - shortened data", cex.lab=1.3, cex.main=1.5, cex.axis=1.2, ylim=c(0,1), cex=1.5);grid()
#abline(a = 0, b = 1, col = "grey", lwd = 2.5)
#plot(postprobsfulldata[,81],postprobs81Sep2019[,81], col="purple", pch=15, main="September 2019", xlab = "Posterior prob. - full data", ylab = "Posterior prob. - shortened data", cex.lab=1.3, cex.main=1.5, cex.axis=1.2, ylim=c(0,1), cex=1.5);grid()
#abline(a = 0, b = 1, col = "grey", lwd = 2.5)
#plot(postprobsfulldata[,82],postprobs82Oct2019[,82], col="purple", pch=15, main="October 2019", xlab = "Posterior prob. - full data", ylab = "Posterior prob. - shortened data", cex.lab=1.3, cex.main=1.5, cex.axis=1.2, ylim=c(0,1), cex=1.5);grid()
#abline(a = 0, b = 1, col = "grey", lwd = 2.5)
#plot(postprobsfulldata[,83],postprobs83Nov2019[,83], col="purple", pch=15, main="November 2019", xlab = "Posterior prob. - full data", ylab = "Posterior prob. - shortened data", cex.lab=1.3, cex.main=1.5, cex.axis=1.2, ylim=c(0,1), cex=1.5);grid()
#abline(a = 0, b = 1, col = "grey", lwd = 2.5)


SlidesSimheatmap<- function(all.infobjects, all.simobjects, adjmat){
  time<- ncol(all.simobjects[[2]][[1]])
  ndept<- nrow(all.simobjects[[2]][[1]])
  pdf("SlidesSimheatmap.pdf", paper="special", width=21,height=12, pointsize=14)
  #plot.new()
  par(mfrow=c(2,4), mar = c(5.1, 6.8, 6.8, 0.8))
  #par(mfrow=c(8,8), mar = c(5.1, 6.8, 4, 0.8))
  #c(bottom, left, top, right) ==> specification order

  X_it<- all.simobjects[[2]][[6]]
  smallxit<- X_it[c(1,3,5,7,9), ]
  bigxit<- X_it[c(2,4,6,8), ]
  bigsmallxit<- X_it[c(2,4,6,8,1,3,5,7,9), ]
  #image(x=1:time, y=1:ndept, t(bigsmallxit), main ="", axes=F, ylab="spatial location", xlab="Time [month]", cex.lab=1.80)

  #abline(h=4.5, col="black", lty=2)
  #custom Y-axis
  #axis(2, at=seq(1, 4, length.out=4), labels=c("u2", "u4", "u6", "u8"), col = "red", col.axis="red", lwd.ticks = 1, las = 1, lwd=0, cex.axis = 1.8, cex.lab=1.8)
  #axis(2, at=seq(5, 9, length.out=5), labels=c("u1", "u3", "u5", "u7", "u9"), col = "blue", col.axis="blue", lwd.ticks = 1, las = 1, lwd=0, cex.axis = 1.8, cex.lab=1.8)
  #custom X-axis
  #axis(1, cex.axis = 1.8)
  #legendary::labelFig("A", adj = c(-0.15, 0.10), font=2, cex=1.8)

  for(i in 1:7){
    sim.object<- all.simobjects[[i+1]]
    inf.object<- all.infobjects[[i+1]]
    Model<- i

    y<- sim.object[[1]]
    e_it<- sim.object[[2]]
    sim.r<- sim.object[[3]]
    sim.s<- sim.object[[4]]
    sim.u<- sim.object[[5]]
    X_it<- sim.object[[6]]

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
      image(x=1:time, y=1:ndept, t(bigsmallmeanxit), main = paste("Model", as.roman(Model)), axes=F, ylab = "Spatial location", xlab = "Time [month]", cex.lab=1.80, cex.main=2.0)
      abline(h=4.5, col="black", lty=2)
      #custom Y-axis
      axis(2, at=seq(1, 4, length.out=4), labels=c("u2", "u4", "u6", "u8"), col = "red", col.axis="red", lwd.ticks = 1, las = 1, lwd=0, cex.axis = 1.8, cex.lab=1.8)
      axis(2, at=seq(5, 9, length.out=5), labels=c("u1", "u3", "u5", "u7", "u9"), col = "blue", col.axis="blue", lwd.ticks = 1, las = 1, lwd=0, cex.axis = 1.8, cex.lab=1.8)
      #custom X-axis
      axis(1, cex.axis = 1.8)
      # Build grid of row/col indices
      grid <- expand.grid(x = seq(nrow(bigsmallxit)),
                          y = seq(ncol(bigsmallxit)))
      out <- transform(grid, z = bigsmallxit[as.matrix(grid)])

      # Select true outbreak cells
      outbreakcell <- out$z == 1

      # Overlay truth as colored dots at cell centers
      points(out$y[outbreakcell],   # time axis
             out$x[outbreakcell],   # location axis
             pch = 16, cex = 2.0, col = "magenta")
#      legendary::labelFig(LETTERS[Model], adj = c(-0.15, 0.10), font=2, cex=2.0)
    }
  }
  par(mar = c(5, 16, 7, 15))
  #c(bottom, left, top, right) ==> specification order

  zseq <- seq(0, 1, length.out = 101)
  xseq <- c(0, 1)
  zmat <- matrix(zseq[-1], nrow = 1)

  image(x = xseq, y = zseq, z = zmat,
        axes = FALSE, xlab = "", ylab = "", main = "")

  axis(4, at = seq(0, 1, 0.2), labels = seq(0, 1, 0.2), las = 1, cex.axis=2.0)
  mtext("Posterior probability of outbreak", side = 4, line = 5.0, cex = 1.7)
  box()
  add_legend("topright", legend=c("Truth", "Small cities", "Large cities"), lty=c(NA, 1, 1),
             pch=c(16, NA, NA), col=c("magenta", "blue", "red"),
             horiz=TRUE, bty='n', cex=2.0)
  dev.off()
}



SlidesAppheatmap<- function(all.infobjects, realdata, adjmat){
  y<- realdata[[1]]
  e_it<- realdata[[2]]
  countries<- realdata[[3]]
  time<- ncol(y)
  ndept<- nrow(y)
  #listofOutProbs<- list()

  pdf("SlidesAppheatmap.pdf", paper="special", width=21,height=12, pointsize=14)
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
      mean_Xit<- mean_Xit[ndept:1, ]
      par(mar = c(4, 7.5, 4, 1))
      image(x=1:time, y=1:ndept, t(mean_Xit), main = paste("Model", as.roman(Model)), axes=F, ylab = "", xlab = "Time [month/year]", cex.lab=1.8, cex.main=2.5)
      #custom Y-axis
      axis(2, at=seq(1, length(countries), length.out=length(countries)), labels=rev(countries), lwd.ticks = 1, las = 1, lwd=0, cex.axis = 1.40)
      #custom X-axis
      years<- 2013:2019
      axis(1, at = seq(1, 84, by = 12), labels = years, cex.axis = 1.8)
#      legendary::labelFig(LETTERS[Model], adj = c(-0.15, 0.05), font=2, cex=2.0)
    }
    #listofOutProbs[[Model]]<- sum_Xit/length(thinning)
  }

  par(mar = c(4, 15, 4, 14))
  zseq <- seq(0, 1, length.out = 101)
  xseq <- c(0, 1)
  zmat <- matrix(zseq[-1], nrow = 1)

  image(x = xseq, y = zseq, z = zmat,
        axes = FALSE, xlab = "", ylab = "", main = "")

  axis(4, at = seq(0, 1, 0.2), labels = seq(0, 1, 0.2), las = 1, cex.axis=2.0)
  mtext("Posterior probability of outbreak", side = 4, line = 5.0, cex = 2.0)
  box()
  dev.off()
  #return(listofOutProbs)
}
