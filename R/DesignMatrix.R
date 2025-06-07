DesignMatrixModel0<- function(y, adjacencymatrix){
  ndept<- nrow(y)
  time<- ncol(y)
  z_it<- matrix(0, ndept, time)
  z_it2<- matrix(0, ndept, time)
  return(list(z_it, z_it2))
}

DesignMatrixModel1<- function(y, adjacencymatrix){
  y<- ifelse(is.na(y), 99, y)
  ndept<- nrow(y)
  time<- ncol(y)
  z_it<- matrix(0, ndept, time)
  for(t in 2:time){
    for(i in 1:ndept){
      if(y[i, t-1] > 0)
        z_it[i, t]<- 1
    }
  }
  z_it<- ifelse(is.na(z_it),1,z_it)
  z_it2<- matrix(0, ndept, time)
  return(list(z_it, z_it2))
}

DesignMatrixModel2<- function(y, adjacencymatrix){
  y<- ifelse(is.na(y), 99, y)
  ndept<- nrow(y)
  time<- ncol(y)
  z_it<- matrix(0, ndept, time)

  for(t in 2:time){
    for(i in 1:ndept){
      z_it[i, 1] = 0
      indexes = which(adjacencymatrix[i, ] > 0 & 1:ndept != i)

      #for(b in 1:length(indexes)){
      # neighbor_index <- indexes[b]
      if(length(indexes) != 0){
        if(y[i, t-1] > 0 || any(y[indexes, t-1] > 0)){
          z_it[i, t] = 1
        }else{
          z_it[i, t] = 0
        }
        #}
      }else{
        if(y[i, t-1] > 0){
          z_it[i, t] = 1
        }else{
          z_it[i, t] = 0
        }
      }
    }
  }
  z_it<- ifelse(is.na(z_it),1,z_it)
  z_it2<- matrix(0, ndept, time)
  return(list(z_it, z_it2))
}

DesignMatrixModel3<- function(y, adjacencymatrix){
  y<- ifelse(is.na(y), 99, y)
  ndept<- nrow(y)
  time<- ncol(y)
  z_it<- matrix(0, ndept, time)
  z_it2<- matrix(0, ndept, time)

  for(t in 2:time){
    for(i in 1:ndept){
      z_it[i, 1] = 0
      if(y[i, t-1] > 0){
        z_it[i, t]<- 1
      }else{
        z_it[i, t]<- 0
      }
    }
  }

  for(t in 2:time){
    for(i in 1:ndept){
      indexes = which(adjacencymatrix[i, ] > 0 & 1:ndept != i)
      z_it2[i, 1] = 0
      #flag = 0

      if(length(indexes)>0){
        #for(b in 1:length(indexes)){
        # neighbor_index <- indexes[b]
        if(any(y[indexes, t-1] > 0)){
          z_it2[i, t] = 1
          # flag = 1
          #break;
        }else{
          z_it2[i, t]= 0
        }
        #}
      }else{
        z_it2[i, t] = 0
      }
    }
  }
  z_it<- ifelse(is.na(z_it),1,z_it)
  z_it2<- ifelse(is.na(z_it2),1,z_it2)
  return(list(z_it, z_it2))
}

DesignMatrixModel4<- function(y, adjacencymatrix){
  y<- ifelse(is.na(y), -1, y)
  ndept<- nrow(y)
  time<- ncol(y)
  z_it<- matrix(0, ndept, time)
  for(t in 2:time){
    for(i in 1:ndept){
      z_it[i, t]<- log(y[i, t-1] + 1)
    }
  }
  z_it<- ifelse(is.nan(z_it), 0, z_it)
  z_it<- ifelse(is.infinite(z_it), 0, z_it)
  z_it2<- matrix(0, ndept, time)
  return(list(z_it, z_it2))
}

DesignMatrixModel5<- function(y, adjacencymatrix){
  y<- ifelse(is.na(y), -1, y)
  ndept<- nrow(y)
  time<- ncol(y)
  z_it<- matrix(0, ndept, time)

  for(t in 2:time){
    for(i in 1:ndept){
      indexes <- which(adjacencymatrix[i, ] > 0 & 1:ndept != i)

      sum_neighbors <- 0
      for(b in indexes){
        sum_neighbors<- sum_neighbors + y[b, t-1]
      }
      z_it[i, t]<- sum_neighbors
    }
  }
  for(t in 2:time){
    for(i in 1:ndept){
      z_it[i, t] = log(y[i, t-1] + z_it[i, t] + 1)
    }
  }
  z_it<- ifelse(is.nan(z_it), 0, z_it)
  z_it<- ifelse(is.infinite(z_it), 0, z_it)
  z_it2<- matrix(0, ndept, time)
  return(list(z_it, z_it2))
}

DesignMatrixModel6 <- function(y, adjacencymatrix) {
  y<- ifelse(is.na(y), -1, y)
  ndept <- nrow(y)
  time <- ncol(y)
  z_it <- matrix(0, ndept, time)
  z_it2 <- matrix(0, ndept, time)

  for (t in 2:time) {
    for (i in 1:ndept) {
      z_it[i, t] <- log(y[i, t-1] + 1)
    }
  }

  for(t in 2:time){
    for(i in 1:ndept){
      indexes <- which(adjacencymatrix[i, ] > 0 & 1:ndept != i)

      sum_neighbors <- 0
      for(b in indexes){
        sum_neighbors<- sum_neighbors + y[b, t-1]
      }
      z_it2[i, t]<- log(sum_neighbors + 1)
    }
  }
  z_it<- ifelse(is.nan(z_it), 0, z_it)
  z_it<- ifelse(is.infinite(z_it), 0, z_it)
  z_it2<- ifelse(is.nan(z_it2), 0, z_it2)
  z_it2<- ifelse(is.infinite(z_it2), 0, z_it2)
  return(list(z_it, z_it2))
}

DesignMatrixModel7<- function(y, adjacencymatrix){
  ndept<- nrow(y)
  time<- ncol(y)
  z_it<- matrix(1, ndept, time)
  z_it2<- matrix(0, ndept, time)
  return(list(z_it, z_it2))
}
