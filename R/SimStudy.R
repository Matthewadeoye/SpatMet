library(BB)
# Define adjacency matrix
sim_adjmat <- matrix(0, nrow = 9, ncol = 9)
uppertriang <- c(1,0,1,0,0,0,0,0,1,0,1,0,0,0,0,0,0,1,0,0,0,1,0,1,0,0,1,0,1,0,0,0,1,1,0,1)
gdata::upperTriangle(sim_adjmat, byrow=TRUE) <- uppertriang
gdata::lowerTriangle(sim_adjmat, byrow=FALSE) <- uppertriang

# Define parameters
beta<- params[["beta"]]
gamma<- params[["gamma"]]
rho<- params[["rho"]]
epsilon<- params[["epsilon"]]
n_patches<- nrow(sim_adjmat)
max_time<- params[["max_time"]]
dt<- params[["dt"]]
amplitude<- params[["amplitude"]]

set.seed(1234)
# Initialize state variables
N <- c(rbind(rpois(5,500), rpois(5,1000)))
N <- N[-10]

SIRS<- function(s, M_matrix, N, params) {
  neq<- nrow(M_matrix)
  eq <- rep(NA, neq)
  I<- rep(NA, neq)
  for(i in 1:neq){
    I[i]<- N[i]-s[i]
  }
  for(i in 1:neq){
    eq[i] <- abs(s[i] * params[["rho"]]/(params[["rho"]]+params[["gamma"]]) * (M_matrix[i,] %*% I) - ((params[["rho"]]*params[["gamma"]])/(params[["gamma"]]+params[["rho"]]) * I[i]))
  }
  return(eq)
}

simulate.SIRS<- function(init.condition, params){
  beta<- params[["beta"]]
  gamma<- params[["gamma"]]
  rho<- params[["rho"]]

  #Obtain steady state conditions
  M_matrix <- matrix(0, nrow=n_patches, ncol=n_patches)
  for (i in 1:n_patches) {
    neighbors <- which(sim_adjmat[i, ] == 1)
    pop_neighbors <- sum(N[neighbors])

    for (j in 1:n_patches) {
      if (i == j) {
        # Within sub-population transmission
        M_matrix[i, j] <- (beta * (1 - pi * length(neighbors)) /
                             ((1 - pi * length(neighbors)) * N[i] + pi * pop_neighbors))
      } else if (j %in% neighbors) {
        # Transmission from neighboring patches
        M_matrix[i, j] <- (beta * pi  /
                             ((1 - pi * length(neighbors)) * N[i] + pi * pop_neighbors))
      }
    }
  }

  chk<- BBsolve(par = init.condition, fn = function(s) SIRS(s, M_matrix, N, params))
  return(chk$par)
}


library(lhs)
B0<- seq(0, 500, by=50)

set.seed(0)
X <- randomLHS(n = length(B0), k = n_patches)
Y <- X
lens <- rep(length(B0), n_patches)
for (i in 1:n_patches){
  Y_index <- 1 + floor(Y[,i] * length(B0))
  Y[,i] <- B0[Y_index]
}
colnames(Y)<- paste0("S",1:9)
Y

Sstdstates<- matrix(NA, nrow = nrow(Y), ncol = n_patches)
for(i in 1:nrow(Y)){
  Sstdstates[i, ]<- simulate.SIRS(Y[i,], params = params)
}
Sstdstates

N_Mat<- matrix(rep(N, nrow(Sstdstates)), nrow = nrow(Sstdstates), ncol = n_patches, byrow = T)
Istdstates<- (rho/(gamma+rho))*(N_Mat-Sstdstates)
Rstdstates<- N_Mat - Sstdstates

Sstdstates/N_Mat
Istdstates/N_Mat
Rstdstates/N_Mat
Veccolors<- randomcoloR::randomColor(nrow(Y))

#Phase plot
par(mfrow=c(3,3))
for(j in 1:n_patches){
  plot(0, type = "n", xlim = c(0, max(Rstdstates)), ylim = c(0, max(Istdstates)), xlab = "R", ylab = "I", main = paste("Subpopulation ", j), cex.lab = 1.7, cex.axis = 1.7, cex.main=2.0)
  points(Rstdstates[,j],Istdstates[,j], col= Veccolors, pch=3, cex=1)
}
