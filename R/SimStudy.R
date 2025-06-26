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
py<- params[["py"]]
n_patches<- nrow(sim_adjmat)
max_time<- params[["max_time"]]
dt<- params[["dt"]]
amplitude<- params[["amplitude"]]

set.seed(1234)
# Initialize state variables
#N <- c(rbind(rpois(5,500), rpois(5,1000)))
N<- rep(c(500,1000), 5)
N <- N[-10]

SIRS_constrained <- function(S, M_matrix, N, params) {
  neq <- length(S)
  eq <- numeric(neq)

  S_constrained <- pmin(pmax(S, 0), N)

  I <- N - S_constrained
  R <- N - S_constrained - I

  for (i in 1:neq) {
    eq[i] <- S_constrained[i] * params[["rho"]] / (params[["rho"]] + params[["gamma"]]) * (M_matrix[i,] %*% I) -
        ((params[["rho"]] * params[["gamma"]]) / (params[["gamma"]] + params[["rho"]]) * I[i])
  }
  return(eq)
}

solveSteadystate <- function(init.condition, params, N, sim_adjmat){
  beta <- params[["beta"]]
  gamma <- params[["gamma"]]
  rho <- params[["rho"]]
  py <- params[["py"]]

  #transmission matrix M_matrix
  M_matrix <- matrix(0, nrow = n_patches, ncol = n_patches)
  for (i in 1:n_patches) {
    neighbors <- which(sim_adjmat[i, ] == 1)
    pop_neighbors <- sum(N[neighbors])
    for (j in 1:n_patches) {
      if (i == j) {
        M_matrix[i, j] <- (beta * (1 - py * length(neighbors))) /
          ((1 - py * length(neighbors)) * N[i] + py * pop_neighbors)
      } else if (j %in% neighbors) {
        M_matrix[i, j] <- (beta * py) /
          ((1 - py * length(neighbors)) * N[i] + py * pop_neighbors)
      }
    }
  }
  res <- BBsolve(par = init.condition, fn = function(S) SIRS_constrained(S, M_matrix, N, params))
  if(res$convergence == 0){
  return(pmin(pmax(res$par, 0), N))
  }else{
    return(N)
  }
}

library(lhs)
B0<- seq(0, 1000, by=10)

set.seed(0)
X <- randomLHS(n = length(B0), k = n_patches)
N_Mat<- matrix(rep(N, nrow(X)), nrow = nrow(X), ncol = n_patches, byrow = T)
Y<- X*N_Mat
colnames(Y)<- paste0("S",1:9)

Sstdstates<- matrix(NA, nrow = nrow(Y), ncol = n_patches)
for(i in 1:nrow(Y)){
  Sstdstates[i, ]<- solveSteadystate(init.condition = Y[i,], params = params, N = N, sim_adjmat = sim_adjmat)
}
Sstdstates

Istdstates<- (rho/(gamma+rho))*(N_Mat-Sstdstates)
Rstdstates<- N_Mat - Sstdstates - Istdstates

Sstdstates/N_Mat
Istdstates/N_Mat
Rstdstates/N_Mat
Sstdstates
Istdstates
Rstdstates
#Veccolors<- randomcoloR::randomColor(nrow(Y))
Veccolors<- rainbow(10)

#Phase plot
par(mfrow=c(3,3))
for(j in 1:n_patches){
  plot(0, type = "n", xlim = c(0, 1), ylim = c(0, 1), xlab = "R", ylab = "I", main = paste("Subpopulation ", j), cex.lab = 1.7, cex.axis = 1.0, cex.main=2.0)
  points(Rstdstates[,j]/N_Mat[,j],Istdstates[,j]/N_Mat[,j], col= Veccolors, pch=3, cex=1.7)
  grid()
}
