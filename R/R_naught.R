library(gdata)

# Define adjacency matrix
sim_adjmat <- matrix(0, nrow = 9, ncol = 9)
uppertriang <- c(1,0,1,0,0,0,0,0,1,0,1,0,0,0,0,0,0,1,0,0,0,1,0,1,0,0,1,0,1,0,0,0,1,1,0,1)
gdata::upperTriangle(sim_adjmat, byrow=TRUE) <- uppertriang
gdata::lowerTriangle(sim_adjmat, byrow=FALSE) <- uppertriang

# Define parameters
beta<- params[["beta"]]
gamma<- params[["gamma"]]
rho<- params[["rho"]]
pi<- params[["pi"]]
n_patches<- nrow(sim_adjmat)
max_time<- params[["max_time"]]
dt<- params[["dt"]]
amplitude<- params[["amplitude"]]

set.seed(1234)
# Initialize state variables
#N <- c(rbind(rpois(5,500), rpois(5,1000)))
N<- rep(c(500,1000), 5)
N <- N[-10]

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

S<- c(round(solve(M_matrix) %*% rep(gamma, n_patches)))
I<- c(round((rho/(gamma+rho))*(N-S)))
R<- N-S-I

T_matrix <- matrix(0, nrow=n_patches, ncol=n_patches)
for (i in 1:n_patches) {
  neighbors <- which(sim_adjmat[i, ] == 1)
  pop_neighbors <- sum(N[neighbors])

  for (j in 1:n_patches) {
    if (i == j) {
      # Within sub-population transmission
      T_matrix[i, j] <- beta * (1 - pi * length(neighbors)) * N[i] /
        ((1 - pi * length(neighbors)) * N[i] + pi * pop_neighbors)
    } else if (j %in% neighbors) {
      # Transmission from neighboring patches
      T_matrix[i, j] <- beta * pi * N[i] /
        ((1 - pi * length(neighbors)) * N[i] + pi * pop_neighbors)
    }
  }
}

SIGMA_matrix <- diag(rep(gamma, n_patches))

NGM_matrix <- T_matrix %*% solve(SIGMA_matrix)

R0_est <- max(eigen(NGM_matrix)$values)

print(R0_est)
