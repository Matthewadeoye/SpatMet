library(gdata)
library(randomcoloR)

# Define adjacency matrix
sim_adjmat <- matrix(0, nrow = 9, ncol = 9)
uppertriang <- c(1,0,1,0,0,0,0,0,1,0,1,0,0,0,0,0,0,1,0,0,0,1,0,1,0,0,1,0,1,0,0,0,1,1,0,1)
gdata::upperTriangle(sim_adjmat, byrow=TRUE) <- uppertriang
gdata::lowerTriangle(sim_adjmat, byrow=FALSE) <- uppertriang

# Define parameters
beta <- 0.01
gamma <- 0.00165
rho <- 0.005
epsilon <- 0.01
n_patches <- 9
max_time <- 365*10
dt <- 1
amplitude <- 1.0

# Initialize state variables
S <- c(rbind(rpois(5,500), rpois(5,1000)))
S <- S[-10]
I <- rpois(n_patches, 5)
R <- rep(0, n_patches)

results <- matrix(c(S, I, R, 0), nrow = 1)
colnames(results) <- c(paste0("S", 1:n_patches), paste0("I", 1:n_patches), paste0("R", 1:n_patches), "time")

resultsIncidence <- matrix(0, nrow = 1, ncol = n_patches)

time <- 0

while (time < max_time & any(I > 0)) {

  dS <- rep(0, n_patches)
  dI <- rep(0, n_patches)
  dR <- rep(0, n_patches)
  ObservedIncidence<- rep(0, n_patches)

  for (p in 1:n_patches) {

    #Intra & Inter patch infections
    neighbors<- which(sim_adjmat[p, ] == 1)
    pop_neighbors<- sum(c(S[neighbors],I[neighbors],R[neighbors]))
    if (S[p] > 0 & I[p] > 0 & any(I[neighbors] > 0)) {
      prob_inf <- 1 - exp(-(beta * (1 + amplitude * sin(2 * pi * (time+dt)/365))) * (((1 - epsilon * length(neighbors)) * I[p] + epsilon * sum(I[neighbors]))/((1 - epsilon * length(neighbors)) * (S[p]+I[p]+R[p]) + epsilon * pop_neighbors)) * dt)
      new_infections <- rbinom(1, S[p], prob = min(1, prob_inf))
      dS[p] <- dS[p] - new_infections
      dI[p] <- dI[p] + new_infections
    }

    # Recovery
    if (I[p] > 0) {
      prob_recovery <- 1 - exp(-gamma * dt)
      new_recoveries <- rbinom(1, I[p], prob = min(1, prob_recovery))
      ObservedIncidence[p] <- new_recoveries
      dI[p] <- dI[p] - new_recoveries
      dR[p] <- dR[p] + new_recoveries
    }

    # Immunity loss
    if (R[p] > 0) {
      prob_immunity_loss <- 1 - exp(-rho * dt)
      lost_immunity <- rbinom(1, R[p], prob = min(1, prob_immunity_loss))
      dR[p] <- dR[p] - lost_immunity
      dS[p] <- dS[p] + lost_immunity
    }
  }

  # Update state variables
  S <- S + dS
  I <- I + dI
  R <- R + dR

  time <- time + dt

  results <- rbind(results, c(S, I, R, time))
  resultsIncidence <- rbind(resultsIncidence, c(ObservedIncidence))
  colnames(resultsIncidence) <- paste0("subpop", 1:n_patches)
}

# Plot results
plot(0, type = "n", xlim = c(0, max_time+1), ylim = c(0, max(results[,10:18])), xlab = "Time [days]", ylab = "Population", main ="", cex.lab = 1.8, cex.axis = 1.7)

labels <- c()
Linetypes <- rep(c("dotted", "dashed", "dotdash", "longdash", "twodash"), length.out = n_patches)
colors <- rep(c("blue", "red"), length.out = n_patches)

for (i in 10:(2*n_patches)) {
  lines(results[,ncol(results)], results[,i], col = colors[i-n_patches], lty = Linetypes[i-n_patches])
  labels <- c(labels, colnames(results)[i])
}
#legend("topleft", legend = labels, col = colors, lty = Linetypes, cex = 0.9)
grid()


# Plot daily incidence data
plot(0, type = "n", xlim = c(0, max_time+1), ylim = c(0, max(resultsIncidence)), xlab = "Time [days]", ylab = "Population", main ="", cex.lab = 1.8, cex.axis = 1.7)

labels <- c()
Linetypes <- rep(c("dotted", "dashed", "dotdash", "longdash", "twodash"), length.out = n_patches)
colors <- rep(c("blue", "red"), length.out = n_patches)

for (i in 1:n_patches) {
  lines(1:nrow(resultsIncidence), resultsIncidence[,i], col = colors[i], lty = Linetypes[i])
  labels <- c(labels, colnames(resultsIncidence)[i])
}
#legend("topleft", legend = labels, col = colors, lty = Linetypes, cex = 0.9)
grid()


summingfunction<- function(incidencedata){
  incidencedata<- incidencedata[-1,]
  approxYear<- ceiling(nrow(incidencedata)/365)
  nmonths<- approxYear*12
  summedresult<- matrix(NA, nrow = nmonths, ncol = ncol(incidencedata))
  dummy<- 1
  for(i in 1:nmonths){
    currentmonth<- i%%12
    if(currentmonth %in% c(1,3,5,7,8,10,0)){
      monthlength<- 31
    }else if(currentmonth %in% c(4,6,9,11)){
      monthlength<- 30
    }else{
      monthlength<- 28
    }
    dummy<- dummy + monthlength
    summedresult[i, ]<- colSums(incidencedata[(dummy-monthlength):(dummy-1), ])
  }
  colnames(summedresult) <- paste0("subpop", 1:ncol(summedresult))
  return(summedresult)
}

monthlyincidence<- summingfunction(resultsIncidence)

# Plot monthly incidence
plot(0, type = "n", xaxt="n", xlim = c(0, nrow(monthlyincidence)), ylim = c(0, max(monthlyincidence)), xlab = "Time [month/year]", ylab = "Population", main ="", cex.lab = 1.8, cex.axis = 1.7)

labels <- c()
Linetypes <- rep(c("dotted", "dashed", "dotdash", "longdash", "twodash"), length.out = n_patches)
colors <- rep(c("blue", "red"), length.out = n_patches)

for (i in 1:n_patches) {
  lines(1:nrow(monthlyincidence), monthlyincidence[,i], col = colors[i], lty = Linetypes[i])
  labels <- c(labels, colnames(monthlyincidence)[i])
}
#legend("topleft", legend = labels, col = colors, lty = Linetypes, cex = 0.9)
#years<- 2013:2019
#axis(1, at = seq(1, 84, by = 12), labels = years, cex.axis = 1.8)
grid()
