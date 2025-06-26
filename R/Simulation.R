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
py<- params[["py"]]
n_patches<- nrow(sim_adjmat)
max_time<- params[["max_time"]]
dt<- params[["dt"]]
amplitude<- params[["amplitude"]]

#set.seed(1234)
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
      M_matrix[i, j] <- (beta * (1 - py * length(neighbors)) /
                           ((1 - py * length(neighbors)) * N[i] + py * pop_neighbors))
    } else if (j %in% neighbors) {
      # Transmission from neighboring patches
      M_matrix[i, j] <- (beta * py  /
                           ((1 - py * length(neighbors)) * N[i] + py * pop_neighbors))
    }
  }
}


S<- c(round(solve(M_matrix) %*% rep(gamma, n_patches)))
I<- c(round((rho/(gamma+rho))*(N-S)))
R<- N-S-I

S<- N-5
I<- rep(5, n_patches)
R<- rep(0, n_patches)

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
      prob_inf <- 1 - exp(-(beta * (1 + amplitude * sin(2 * pi * (time+dt)/365))) * (((1 - py * length(neighbors)) * I[p] + py * sum(I[neighbors]))/((1 - py * length(neighbors)) * (S[p]+I[p]+R[p]) + py * pop_neighbors)) * dt)
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

par(mfrow=c(1,3))

# Plot results
plot(0, type = "n", xlim = c(0, max_time+1), ylim = c(0, max(results[,10:18])), xlab = "Time [day]", ylab = "Population", main ="Prevalence", cex.lab = 1.7, cex.axis = 1.7, cex.main=2.0)

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
plot(0, type = "n", xlim = c(0, max_time+1), ylim = c(0, max(resultsIncidence)), xlab = "Time [day]", ylab = "Population", main ="Daily incidence", cex.lab = 1.7, cex.axis = 1.7, cex.main=2.0)

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
  mindex<- numeric(nmonths)
  for(i in 1:nmonths){
    currentmonth<- (i - 1) %% 12 + 1
    if(currentmonth %in% c(1,3,5,7,8,10,12)){
      monthlength<- 31
    }else if(currentmonth %in% c(4,6,9,11)){
      monthlength<- 30
    }else{
      monthlength<- 28
    }
    dummy<- dummy + monthlength
    mindex[i]<- dummy-1
    summedresult[i, ]<- colSums(incidencedata[(dummy-monthlength):(dummy-1), ])
  }
  colnames(summedresult) <- paste0("subpop", 1:ncol(summedresult))
  return(list(summedresult, mindex))
}

summed<- summingfunction(resultsIncidence)
monthlyincidence<- summed[[1]]
#mindex<- summed[[2]]
#noinit<- results[-1,]
#monthlysusceptibles<- noinit[mindex,1:n_patches]

# Plot monthly incidence
plot(0, type = "n", xlim = c(0, nrow(monthlyincidence)), ylim = c(0, max(monthlyincidence)), xlab = "Time [month]", ylab = "Population", main ="Monthly incidence", cex.lab = 1.7, cex.axis = 1.7, cex.main=2.0)

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

#Phase plot
par(mfrow=c(3,3))
for(j in 1:n_patches){
  plot(results[,18+j]/N[j], results[,9+j]/N[j], xlim=c(0,1), ylim=c(0,1), type = "l", xlab = "R", ylab = "I", main = paste("Subpopulation ", j))
}


n_time <- nrow(results)
#colors <- viridis::viridis(n_time - 1)
seg_len <- floor(n_time/3)
segcolors <- c("blue", "yellow", "red")
par(mfrow = c(3, 3))
for (j in 1:n_patches) {
  I_vals <- results[, 9 + j] / N[j]
  R_vals <- results[, 18 + j] / N[j]
  plot(NA, xlim = c(0, 1), ylim = c(0, 1), xlab = "R", ylab = "I", main = paste("Subpopulation", j))
  grid()
  for (i in 1:(n_time - 1)) {
    segment_color <- if (i <= seg_len) {
      segcolors[1]
    } else if (i <= 2 * seg_len) {
      segcolors[2]
    } else {
      segcolors[3]
    }
    segments(
      x0 = R_vals[i],
      y0 = I_vals[i],
      x1 = R_vals[i + 1],
      y1 = I_vals[i + 1],
      col = segment_color
    )
  }
}


#Model fitting
library(DetectOutbreaks)
y<- t(monthlyincidence[-(1:36),])
e_it<- matrix(rep(N, ncol(y)), nrow = n_patches, ncol = ncol(y), byrow = F)

#Mod0<- DetectOutbreaks::infer(y=y, e_it = e_it, adjmat = sim_adjmat, Model = 0, verbose = T)
#Mod1<- DetectOutbreaks::infer(y=y, e_it = e_it, adjmat = sim_adjmat, Model = 1, verbose = T)
#Mod2<- DetectOutbreaks::infer(y=y, e_it = e_it, adjmat = sim_adjmat, Model = 2, verbose = T)
#Mod3<- DetectOutbreaks::infer(y=y, e_it = e_it, adjmat = sim_adjmat, Model = 3, verbose = T)
#Mod4<- DetectOutbreaks::infer(y=y, e_it = e_it, adjmat = sim_adjmat, Model = 4, verbose = T)
#Mod5<- DetectOutbreaks::infer(y=y, e_it = e_it, adjmat = sim_adjmat, Model = 5, verbose = T)
#Mod6<- DetectOutbreaks::infer(y=y, e_it = e_it, adjmat = sim_adjmat, Model = 6, verbose = T)
#Mod7<- DetectOutbreaks::infer(y=y, e_it = e_it, adjmat = sim_adjmat, Model = 7, verbose = T)

#simdata<- list(y, e_it)
#allmods<- list(Mod0,Mod1,Mod2,Mod3,Mod4,Mod5,Mod6,Mod7)
#Posteriorpredictive(allmods, simdata,sim_adjmat)


#applicationCases <- readRDS("C:/Users/Matthew Adeoye/Downloads/applicationCases.RDS")
#applicationPopulation <- readRDS("C:/Users/Matthew Adeoye/Downloads/applicationPopulation.RDS")
#adjmat <- readRDS("C:/Users/Matthew Adeoye/Downloads/applicationAdjMat.RDS")
#y<- applicationCases[,-(1:168)]
#e_it<- applicationPopulation[,-(1:168)]
#y<- ApplicationCounts
#e_it<- ApplicationPopulation
#sim_adjmat<- ApplicationAdjMat

R<- -1 * sim_adjmat
diag(R)<- -rowSums(R, na.rm = T)
rankdef<- nrow(R)-qr(R)$rank
nstate<- 2
ndept<- nrow(y)
time<- ncol(y)

original.y<- y

#flag missing data for inference
y<- ifelse(is.na(y), -1, y)

RW1PrecMat<- matrix(0, nrow=12, ncol=12)
RW1PrecMat[1, ]<- c(2,-1, rep(0, 12-3), -1)
RW1PrecMat[2, ]<- c(-1,2,-1, rep(0, 12-3))
RW1PrecMat[3, ]<- c(0, -1,2,-1, rep(0, 12-4))
RW1PrecMat[(12-1), ]<- c(rep(0, 12-3), -1,2,-1)
RW1PrecMat[12, ]<- c(-1, rep(0, 12-3), -1, 2)
for(i in 3:(12-3)){
  RW1PrecMat[i+1, ((i):(i+2))]<- c(-1,2,-1)
}
strs<- RW1PrecMat

Model<- 0
npar<- 0

design_matrix_func <- get(paste0("DesignMatrixModel", Model))
z_it <- design_matrix_func(original.y, sim_adjmat)[[1]]
z_it2 <- design_matrix_func(original.y, sim_adjmat)[[2]]

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel :: detectCores())
app_data = list(ndept=ndept, time=time, nstate=nstate, rankdef=rankdef, y=y, e_it=e_it, R=R,
                SMat = strs, Model = Model, npar = npar, z_it = z_it, z_it2 = z_it2)
parameters =c("G12", "G21", "kappa_u" , "kappa_r", "kappa_s", "u", "uconstrained", "r", "sraw", "s", "B")

nchains =4
n_warmups =1500
n_iter =3000
n_thin =1

initials <- list(G12 = 0.1, G21 = 0.3, u = rep(0, ndept-1), r = rep(0, time), sraw = rep(0, 11), kappa_u=20, kappa_r=20, kappa_s=20, B=rep(0.1, npar))
initials_list <- lapply(1:nchains, function(x) initials)

#nuts_fit0 = stan(file = "C:/Users/Matthew Adeoye/Documents/PhD Statistics/Year 3/SpatMet/inst/stan/modelfile.stan" , # Stan program
#                data = app_data,
#                pars = parameters,
#                init = initials_list,
#                chains = nchains,
#                warmup = n_warmups,
#                iter = n_iter,
#                thin = n_thin,
#                seed =13219)

#model_code <- "C:/Users/Matthew Adeoye/Documents/PhD Statistics/Year 3/SpatMet/inst/stan/modelfile.stan"
#model_env <- cmdstanr::cmdstan_model(stan_file = model_code, compile = TRUE)
#cmdstanrfit<- model_env$sample(data = list(ndept=ndept, time=time, nstate=nstate, rankdef=rankdef, y=y, e_it=e_it, R=R,
#                               SMat = strs, Model = Model, npar = npar, z_it = z_it, z_it2 = z_it2),
#                   init = initials_list, chains = nchains, iter_warmup = 1500,
#                   iter_sampling = 1500, parallel_chains = nchains,
#                   seed=1234, adapt_delta = 0.90)

#bridgesampling::bridge_sampler(samples=nuts_fit0, maxiter=2000)
#DetectOutbreaks::ModelEvidence(y=y, e_it = e_it, adjmat = sim_adjmat, Model = 0,inf.object = Mod0, num_samples = 50000)
#mod0=-2881.156, nutsfit0=-2901.412
#mod3=-2884.392, nutsfit3=-2950.459

#appData7<- forecast(Model=7, inf.object=CorrectedGPU_7, forecastlength=24, y=applicationCases[,-(1:168)], popdata=applicationPopulation[,-(1:168)], adj.matrix=adjmat)
#DetectOutbreaks::sim.plot(cbind(applicationCases[,-(1:168)],appData7[[1]]))



#ISestimates<- numeric(20);Bridgeestimates<- numeric(20);for(i in 1:20){
#  cFit<- DetectOutbreaks::infer(y=y, e_it = e_it, adjmat = sim_adjmat, Model = 0, verbose = T)
#  ISestimates[i]<- DetectOutbreaks::ModelEvidence(y=y, e_it = e_it, adjmat = sim_adjmat, Model = 0,inf.object = cFit, num_samples = 50000)
#  Bridgeestimates[i]<- ModelEvidence.Bridge3(y=y, e_it = e_it, adjmat = sim_adjmat, Model = 0,inf.object = cFit, num_samples = 50000)
#}
