functions{

  int intPower(int a, int b){
    int res = 1;
    for (i in 1:b) {
    res = a * res;
    }
    return(res);
}

  vector logVecMatMult(vector logV, matrix logM){
    int S = dims(logV)[1];
    vector[S] res;
    for(s in 1:S){
      res[s] = log_sum_exp(logV + logM[, s]);
    }
    return(res);
  }

  matrix JointTransitionMatrix(matrix gamma, int K){
    int S = intPower(2, K);
    matrix[S, S] Gamma;
      for(a in 0:(S - 1)){
      for(b in 0:(S - 1)){
      real prob = 1;
      for(k in 1:K){
       int from_k = (a %/% intPower(2, (k - 1))) % 2;
       int to_k = (b %/% intPower(2, (k - 1))) % 2;
        prob = prob * gamma[from_k + 1, to_k + 1];
      }
      Gamma[a + 1, b + 1] = prob;
    }
  }
  return(Gamma);
  }

  //Stationary distribution
  vector stationarydist(matrix Gamma){
    int ncol = dims(Gamma)[2];

    matrix[ncol, ncol] mT = transpose(Gamma);

    complex_vector[ncol] E_values = eigenvalues(mT); vector[ncol] NE_values = get_real(E_values);
    complex_matrix[ncol, ncol] E_vectors = eigenvectors(mT);

    real min_dist = positive_infinity();
    int index = 1;
    for(i in 1:ncol){
    real dist = abs(NE_values[i] - 1);
    if(dist < min_dist){
    min_dist = dist;
    index = i;
    }
  }
    complex_vector[ncol] stationary_distribution = E_vectors[, index];
    vector[ncol] Nstationary_distribution = get_real(stationary_distribution);
    Nstationary_distribution /= sum(Nstationary_distribution);
    return(Nstationary_distribution);
}

  //colsums function
  array[] int colSums(array[,] int M){
    int ncol = dims(M)[2];
    array[ncol] int sums;

    for(i in 1:ncol){
      sums[i] = sum(M[,i]);
    }
    return(sums);
  }

  // Intrinsic GMRF density
  real IGMRF1(vector uconstrained, real kappa_u, matrix R, int rankdef) {
    int n = rows(R);
    return (((n - rankdef) / 2.0) * (log(kappa_u) - log(2.0 * pi())) - (kappa_u / 2.0) * quad_form(R, uconstrained));
  }

  // Random walk density
  real randomwalk2(vector r, real kappa_r) {
    int time = dims(r)[1];
    real res = 0;
    for (i in 3:time) {
      res += (r[i-2] - (2 * r[i-1]) + r[i])^2;
    }
    return (((time - 2) / 2.0) * log(kappa_r) - (kappa_r / 2.0) * res);
  }

  // Seasonal components' density
    real seasonalComp(vector s, real kappa_s, matrix SMat) {
    int n = rows(SMat);
    return (((n - 1) / 2.0) * (log(kappa_s) - log(2.0 * pi())) - (kappa_s / 2.0) * quad_form(SMat, s));
}

  //Gamma matrix function
    matrix G(real G12, real G21) {
    matrix[2, 2] m;
    m[1, 1] = 1 - G12;
    m[1, 2] = G12;
    m[2, 1] = G21;
    m[2, 2] = 1 - G21;
    return m;
  }

  real my_dpois(int y, real lambda){
  real res = -lambda + y * log(lambda) - log(tgamma(y+1));
  return res;
}

  //multi-strain loglikelihood via forward filtering
  real Stan_Loglikelihood(array[] int y, int ndept, int time, int nstrain, vector a_k, vector r, vector s, vector u, matrix gamma, matrix e_it, vector B, int Model, matrix Bits, int independentChains){

    //Model0
  if(Model == 0){
    real allLoglikelihood = 0;
  for(i in 1:ndept) {
    for(t in 1:time) {
      int month_index = (t - 1) % 12 + 1;
      for(k in 1:nstrain){
      if(y[(k-1) * ndept * time + (i-1) * time + t] == -1){
        allLoglikelihood += 0;
      }else{
        allLoglikelihood += poisson_lpmf(y[(k-1) * ndept * time + (i-1) * time + t] | e_it[i, t] * exp(a_k[k] + r[t] + s[month_index] + u[i]));
        }
      }
    }
  }
      return allLoglikelihood;
}

 //If Model specification is not Model 0.
  else{
    if(independentChains==0){
  int nstate = intPower(2, nstrain);
  matrix[nstate, nstate] jointTPM = JointTransitionMatrix(gamma, nstrain);
  matrix[time, nstate] Alpha;
  vector[nstate] init_density = stationarydist(jointTPM);
  vector[ndept] log_forwards;

  for (i in 1:ndept){
    vector[nstate] prodEmission;
  // Initialization of the first time step for each department
  for(n in 1:nstate){
    prodEmission[n] = 0;
  for(k in 1:nstrain){
    vector[nstrain] newB = rep_vector(0.0, nstrain);
    newB[k] = B[k];
  if(y[(k-1) * ndept * time + (i-1) * time + 1] == -1){
    prodEmission[n] += 0;
    }else{
    prodEmission[n] += poisson_lpmf(y[(k-1) * ndept * time + (i-1) * time + 1] | e_it[i, 1] * exp(a_k[k] + r[1] + s[1] + u[i] + dot_product(newB, Bits[n, ])));
    }
  }
}
  Alpha[1] = transpose(log(init_density) + prodEmission);

  // Dynamic programming loop for the remaining time steps
      for(t in 2:time) {
        int month_index = (t - 1) % 12 + 1;
        for(n in 1:nstate){
          prodEmission[n] = 0;
        for(k in 1:nstrain){
          vector[nstrain] newB = rep_vector(0.0, nstrain);
          newB[k] = B[k];
        if(y[(k-1) * ndept * time + (i-1) * time + t] == -1){
          prodEmission[n] += 0;
          }else{
         prodEmission[n] += poisson_lpmf(y[(k-1) * ndept * time + (i-1) * time + t] | e_it[i, t] * exp(a_k[k] + r[t] + s[month_index] + u[i] + dot_product(newB, Bits[n, ])));
         }
        }
      }
        Alpha[t] = transpose(logVecMatMult(transpose(Alpha[t-1, ]), log(jointTPM)) + prodEmission);
    }
        log_forwards[i] = log_sum_exp(Alpha[time, ]);
  }
      real fullLogLikelihood = sum(log_forwards);
      return fullLogLikelihood;
    }//If independent Markov chains.
else{
  int nstate = 2;
  matrix[time, nstate] Alpha;
  vector[nstate] prodEmission;
  vector[nstate] init_density = stationarydist(gamma);
  vector[ndept] log_forwards;

  for (i in 1:ndept){
    // Initialization of the first time step for each department
    for(n in 1:nstate){
      prodEmission[n] = 0;
    for(k in 1:nstrain){
      if(y[(k-1) * ndept * time + (i-1) * time + 1] == -1){
        prodEmission[n] += 0;
      }else{
          prodEmission[n] += poisson_lpmf(y[(k-1) * ndept * time + (i-1) * time + 1] | e_it[i, 1] * exp(a_k[k] + r[1] + s[1] + u[i] + B[k] * (n-1)));
        }
    }
  }
     Alpha[1] = transpose(log(init_density) + prodEmission);
    // Dynamic programming loop for the remaining time steps
    for(t in 2:time) {
      int month_index = (t - 1) % 12 + 1;
      for(n in 1:nstate){
        prodEmission[n] = 0;
        for(k in 1:nstrain){
          if(y[(k-1) * ndept * time + (i-1) * time + t] == -1){
            prodEmission[n] += 0;
          }else{
            prodEmission[n] += poisson_lpmf(y[(k-1) * ndept * time + (i-1) * time + t] | e_it[i, t] * exp(a_k[k] + r[t] + s[month_index] + u[i] + B[k] * (n-1)));
          }
        }
      }
      Alpha[t] = transpose(logVecMatMult(transpose(Alpha[t-1, ]), log(gamma)) + prodEmission);
    }
    log_forwards[i] = log_sum_exp(Alpha[time, ]);
  }
  real fullLogLikelihood = sum(log_forwards);
  return fullLogLikelihood;
}
   return 0;
    }
  }
}

data {
  int<lower=1> ndept;                 // Number of departments
  int<lower=1> time;                  // Time
  int<lower=1> nstate;                // Number of states
  int<lower=1> rankdef;               // Rank deficiency of structure matrix (R)
  int<lower=1> nstrain;               // Number of strains
  int<lower=0> npar;                  // Number of regression parameters
  array[ndept*time*nstrain] int y;  // data matrix
  matrix[ndept, time] e_it;           // initial Susceptibles
  matrix[ndept, ndept] R;             // Structure matrix (IGMRF1)
  int<lower=0, upper=1> Model;        // Model's functional form
  matrix[12, 12] SMat;                //Structure matrix (seasonal_comp)
  matrix[nstate, nstrain] Bits;       //Bits matrix
  int<lower=0> independentChains;     // 0=> False, 1=> True
}

parameters {
  real<lower=0, upper=1> G12;            // transition to hyperendemic
  real<lower=0, upper=1> G21;            // transition to endemic
  real<lower=0> kappa_u;                 // spatial precision parameter
  real<lower=0> kappa_r;                 // trend precision parameter
  real<lower=0> kappa_s;                 // seasonal precision parameter
  vector[ndept-1] u;                     // Spatial components
  vector[time] rraw;                     // Trend components
  vector[11] sraw;                       // cyclic seasonal components
  vector<lower=0>[npar] B;            // autoregressive parameters
  vector[nstrain] a_k;                   // background intercepts
}

transformed parameters {
  real sumC = sum(u[1:(ndept-1)]);
  vector[ndept] uconstrained;
  uconstrained = append_row(u[1:(ndept-1)], -sumC);
  real sumS = sum(sraw[1:11]);
  vector[12] s;
  s = append_row(sraw[1:11], -sumS);

  vector[time] r;
  real MeanR = mean(rraw);
  for(i in 1:time){
  r[i] = rraw[i] - MeanR;
  }
}

model {
  target += beta_lpdf(G12 | 2, 2);
  target += beta_lpdf(G12 | 2, 2);
  target += gamma_lpdf(kappa_u | 1, 0.01);
  target += gamma_lpdf(kappa_r | 1, 0.0001);
  target += gamma_lpdf(kappa_s | 1, 0.001);
  if(Model==1){
  target += gamma_lpdf(B | 2, 2);
  }
  target += IGMRF1(uconstrained, kappa_u, R, rankdef);
  target += randomwalk2(r, kappa_r);
  target += seasonalComp(s, kappa_s, SMat);

  // Likelihood
  target += Stan_Loglikelihood(y, ndept, time, nstrain, a_k, r, s, uconstrained, G(G12, G21), e_it, B, Model, Bits, independentChains);
}

generated quantities{
  real log_lik = Stan_Loglikelihood(y, ndept, time, nstrain, a_k, r, s, uconstrained, G(G12, G21), e_it, B, Model, Bits, independentChains);
  real state1_stationary_dist = stationarydist(G(G12, G21))[2];
  vector[nstate] stationaryDistribution = stationarydist(JointTransitionMatrix(G(G12, G21), nstrain));
}
