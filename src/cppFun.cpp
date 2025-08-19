#include <Rcpp.h>
#include <RcppEigen.h>
#include <cmath>
#include <algorithm>

using namespace Rcpp;

// [[Rcpp::depends(RcppEigen)]

// [[Rcpp::export]]
double logSumExp_cpp(NumericVector x) {
  double max_val = max(x);
  return log(sum(exp(x - max_val))) + max_val;
}

// [[Rcpp::export]]
int intPower(int a, int b){
  int res = 1;
  for (int i = 0; i < b; ++i) {
    res = a * res;
  }
  return res;
}

// [[Rcpp::export]]
NumericVector logVecMatMult(NumericVector logV, NumericMatrix M) {
  int S = logV.size();
  NumericVector res(S);
  for (int p = 0; p < S; ++p) {
    NumericVector temp(S);
    for (int s = 0; s < S; ++s) {
      temp[s] = logV[s] + log(M(s, p));
    }
    res[p] = logSumExp_cpp(temp);
  }
  return res;
}

// [[Rcpp::export]]
NumericMatrix JointTransitionMatrix_cpp(NumericMatrix gamma, int K) {
  int S = intPower(2, K);
  NumericMatrix jointGamma(S, S);
  for (int a = 0; a < S; ++a) {
    for (int b = 0; b < S; ++b) {
      double prob = 1.0;
      for (int k = 0; k < K; ++k) {
        int from_k = (a / intPower(2, k)) % 2;
        int to_k   = (b / intPower(2, k)) % 2;
        prob *= gamma(from_k, to_k);
      }
      jointGamma(a, b) = prob;
    }
  }
  return jointGamma;
}

// [[Rcpp::export]]
NumericVector state_dist_cpp2(NumericMatrix Gamma) {
  int n = Gamma.ncol();
  Eigen::MatrixXd m(n, n);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      m(i, j) = Gamma(i, j);
    }
  }

  Eigen::MatrixXd mT = m.transpose();

  // Eigen decomposition
  Eigen::EigenSolver<Eigen::MatrixXd> es(mT);
  Eigen::VectorXd eigenvalues = es.eigenvalues().real();
  Eigen::MatrixXd eigenvectors = es.eigenvectors().real();

  // Index of eigenvalue close to 1
  int index = 0;
  double min_diff = std::abs(eigenvalues(0) - 1.0);
  for (int i = 1; i < eigenvalues.size(); ++i) {
    double diff = std::abs(eigenvalues(i) - 1.0);
    if (diff < min_diff) {
      min_diff = diff;
      index = i;
    }
  }

  // corresponding eigenvector
  Eigen::VectorXd stationary_distribution = eigenvectors.col(index);

  // Normalize stationary distribution
  stationary_distribution = stationary_distribution / stationary_distribution.sum();

  // Convert result to NumericVector
  NumericVector result(stationary_distribution.data(),
                       stationary_distribution.data() + stationary_distribution.size());

  return result;
}

// [[Rcpp::export]]
NumericVector state_dist_cpp(double G12, double G21) {
  Eigen::MatrixXd m(2, 2);
  m(0, 0) = 1 - G12;
  m(0, 1) = G12;
  m(1, 0) = G21;
  m(1, 1) = 1 - G21;

  Eigen::MatrixXd b = m.transpose();

  // Eigen decomposition
  Eigen::EigenSolver<Eigen::MatrixXd> es(b);
  Eigen::VectorXd eigenvalues = es.eigenvalues().real();
  Eigen::MatrixXd eigenvectors = es.eigenvectors().real();

  // Index of eigenvalue close to 1
  int index = 0;
  double min_diff = std::abs(eigenvalues(0) - 1.0);
  for (int i = 1; i < eigenvalues.size(); ++i) {
    double diff = std::abs(eigenvalues(i) - 1.0);
    if (diff < min_diff) {
      min_diff = diff;
      index = i;
    }
  }

  // corresponding eigenvector
  Eigen::VectorXd stationary_distribution = eigenvectors.col(index);

  // Normalize stationary distribution
  stationary_distribution = stationary_distribution / stationary_distribution.sum();

  // Convert result to NumericVector
  NumericVector result(stationary_distribution.data(),
                       stationary_distribution.data() + stationary_distribution.size());

  return result;
}

// [[Rcpp::export]]
double randomwalk2_cpp(NumericVector componentR, double PrecisionR){
  int time = componentR.size();
  double Sumres = 0;
  for(int i = 2; i < time; ++i) {
    double res = pow(componentR[i - 2] - (2 * componentR[i - 1]) + componentR[i], 2);
    Sumres += res;
  }
  return (time - 2) / 2.0 * log(PrecisionR) - PrecisionR / 2.0 * Sumres;
}

// [[Rcpp::export]]
double dotproduct_cpp(NumericVector v1, NumericVector v2) {
  int lengthV1 = v1.size();
  double result = 0;
  for (int i = 0; i < lengthV1 ; ++i) {
    result += v1[i] * v2[i];
  }
  return result;
}

// [[Rcpp::export]]
double quadform_cpp(NumericVector v1, NumericMatrix z) {
  int lengthV1 = v1.size();
  NumericVector newvec(lengthV1);
  for (int i = 0; i < lengthV1 ; ++i) {
    double result = 0;
    for(int j = 0; j < lengthV1 ; ++j){
      result += v1[j] * z(j, i);
    }
    newvec[i] = result;
  }
  double finalRes = dotproduct_cpp(newvec, v1);
  return finalRes;
}

// [[Rcpp::export]]
double seasonalComp2_cpp(NumericVector x, double y, NumericMatrix z){
  int n = z.nrow();
  double sumC = 0;
  NumericVector x_new(n-1);
  for(int i = 0; i < n - 1; ++i) {
    x_new[i] = x[i];
    sumC += x[i];
  }
  x_new.push_back(-sumC);
  double result = (n - 1) / 2.0 * (log(y) - log(2 * M_PI)) - y / 2.0 * quadform_cpp(x_new, z);
  return result;
}

// [[Rcpp::export]]
double logIGMRF1_cpp(NumericVector x, double y, NumericMatrix z, int rankdef){
  int n = z.nrow();
  double sumC = 0;
  NumericVector x_new(n-1);
  for(int i = 0; i < n - 1; ++i) {
    x_new[i] = x[i];
    sumC += x[i];
  }
  x_new.push_back(-sumC);
  double result = (n - rankdef) / 2.0 * (log(y) - log(2 * M_PI)) - y / 2.0 * quadform_cpp(x_new, z);
  return result;
}

//Fast PMF computation
// [[Rcpp::export]]
double dpois_cpp(int y, double lambda){
  double res = -lambda + y * log(lambda) - log(Rcpp::internal::factorial(y));
  return res;
}


// Computing log likelihood
// [[Rcpp::export]]
double GeneralLoglikelihood_cpp(NumericMatrix y, NumericVector r, NumericVector s, NumericVector u,
                                NumericMatrix Gamma, NumericMatrix e_it, NumericVector B, int model,
                                NumericMatrix z_it, NumericMatrix z_it2) {
  int ndept = y.nrow();
  int nstate = Gamma.ncol();
  int time = y.ncol();
  double gamma_11 = log(Gamma(0, 0));
  double gamma_12 = log(Gamma(0, 1));
  double gamma_21 = log(Gamma(1, 0));
  double gamma_22 = log(Gamma(1, 1));

  NumericMatrix log_forward_probs(ndept, nstate);

  // Model 0
  if(model == 0) {
    NumericMatrix log_likelihood(ndept, time);

    for(int i = 0; i < ndept; ++i) {
      for(int t = 0; t < time; ++t) {
        if(y(i, t) == -1){
          log_likelihood(i, t) = 0;
        }else{
          log_likelihood(i, t) = R::dpois(y(i, t), e_it(i, t) * exp(r[t] + s[t] + u[i]), true);
        }
      }
    }

    double full_log_likelihood = sum(log_likelihood);
    return full_log_likelihood;
  }

  // Model 1, 2, 4, 5 or 7
  if(model == 1 || model == 2 || model == 4 || model == 5 || model == 7) {
    NumericVector init_density = state_dist_cpp(Gamma(0, 1), Gamma(1, 0));
    init_density = log(init_density);
    NumericMatrix Alphas(time, nstate);
    NumericMatrix alpha(time, nstate);
    NumericMatrix beta(time, nstate);
    NumericMatrix log_forward_probs(ndept, nstate);
    NumericVector rowlogsumexp(ndept);

    for(int i = 0; i < ndept; ++i) {
      if(y(i, 0) == -1){
        alpha(0, 0) = init_density[0];
        alpha(0, 1) = 0;
        beta(0, 0) = 0;
        beta(0, 1) = init_density[1];
      }else{
        alpha(0, 0) = init_density[0] + R::dpois(y(i, 0), e_it(i, 0) * exp(r[0] + s[0] + u[i]), true);
        alpha(0, 1) = 0;
        beta(0, 0) = 0;
        beta(0, 1) = init_density[1] + R::dpois(y(i, 0), e_it(i, 0) * exp(r[0] + s[0] + u[i] + B[0] * z_it(i, 0)), true);
      }
      Alphas(0, 0) = alpha(0, 0);
      Alphas(0, 1) = beta(0, 1);

      for(int t = 1; t < time; ++t) {
        if(y(i, t) == -1){
          alpha(t, 0) = Alphas(t-1, 0) + gamma_11;
          alpha(t, 1) = Alphas(t-1, 1) + gamma_21;
          beta(t, 0) = Alphas(t-1, 0) + gamma_12;
          beta(t, 1) = Alphas(t-1, 1) + gamma_22;
          Alphas(t, 0) = logSumExp_cpp(alpha(t, _));
          Alphas(t, 1) = logSumExp_cpp(beta(t, _));
        }else{
          alpha(t, 0) = Alphas(t-1, 0) + gamma_11 + R::dpois(y(i, t), e_it(i, t) * exp(r[t] + s[t] + u[i]), true);
          alpha(t, 1) = Alphas(t-1, 1) + gamma_21 + R::dpois(y(i, t), e_it(i, t) * exp(r[t] + s[t] + u[i]), true);
          beta(t, 0) = Alphas(t-1, 0) + gamma_12 + R::dpois(y(i, t), e_it(i, t) * exp(r[t] + s[t] + u[i] + B[0] * z_it(i, t)), true);
          beta(t, 1) = Alphas(t-1, 1) + gamma_22 + R::dpois(y(i, t), e_it(i, t) * exp(r[t] + s[t] + u[i] + B[0] * z_it(i, t)), true);
          Alphas(t, 0) = logSumExp_cpp(alpha(t, _));
          Alphas(t, 1) = logSumExp_cpp(beta(t, _));
        }
      }

      log_forward_probs(i, 0) = logSumExp_cpp(alpha(time-1, _));
      log_forward_probs(i, 1) = logSumExp_cpp(beta(time-1, _));
    }

    for(int i = 0; i < ndept; ++i) {
      rowlogsumexp[i] = logSumExp_cpp(log_forward_probs(i, _));
    }

    double full_log_likelihood = sum(rowlogsumexp);
    return full_log_likelihood;
  }

  // Model 3 or 6
  if(model == 3 || model == 6) {
    NumericVector init_density = state_dist_cpp(Gamma(0, 1), Gamma(1, 0));
    init_density = log(init_density);
    NumericMatrix Alphas(time, nstate);
    NumericMatrix alpha(time, nstate);
    NumericMatrix beta(time, nstate);
    NumericMatrix log_forward_probs(ndept, nstate);
    NumericVector rowlogsumexp(ndept);

    for(int i = 0; i < ndept; ++i) {
      if(y(i, 0) == -1){
        alpha(0, 0) = init_density[0];
        alpha(0, 1) = 0;
        beta(0, 0) = 0;
        beta(0, 1) = init_density[1];
      }else{
        alpha(0, 0) = init_density[0] + R::dpois(y(i, 0), e_it(i, 0) * exp(r[0] + s[0] + u[i]), true);
        alpha(0, 1) = 0;
        beta(0, 0) = 0;
        beta(0, 1) = init_density[1] + R::dpois(y(i, 0), e_it(i, 0) * exp(r[0] + s[0] + u[i] + B[0] * z_it(i, 0) + B[1] * z_it2(i, 0)), true);
      }
      Alphas(0, 0) = alpha(0, 0);
      Alphas(0, 1) = beta(0, 1);

      for(int t = 1; t < time; ++t) {
        if(y(i, t) == -1){
          alpha(t, 0) = Alphas(t-1, 0) + gamma_11;
          alpha(t, 1) = Alphas(t-1, 1) + gamma_21;
          beta(t, 0) = Alphas(t-1, 0) + gamma_12;
          beta(t, 1) = Alphas(t-1, 1) + gamma_22;
          Alphas(t, 0) = logSumExp_cpp(alpha(t, _));
          Alphas(t, 1) = logSumExp_cpp(beta(t, _));
        }else{
          alpha(t, 0) = Alphas(t-1, 0) + gamma_11 + R::dpois(y(i, t), e_it(i, t) * exp(r[t] + s[t] + u[i]), true);
          alpha(t, 1) = Alphas(t-1, 1) + gamma_21 + R::dpois(y(i, t), e_it(i, t) * exp(r[t] + s[t] + u[i]), true);
          beta(t, 0) = Alphas(t-1, 0) + gamma_12 + R::dpois(y(i, t), e_it(i, t) * exp(r[t] + s[t] + u[i] + B[0] * z_it(i, t) + B[1] * z_it2(i, t)), true);
          beta(t, 1) = Alphas(t-1, 1) + gamma_22 + R::dpois(y(i, t), e_it(i, t) * exp(r[t] + s[t] + u[i] + B[0] * z_it(i, t) + B[1] * z_it2(i, t)), true);
          Alphas(t, 0) = logSumExp_cpp(alpha(t, _));
          Alphas(t, 1) = logSumExp_cpp(beta(t, _));
        }
      }

      log_forward_probs(i, 0) = logSumExp_cpp(alpha(time-1, _));
      log_forward_probs(i, 1) = logSumExp_cpp(beta(time-1, _));
    }

    for(int i = 0; i < ndept; ++i) {
      rowlogsumexp[i] = logSumExp_cpp(log_forward_probs(i, _));
    }

    double full_log_likelihood = sum(rowlogsumexp);
    return full_log_likelihood;
  }

  return R_NegInf;
}


// Computing log likelihood using cyclic seasonal components
// [[Rcpp::export]]
double GeneralLoglikelihood_cpp2(NumericMatrix y, NumericVector r, NumericVector s, NumericVector u,
                                 NumericMatrix Gamma, NumericMatrix e_it, NumericVector B, int model,
                                 NumericMatrix z_it, NumericMatrix z_it2) {
  int ndept = y.nrow();
  int nstate = Gamma.ncol();
  int time = y.ncol();
  double gamma_11 = log(Gamma(0, 0));
  double gamma_12 = log(Gamma(0, 1));
  double gamma_21 = log(Gamma(1, 0));
  double gamma_22 = log(Gamma(1, 1));

  NumericMatrix log_forward_probs(ndept, nstate);

  // Model 0
  if(model == 0) {
    NumericMatrix log_likelihood(ndept, time);

    for(int i = 0; i < ndept; ++i) {
      for(int t = 0; t < time; ++t) {
        if(y(i, t) == -1){
          log_likelihood(i, t) = 0;
        }else{
          int month_index = t % 12;
          log_likelihood(i, t) = R::dpois(y(i, t), e_it(i, t) * exp(r[t] + s[month_index] + u[i]), true);
        }
      }
    }

    double full_log_likelihood = sum(log_likelihood);
    return full_log_likelihood;
  }

  // Model 1, 2, 4, 5 or 7
  if(model == 1 || model == 2 || model == 4 || model == 5 || model == 7) {
    NumericVector init_density = state_dist_cpp(Gamma(0, 1), Gamma(1, 0));
    init_density = log(init_density);
    NumericMatrix Alphas(time, nstate);
    NumericMatrix alpha(time, nstate);
    NumericMatrix beta(time, nstate);
    NumericMatrix log_forward_probs(ndept, nstate);
    NumericVector rowlogsumexp(ndept);

    for(int i = 0; i < ndept; ++i) {
      if(y(i, 0) == -1){
        alpha(0, 0) = init_density[0];
        alpha(0, 1) = 0;
        beta(0, 0) = 0;
        beta(0, 1) = init_density[1];
      }else{
        alpha(0, 0) = init_density[0] + R::dpois(y(i, 0), e_it(i, 0) * exp(r[0] + s[0] + u[i]), true);
        alpha(0, 1) = 0;
        beta(0, 0) = 0;
        beta(0, 1) = init_density[1] + R::dpois(y(i, 0), e_it(i, 0) * exp(r[0] + s[0] + u[i] + B[0] * z_it(i, 0)), true);
      }
      Alphas(0, 0) = alpha(0, 0);
      Alphas(0, 1) = beta(0, 1);

      for(int t = 1; t < time; ++t) {
        if(y(i, t) == -1){
          alpha(t, 0) = Alphas(t-1, 0) + gamma_11;
          alpha(t, 1) = Alphas(t-1, 1) + gamma_21;
          beta(t, 0) = Alphas(t-1, 0) + gamma_12;
          beta(t, 1) = Alphas(t-1, 1) + gamma_22;
          Alphas(t, 0) = logSumExp_cpp(alpha(t, _));
          Alphas(t, 1) = logSumExp_cpp(beta(t, _));
        }else{
          int month_index = t % 12;
          alpha(t, 0) = Alphas(t-1, 0) + gamma_11 + R::dpois(y(i, t), e_it(i, t) * exp(r[t] + s[month_index] + u[i]), true);
          alpha(t, 1) = Alphas(t-1, 1) + gamma_21 + R::dpois(y(i, t), e_it(i, t) * exp(r[t] + s[month_index] + u[i]), true);
          beta(t, 0) = Alphas(t-1, 0) + gamma_12 + R::dpois(y(i, t), e_it(i, t) * exp(r[t] + s[month_index] + u[i] + B[0] * z_it(i, t)), true);
          beta(t, 1) = Alphas(t-1, 1) + gamma_22 + R::dpois(y(i, t), e_it(i, t) * exp(r[t] + s[month_index] + u[i] + B[0] * z_it(i, t)), true);
          Alphas(t, 0) = logSumExp_cpp(alpha(t, _));
          Alphas(t, 1) = logSumExp_cpp(beta(t, _));
        }
      }

      log_forward_probs(i, 0) = logSumExp_cpp(alpha(time-1, _));
      log_forward_probs(i, 1) = logSumExp_cpp(beta(time-1, _));
    }

    for(int i = 0; i < ndept; ++i) {
      rowlogsumexp[i] = logSumExp_cpp(log_forward_probs(i, _));
    }

    double full_log_likelihood = sum(rowlogsumexp);
    return full_log_likelihood;
  }

  // Model 3 or 6
  if(model == 3 || model == 6) {
    NumericVector init_density = state_dist_cpp(Gamma(0, 1), Gamma(1, 0));
    init_density = log(init_density);
    NumericMatrix Alphas(time, nstate);
    NumericMatrix alpha(time, nstate);
    NumericMatrix beta(time, nstate);
    NumericMatrix log_forward_probs(ndept, nstate);
    NumericVector rowlogsumexp(ndept);

    for(int i = 0; i < ndept; ++i) {
      if(y(i, 0) == -1){
        alpha(0, 0) = init_density[0];
        alpha(0, 1) = 0;
        beta(0, 0) = 0;
        beta(0, 1) = init_density[1];
      }else{
        alpha(0, 0) = init_density[0] + R::dpois(y(i, 0), e_it(i, 0) * exp(r[0] + s[0] + u[i]), true);
        alpha(0, 1) = 0;
        beta(0, 0) = 0;
        beta(0, 1) = init_density[1] + R::dpois(y(i, 0), e_it(i, 0) * exp(r[0] + s[0] + u[i] + B[0] * z_it(i, 0) + B[1] * z_it2(i, 0)), true);
      }
      Alphas(0, 0) = alpha(0, 0);
      Alphas(0, 1) = beta(0, 1);

      for(int t = 1; t < time; ++t) {
        if(y(i, t) == -1){
          alpha(t, 0) = Alphas(t-1, 0) + gamma_11;
          alpha(t, 1) = Alphas(t-1, 1) + gamma_21;
          beta(t, 0) = Alphas(t-1, 0) + gamma_12;
          beta(t, 1) = Alphas(t-1, 1) + gamma_22;
          Alphas(t, 0) = logSumExp_cpp(alpha(t, _));
          Alphas(t, 1) = logSumExp_cpp(beta(t, _));
        }else{
          int month_index = t % 12;
          alpha(t, 0) = Alphas(t-1, 0) + gamma_11 + R::dpois(y(i, t), e_it(i, t) * exp(r[t] + s[month_index] + u[i]), true);
          alpha(t, 1) = Alphas(t-1, 1) + gamma_21 + R::dpois(y(i, t), e_it(i, t) * exp(r[t] + s[month_index] + u[i]), true);
          beta(t, 0) = Alphas(t-1, 0) + gamma_12 + R::dpois(y(i, t), e_it(i, t) * exp(r[t] + s[month_index] + u[i] + B[0] * z_it(i, t) + B[1] * z_it2(i, t)), true);
          beta(t, 1) = Alphas(t-1, 1) + gamma_22 + R::dpois(y(i, t), e_it(i, t) * exp(r[t] + s[month_index] + u[i] + B[0] * z_it(i, t) + B[1] * z_it2(i, t)), true);
          Alphas(t, 0) = logSumExp_cpp(alpha(t, _));
          Alphas(t, 1) = logSumExp_cpp(beta(t, _));
        }
      }

      log_forward_probs(i, 0) = logSumExp_cpp(alpha(time-1, _));
      log_forward_probs(i, 1) = logSumExp_cpp(beta(time-1, _));
    }

    for(int i = 0; i < ndept; ++i) {
      rowlogsumexp[i] = logSumExp_cpp(log_forward_probs(i, _));
    }

    double full_log_likelihood = sum(rowlogsumexp);
    return full_log_likelihood;
  }

  return R_NegInf;
}

// Computing multi-type log likelihood using cyclic seasonal components
// [[Rcpp::export]]
double multGeneralLoglikelihood_cpp2(IntegerVector y, int ndept, int time, int nstrain, NumericVector a_k,
                                     NumericVector r, NumericVector s, NumericVector u, NumericMatrix Gamma,
                                     NumericMatrix e_it, NumericVector B, int model, NumericMatrix Bits,
                                     int independentChains){
 // Rcpp::Rcout << "Model selected: " << model << std::endl;
  // Model 0
  if(model == 0) {
    double full_log_likelihood = 0;
    for(int i = 0; i < ndept; ++i) {
      for(int t = 0; t < time; ++t) {
        int month_index = t % 12;
        for(int k = 0; k < nstrain; ++k){
          if(y[k * ndept * time + i * time + t] == -1){
            full_log_likelihood += 0;
          }else{
            full_log_likelihood += dpois_cpp(y[k * ndept * time + i * time + t], e_it(i, t) * exp(a_k[k] + r[t] + s[month_index] + u[i]));
//            full_log_likelihood += R::dpois(y[k * ndept * time + i * time + t], e_it(i, t) * exp(a_k[k] + r[t] + s[month_index] + u[i]), true);
          }
        }
      }
    }
    return full_log_likelihood;
  }
  // Model 1
  if(model == 1){
    if(independentChains == 0){
    int nstate = intPower(2, nstrain);
    NumericMatrix jointTPM = JointTransitionMatrix_cpp(Gamma, nstrain);
    NumericVector init_density = state_dist_cpp2(jointTPM);
    NumericVector rowlogsumexp(ndept);
    for(int i = 0; i < ndept; ++i) {
      NumericMatrix Alpha(time, nstate);
      NumericMatrix prodEmission(time, nstate);
      for(int n = 0; n < nstate; ++n){
        for(int k = 0; k < nstrain; ++k){
          if(y[k * ndept * time + i * time + 0] == -1){
            prodEmission(0, n) += 0;
          }else{
            prodEmission(0, n) = prodEmission(0, n) + dpois_cpp(y[k * ndept * time + i * time + 0], e_it(i, 0) * exp(a_k[k] + r[0] + s[0] + u[i] + dotproduct_cpp(B, Bits(n, _))));
//            prodEmission(0, n) = prodEmission(0, n) + R::dpois(y[k * ndept * time + i * time + 0], e_it(i, 0) * exp(a_k[k] + r[0] + s[0] + u[i] + dotproduct_cpp(B, Bits(n, _))), true);
          }
        }
      }
      for (int n = 0; n < nstate; ++n){
        Alpha(0, n) = log(init_density[n]) + prodEmission(0, n);
      }
      for(int t = 1; t < time; ++t){
        int month_index = t % 12;
        for(int n = 0; n < nstate; ++n){
          for(int k = 0; k < nstrain; ++k){
            if(y[k * ndept * time + i * time + t] == -1){
              prodEmission(t, n) += 0;
            }else{
              prodEmission(t, n) = prodEmission(t, n) + dpois_cpp(y[k * ndept * time + i * time + t], e_it(i, t) * exp(a_k[k] + r[t] + s[month_index] + u[i] + dotproduct_cpp(B, Bits(n, _))));
//              prodEmission(t, n) = prodEmission(t, n) + R::dpois(y[k * ndept * time + i * time + t], e_it(i, t) * exp(a_k[k] + r[t] + s[month_index] + u[i] + dotproduct_cpp(B, Bits(n, _))), true);
            }
          }
        }
        NumericVector temp = logVecMatMult(Alpha(t-1, _), jointTPM);
        for (int j = 0; j < nstate; ++j){
          Alpha(t, j) = temp[j] + prodEmission(t, j);
        }
      }
          rowlogsumexp[i] = logSumExp_cpp(Alpha(time-1, _));
    }
    double full_log_likelihood = sum(rowlogsumexp);
    return full_log_likelihood;
    }else{
      //Independent chains
      int nstate = Gamma.ncol();
      NumericVector init_density = state_dist_cpp2(Gamma);
      NumericVector rowlogsumexp(ndept);
      for(int i = 0; i < ndept; ++i) {
        NumericMatrix Alpha(time, nstate);
        NumericMatrix prodEmission(time, nstate);
        for(int n = 0; n < nstate; ++n){
          for(int k = 0; k < nstrain; ++k){
            if(y[k * ndept * time + i * time + 0] == -1){
              prodEmission(0, n) += 0;
            }else{
              prodEmission(0, n) = prodEmission(0, n) + dpois_cpp(y[k * ndept * time + i * time + 0], e_it(i, 0) * exp(a_k[k] + r[0] + s[0] + u[i] + B[k] * n));
//              prodEmission(0, n) = prodEmission(0, n) + R::dpois(y[k * ndept * time + i * time + 0], e_it(i, 0) * exp(a_k[k] + r[0] + s[0] + u[i] + B[k] * n), true);
            }
          }
        }
        for (int n = 0; n < nstate; ++n){
          Alpha(0, n) = log(init_density[n]) + prodEmission(0, n);
        }
        for(int t = 1; t < time; ++t){
          int month_index = t % 12;
          for(int n = 0; n < nstate; ++n){
            for(int k = 0; k < nstrain; ++k){
              if(y[k * ndept * time + i * time + t] == -1){
                prodEmission(t, n) += 0;
              }else{
                prodEmission(t, n) = prodEmission(t, n) + dpois_cpp(y[k * ndept * time + i * time + t], e_it(i, t) * exp(a_k[k] + r[t] + s[month_index] + u[i] + B[k] * n));
//               prodEmission(t, n) = prodEmission(t, n) + R::dpois(y[k * ndept * time + i * time + t], e_it(i, t) * exp(a_k[k] + r[t] + s[month_index] + u[i] + B[k] * n), true);
              }
            }
          }
          NumericVector temp = logVecMatMult(Alpha(t-1, _), Gamma);
          for (int j = 0; j < nstate; ++j){
            Alpha(t, j) = temp[j] + prodEmission(t, j);
          }
        }
        rowlogsumexp[i] = logSumExp_cpp(Alpha(time-1, _));
      }
      double full_log_likelihood = sum(rowlogsumexp);
      return full_log_likelihood;
    }
  }
  return R_NegInf;
}

// [[Rcpp::export]]
NumericVector hardconstraint_cpp(NumericVector x){
  int n = x.size();
  double sumX = 0;
  NumericVector x_new(n);
  for(int i = 0; i < n ; ++i) {
    x_new[i] = x[i];
    sumX += x[i];
  }
  x_new.push_back(-sumX);
  return x_new;
}

// [[Rcpp::export]]
NumericVector rcpp_rmvn(int n, NumericVector v, NumericMatrix m){
  // Obtaining namespace of mvnfast package
  Environment pkg = Environment::namespace_env("mvnfast");
  // Picking up rmvn function from mvnfast package
  Function f = pkg["rmvn"];
  // Executing Matrix( m, sparse = TRIE )
  return f(n, v, m);
}

// [[Rcpp::export]]
double rcpp_dmvn(NumericVector x, NumericVector u, NumericMatrix m){
  // Obtaining namespace of mvnfast package
  Environment pkg = Environment::namespace_env("mvnfast");
  // Picking up dmvn function from mvnfast package
  Function f = pkg["dmvn"];
  // Executing dmvn
  return as<double>( f(_["X"] = x, _["mu"] = u, _["sigma"] = m, _["log"] = true));
}

// [[Rcpp::export]]
NumericMatrix rcpp_updateCov(NumericMatrix C, NumericVector x, int n, NumericVector xbar){
  // Obtaining namespace of mvnfast package
  Environment pkg = Environment::namespace_env("onlinePCA");
  // Picking up updateCovariance function from onlinePCA package
  Function f = pkg["updateCovariance"];
  // Executing updateCovariance
  return f(C, x, n, xbar);
}

// [[Rcpp::export]]
NumericMatrix makematrix_cpp(double g12, double g21){
  NumericMatrix Gmat(2, 2);
  Gmat(0, 0) = 1 - g12;
  Gmat(0, 1) = g12;
  Gmat(1, 0) = g21;
  Gmat(1, 1) = 1 - g21;
  return Gmat;
}

// [[Rcpp::export]]
NumericMatrix tcrossprod_cpp(NumericVector x) {
  int n = x.size();
  NumericMatrix result(n, n);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      result(i, j) = x[i] * x[j];
    }
  }
  return result;
}

// [[Rcpp::export]]
NumericVector ColMeans_cpp(NumericMatrix M) {
  int nRows  = M.nrow();
  int nCols  = M.ncol();
  NumericVector means(nCols);
  for (int c = 0; c < nCols; ++c) {
    double sum = 0.0;
    for (int r = 0; r < nRows; ++r) {
      sum += M(r, c);
    }
    means[c] = sum / nRows;
  }
  return means;
}

// [[Rcpp::export]]
NumericMatrix makeDiagMat_cpp(double epsilonS, int n) {
  NumericMatrix mat(n, n);
  for (int i = 0; i < n; i++) {
    mat(i, i) = epsilonS;
  }
  return mat;
}

// [[Rcpp::export]]
NumericMatrix cov_cpp(NumericMatrix X){
  int n = X.nrow();
  int k = X.ncol();

  NumericVector means = ColMeans_cpp(X);
  NumericMatrix centered(n, k);
  for (int j = 0; j < k; j++) {
    for (int i = 0; i < n; i++) {
      centered(i, j) = X(i, j) - means[j];
    }
  }
  NumericMatrix covmat(k, k);
  for (int j1 = 0; j1 < k; j1++) {
    for (int j2 = j1; j2 < k; j2++) {
      double sum = 0.0;
      for (int i = 0; i < n; i++) {
        sum += centered(i, j1) * centered(i, j2);
      }
      double cov_val = sum / (n - 1);
      covmat(j1, j2) = cov_val;
      covmat(j2, j1) = cov_val; //by symmetry
    }
  }
  return covmat;
}

// [[Rcpp::export]]
NumericMatrix add2matrices_cpp(NumericMatrix A, NumericMatrix B) {
  int nr = A.nrow(), nc = A.ncol();
  if (B.nrow() != nr || B.ncol() != nc) {
    stop("Matrices must have the same dimensions");
  }
  NumericMatrix C(nr, nc);
  for (int j = 0; j < nc; j++) {
    for (int i = 0; i < nr; i++) {
      C(i, j) = A(i, j) + B(i, j);
    }
  }
  return C;
}

// [[Rcpp::export]]
NumericMatrix subtract2matrices_cpp(NumericMatrix A, NumericMatrix B) {
  int nr = A.nrow(), nc = A.ncol();
  if (B.nrow() != nr || B.ncol() != nc) {
    stop("Matrices must have the same dimensions");
  }
  NumericMatrix C(nr, nc);
  for (int j = 0; j < nc; j++) {
    for (int i = 0; i < nr; i++) {
      C(i, j) = A(i, j) - B(i, j);
    }
  }
  return C;
}

// [[Rcpp::export]]
NumericMatrix multiply2matrices_cpp(NumericMatrix A, NumericMatrix B) {
  int nrowA = A.nrow();
  int ncolA = A.ncol();
  int nrowB = B.nrow();
  int ncolB = B.ncol();

  if (ncolA != nrowB) {
    stop("Incompatible matrix dimensions for multiplication.");
  }

  NumericMatrix C(nrowA, ncolB);
  for (int i = 0; i < nrowA; i++) {
    for (int j = 0; j < ncolB; j++) {
      double sum = 0.0;
      for (int k = 0; k < ncolA; k++) {
        sum += A(i, k) * B(k, j);
      }
      C(i, j) = sum;
    }
  }
  return C;
}

// [[Rcpp::export]]
NumericMatrix inv_cpp(NumericMatrix x) {
  int n = x.nrow();
  if (n != x.ncol()) stop("Matrix must be square.");

  NumericMatrix A(clone(x));       // make a copy so we don't overwrite input
  NumericMatrix I(n, n);

  // Initialize I to identity matrix
  for (int i = 0; i < n; i++) {
    I(i, i) = 1.0;
  }

  // Forward elimination
  for (int i = 0; i < n; i++) {
    double pivot = A(i, i);
    if (pivot == 0) stop("Matrix is singular.");

    // Normalize pivot row
    for (int j = 0; j < n; j++) {
      A(i, j) /= pivot;
      I(i, j) /= pivot;
    }

    // Eliminate other rows
    for (int k = 0; k < n; k++) {
      if (k != i) {
        double factor = A(k, i);
        for (int j = 0; j < n; j++) {
          A(k, j) -= factor * A(i, j);
          I(k, j) -= factor * I(i, j);
        }
      }
    }
  }

  return I;
}

//MCMC sampling
// [[Rcpp::export]]
NumericMatrix multInfer_cpp(IntegerVector y, NumericMatrix e_it, int nstrain, int Model, NumericMatrix Bits, NumericVector CrudeR,
                            NumericVector CrudeS, NumericVector CrudeU, NumericMatrix RW2PrecMat, NumericMatrix RW1PrecMat,
                            NumericMatrix R, int rankdef, int independentChains, int num_iteration, double meanR,
                            NumericVector SumYk_vec){

  int time = e_it.ncol();
  int ndept = e_it.nrow();
  double optconstantR = std::pow(2.38, 2)/(time-2);
  double optconstantS = std::pow(2.38, 2)/11;
  double optconstantU = std::pow(2.38, 2)/(ndept-1);
  double lambdaR = 1.0;
  double lambdaS = 1.0;
  double lambdaU = 1.0;
  NumericMatrix MC_chain(num_iteration, 5+time+12+ndept+nstrain+nstrain+1);
  //initials
  MC_chain(0, 0) = R::runif(0,1);
  MC_chain(0, 1) = R::runif(0,1);
  MC_chain(0, 2) = R::runif(0,1000);
  MC_chain(0, 3) = R::runif(0,100);
  MC_chain(0, 4) = R::runif(0,30);
  MC_chain(0, 5+time+12+ndept+nstrain+nstrain) = state_dist_cpp(MC_chain(0,0), MC_chain(0,1))[1];
  for(int t = 0; t < time; ++t){
    MC_chain(0, 5+t) = CrudeR[t];
  }
  for(int s = 0; s < 12; ++s){
    MC_chain(0, 5+time+s) = CrudeS[s];
  }
  for(int u = 0; u < ndept; ++u){
    MC_chain(0, 5+time+12+u) = CrudeU[u];
  }
  for(int b = 0; b < nstrain ; ++b){
    MC_chain(0, 5+time+12+ndept+b) = 0.0;
  }
  for(int a = 0; a < nstrain ; ++a){
    MC_chain(0, 5+time+12+ndept+nstrain+a) = meanR;
  }

  double epsilonR = 0.007;
  NumericMatrix XnR(5, time);
  NumericVector XnbarR(time);
  NumericVector XnbarPrevR(time);
  NumericMatrix currentzigmaR(time, time);
  NumericMatrix tempzigmaR(time, time);

  currentzigmaR = makeDiagMat_cpp(0.1,time);


  double epsilonS = 0.007;
  NumericMatrix XnS(5, 11);
  NumericVector XnbarS(11);
  NumericVector XnbarPrevS(11);
  NumericMatrix currentzigmaS(11, 11);
  NumericMatrix tempzigmaS(11, 11);

  currentzigmaS = makeDiagMat_cpp(0.1,11);

  double epsilonU = 0.007;
  NumericMatrix XnU(5, ndept-1);
  NumericVector XnbarU(ndept-1);
  NumericVector XnbarPrevU(ndept-1);
  NumericMatrix currentzigmaU(ndept-1, ndept-1);
  NumericMatrix tempzigmaU(ndept-1, ndept-1);

  currentzigmaU = makeDiagMat_cpp(0.08,ndept-1);

  int Rsize = time/8;

  NumericMatrix invRconditionalcovA = inv_cpp(RW2PrecMat(Range(0,Rsize-1), Range(0,Rsize-1)));
  NumericMatrix invRconditionalcovB = inv_cpp(RW2PrecMat(Range(Rsize, 2*Rsize-1), Range(Rsize, 2*Rsize-1)));
  NumericMatrix invRconditionalcovC = inv_cpp(RW2PrecMat(Range(2*Rsize,3*Rsize-1), Range(2*Rsize, 3*Rsize-1)));
  NumericMatrix invRconditionalcovD = inv_cpp(RW2PrecMat(Range(3*Rsize,4*Rsize-1), Range(3*Rsize,4*Rsize-1)));
  NumericMatrix invRconditionalcovE = inv_cpp(RW2PrecMat(Range(4*Rsize,5*Rsize-1), Range(4*Rsize,5*Rsize-1)));
  NumericMatrix invRconditionalcovf = inv_cpp(RW2PrecMat(Range(5*Rsize,6*Rsize-1), Range(5*Rsize,6*Rsize-1)));
  NumericMatrix invRconditionalcovg = inv_cpp(RW2PrecMat(Range(6*Rsize,7*Rsize-1), Range(6*Rsize,7*Rsize-1)));
  NumericMatrix invRconditionalcovH = inv_cpp(RW2PrecMat(Range(7*Rsize,time-1), Range(7*Rsize,time-1)));

  NumericMatrix currentRW2PrecMat(time, time);

  NumericVector acceptedS(12);
  NumericVector acceptedU(ndept);
  NumericVector acceptedR(time);
  NumericVector acceptedB(nstrain);
  NumericVector accepteda_k(nstrain);

  //Start MCMC sampling
  for(int i = 1; i < num_iteration; ++i){

    //Kappa_r update
    NumericMatrix matR = MC_chain(Range(i-1,i-1), Range(5, 5+(time-1)));
    NumericVector currentR = matR(0, _);
    MC_chain(i, 2) = R::rgamma(1 + (time-2)/2.0, 1.0/(0.0001 + quadform_cpp(currentR, RW2PrecMat)/2.0));

    //Kappa_s update
    NumericMatrix matS = MC_chain(Range(i-1, i-1), Range(5+time, 5+time+11));
    NumericVector currentS = matS(0, _);
    MC_chain(i, 3) = R::rgamma(1 + 11/2.0, 1.0/(0.001 + quadform_cpp(currentS, RW1PrecMat)/2.0));

    //kappa_u update
    NumericMatrix matU = MC_chain(Range(i-1,i-1), Range(5+time+12, 5+time+12+(ndept-1)));
    NumericVector currentU = matU(0, _);
    MC_chain(i, 4) = R::rgamma(1 + (ndept-rankdef)/2.0, 1.0/(0.01 + quadform_cpp(currentU, R)/2.0));

    //store current values Gamma matrix, Betas and a_k's
    NumericMatrix Gmat = makematrix_cpp(MC_chain(i-1, 0), MC_chain(i-1, 1));

    NumericMatrix matB = MC_chain(Range(i-1,i-1), Range(5+time+12+ndept, 5+time+12+ndept+(nstrain-1)));
    NumericVector currentB = matB(0, _);

    NumericMatrix mata_k = MC_chain(Range(i-1,i-1), Range(5+time+12+ndept+nstrain, 5+time+12+ndept+nstrain+(nstrain-1)));
    NumericVector currenta_k = mata_k(0, _);

    //spatial components update (U)
    NumericVector tempU = currentU;
    tempU.erase(tempU.begin()+ndept-1);
    NumericVector proposedUcomps = rcpp_rmvn(1,tempU,currentzigmaU);
    NumericVector newproposedUcomps = hardconstraint_cpp(proposedUcomps);

    double priorcurrentUcomps = logIGMRF1_cpp(currentU, MC_chain(i, 4), R, rankdef);
    double priorproposedUcomps = logIGMRF1_cpp(newproposedUcomps, MC_chain(i, 4), R, rankdef);

    double proposalproposedcompsU = rcpp_dmvn(proposedUcomps, tempU, currentzigmaU);
    double proposalcurrentcompsU = rcpp_dmvn(tempU, proposedUcomps, currentzigmaU);

    double likelihoodcurrentU = multGeneralLoglikelihood_cpp2(y,ndept,time,nstrain,currenta_k, currentR,currentS,currentU,Gmat,e_it, currentB, Model, Bits, independentChains);
    double likelihoodproposedU = multGeneralLoglikelihood_cpp2(y,ndept,time,nstrain,currenta_k,currentR,currentS,newproposedUcomps,Gmat,e_it, currentB, Model, Bits, independentChains);

    double mh_ratioU = exp(likelihoodproposedU + priorproposedUcomps + proposalcurrentcompsU - likelihoodcurrentU - priorcurrentUcomps - proposalproposedcompsU);

    if(R::runif(0, 1) < mh_ratioU) {
      for (int u = 0; u < ndept; ++u) {
        MC_chain(i, 5 + time + 12 + u) = newproposedUcomps[u];
        acceptedU[u] = newproposedUcomps[u];
      }
    }else{
      for (int u = 0; u < ndept; ++u) {
        MC_chain(i, 5 + time + 12 + u) = currentU[u];
        acceptedU[u] = currentU[u];
      }
    }

    //update trend components (R)
    double tempKappaR = 1.0 / MC_chain(i, 2);
    currentRW2PrecMat = MC_chain(i, 2) * RW2PrecMat;
    NumericMatrix RconditionalcovA = tempKappaR * invRconditionalcovA;
    NumericMatrix RconditionalcovB = tempKappaR * invRconditionalcovB;
    NumericMatrix RconditionalcovC = tempKappaR * invRconditionalcovC;
    NumericMatrix RconditionalcovD = tempKappaR * invRconditionalcovD;
    NumericMatrix RconditionalcovE = tempKappaR * invRconditionalcovE;
    NumericMatrix Rconditionalcovf = tempKappaR * invRconditionalcovf;
    NumericMatrix Rconditionalcovg = tempKappaR * invRconditionalcovg;
    NumericMatrix RconditionalcovH = tempKappaR * invRconditionalcovH;

    List covBlocks = List::create(RconditionalcovA, RconditionalcovB, RconditionalcovC, RconditionalcovD,
                                  RconditionalcovE, Rconditionalcovf, Rconditionalcovg, RconditionalcovH);

      for(int j = 0; j < 8 ; ++j){
        NumericVector proposedRcomps(time);
        NumericVector Rconditionalmean(Rsize);
        NumericVector blockproposedRcomps(Rsize);
        NumericVector lastproposedRcomps(time-7*Rsize);
        NumericVector lastproposedRconditionalmean(time-7*Rsize);

        if(j==0){
          NumericMatrix othercompsMat = MC_chain(Range(i-1, i-1), Range(5+Rsize, 5+time-1));
          NumericMatrix tempCov = covBlocks[j];
          NumericMatrix newtempCov = -1 * tempCov;
          NumericMatrix tempRW2Mat = currentRW2PrecMat(Range(0,Rsize-1), Range(Rsize, time-1));
          NumericMatrix RconditionalmeanMat = multiply2matrices_cpp(multiply2matrices_cpp(newtempCov, tempRW2Mat), transpose(othercompsMat));
          Rconditionalmean = RconditionalmeanMat(_, 0);
          blockproposedRcomps = rcpp_rmvn(1, Rconditionalmean, covBlocks[j]);
          //Rcpp::Rcout << "Succesful" << i << std::endl;
          for(int t = 0; t < Rsize; ++t){
          proposedRcomps[t] = blockproposedRcomps[t];
          }
          for(int tt = Rsize; tt < time; ++tt){
            proposedRcomps[tt] = MC_chain(i-1, 5+tt);
          }
        }
        else if(j!=0 && j!=7){
          NumericMatrix Before_othercompsMat = MC_chain(Range(i, i), Range(5, 5+j*Rsize-1));
          NumericMatrix After_othercompsMat = MC_chain(Range(i, i), Range(5+(j+1)*Rsize, 5+time-1));
          NumericMatrix tempCov = covBlocks[j];
          NumericMatrix newtempCov = -1 * tempCov;
          NumericMatrix firsttempRW2Mat = currentRW2PrecMat(Range(j*Rsize, (j+1)*Rsize-1), Range(0, j*Rsize-1));
          NumericMatrix secondtempRW2Mat = currentRW2PrecMat(Range(j*Rsize, (j+1)*Rsize-1), Range((j+1)*Rsize, time-1));
          NumericMatrix RconditionalmeanMat1 = multiply2matrices_cpp(firsttempRW2Mat, transpose(Before_othercompsMat));
          NumericMatrix RconditionalmeanMat2 = multiply2matrices_cpp(secondtempRW2Mat, transpose(After_othercompsMat));
          //Rcpp::Rcout << "Succesful" << j << std::endl;
          NumericMatrix RconditionalmeanMat = multiply2matrices_cpp(newtempCov, add2matrices_cpp(RconditionalmeanMat1, RconditionalmeanMat2));

         // RconditionalmeanMat = add2matrices_cpp(multiply2matrices_cpp(multiply2matrices_cpp(newtempCov, firsttempRW2Mat), transpose(Before_othercompsMat)), multiply2matrices_cpp(secondtempRW2Mat, transpose(After_othercompsMat)));
          Rconditionalmean = RconditionalmeanMat(_, 0);
          blockproposedRcomps = rcpp_rmvn(1, Rconditionalmean, covBlocks[j]);
          for(int t = 0; t < j*Rsize; ++t){
          proposedRcomps[t] = MC_chain(i, 5+t);
          }
          for(int tt = j*Rsize; tt < (j+1)*Rsize; ++tt){
            proposedRcomps[tt] = blockproposedRcomps[tt-j*Rsize];
          }
          for(int ttt = (j+1)*Rsize; ttt < time; ++ttt){
            proposedRcomps[ttt] = MC_chain(i, 5+ttt);
          }
        }
        else if(j==7){
          NumericMatrix lastproposedRconditionalmeanMat(time-7*Rsize, 1);
          NumericMatrix othercompsMat = MC_chain(Range(i, i), Range(5, 5+(7*Rsize)-1));
          NumericMatrix tempCov = covBlocks[j];
          NumericMatrix newtempCov = -1 * tempCov;
          NumericMatrix tempRW2Mat = currentRW2PrecMat(Range(7*Rsize, time-1), Range(0, 7*Rsize-1));
          lastproposedRconditionalmeanMat = multiply2matrices_cpp(multiply2matrices_cpp(newtempCov, tempRW2Mat), transpose(othercompsMat));
          //Rcpp::Rcout << "Succesful" << j << std::endl;
          lastproposedRconditionalmean = lastproposedRconditionalmeanMat(_, 0);
          blockproposedRcomps = rcpp_rmvn(1, lastproposedRconditionalmean, covBlocks[j]);
          for(int t = 0; t < 7*Rsize; ++t){
            proposedRcomps[t] = MC_chain(i, 5+t);
          }
          for(int tt = 7*Rsize; tt < time; ++tt){
            proposedRcomps[tt] = blockproposedRcomps[tt-7*Rsize];
          }
        }

        double likelihoodcurrentR = multGeneralLoglikelihood_cpp2(y,ndept,time,nstrain,currenta_k, currentR,currentS,acceptedU,Gmat,e_it, currentB, Model, Bits, independentChains);
        double likelihoodproposedR = multGeneralLoglikelihood_cpp2(y,ndept,time,nstrain,currenta_k,proposedRcomps,currentS,acceptedU,Gmat,e_it, currentB, Model, Bits, independentChains);

        double mh_ratioR = exp(likelihoodproposedR - likelihoodcurrentR);

          if(R::runif(0,1) < mh_ratioR){
            for(int t = 0; t < time; ++t){
            MC_chain(i, 5+t) = proposedRcomps[t];
            acceptedR[t] = proposedRcomps[t];
            }
          }
          else{
            if(j==0){
              for(int t = 0; t < time; ++t){
              MC_chain(i, 5+t) = MC_chain(i-1, 5+t);
              acceptedR[t] = MC_chain(i-1, 5+t);
              }
            }
            else if(j!=0){
              for(int t = 0; t < time; ++t){
                MC_chain(i, 5+t) = MC_chain(i, 5+t);
                acceptedR[t] = MC_chain(i, 5+t);
              }
            }
          }
      }
      //Rcpp::Rcout << "Condprior" << i << std::endl;

    //seasonal components update (S)
    NumericVector tempS = currentS;
    tempS.erase(tempS.begin()+11);
    NumericVector proposedScomps = rcpp_rmvn(1,tempS,currentzigmaS);
    NumericVector newproposedScomps = hardconstraint_cpp(proposedScomps);

    double priorcurrentScomps = seasonalComp2_cpp(currentS, MC_chain(i, 3), RW1PrecMat);
    double priorproposedScomps = seasonalComp2_cpp(newproposedScomps, MC_chain(i, 3), RW1PrecMat);

    double proposalproposedcompsS = rcpp_dmvn(proposedScomps, tempS, currentzigmaS);
    double proposalcurrentcompsS = rcpp_dmvn(tempS, proposedScomps, currentzigmaS);

    double likelihoodcurrentS = multGeneralLoglikelihood_cpp2(y,ndept,time,nstrain,currenta_k, acceptedR,currentS,acceptedU,Gmat,e_it, currentB, Model, Bits, independentChains);
    double likelihoodproposedS = multGeneralLoglikelihood_cpp2(y,ndept,time,nstrain,currenta_k,acceptedR,newproposedScomps,acceptedU,Gmat,e_it, currentB, Model, Bits, independentChains);

    double mh_ratioS = exp(likelihoodproposedS + priorproposedScomps - likelihoodcurrentS - priorcurrentScomps);

    if(R::runif(0, 1) < mh_ratioS) {
      for (int s = 0; s < 12; ++s) {
        MC_chain(i, 5 + time + s) = newproposedScomps[s];
        acceptedS[s] = newproposedScomps[s];
      }
    }else{
      for (int s = 0; s < 12; ++s) {
        MC_chain(i, 5 + time + s) = currentS[s];
        acceptedS[s] = currentS[s];
      }
    }
    //Rcpp::Rcout << "S update" << i << std::endl;

    // update Betas
    NumericVector proposedB(nstrain);
    double priorcurrentB = 0.0;
    double priorproposedB = 0.0;

    if(Model == 0){
      proposedB = rep(0.0, nstrain);
      priorcurrentB = -99999.0;
      priorproposedB = -99999.0;
    }else{
    for (int j = 0; j < nstrain; ++j) {
      proposedB[j] = std::abs(R::rnorm(currentB[j], 0.1));
      priorcurrentB += R::dgamma(currentB[j], 2, 1.0 / 2.0, true);
      priorproposedB += R::dgamma(proposedB[j], 2, 1.0 / 2.0, true);
      }
    }

    double likelihoodcurrentB = multGeneralLoglikelihood_cpp2(y, ndept, time, nstrain, currenta_k, acceptedR, acceptedS, acceptedU,Gmat, e_it, currentB, Model, Bits, independentChains);
    double likelihoodproposedB = multGeneralLoglikelihood_cpp2(y, ndept, time, nstrain, currenta_k, acceptedR, acceptedS, acceptedU,Gmat, e_it, proposedB, Model, Bits, independentChains);

    double mh_ratioB = exp(likelihoodproposedB + priorproposedB
                             - likelihoodcurrentB - priorcurrentB);

    if(R::runif(0, 1) < mh_ratioB) {
      for (int b = 0; b < nstrain; ++b) {
        MC_chain(i, 5 + time + 12 + ndept + b) = proposedB[b];
        acceptedB[b] = proposedB[b];
      }
    } else {
      for (int b = 0; b < nstrain; ++b) {
        MC_chain(i, 5 + time + 12 + ndept + b) = currentB[b];
        acceptedB[b] = currentB[b];
      }
    }
    //Rcpp::Rcout << "Beta update" << i << std::endl;

    // update a_k's
   // NumericVector proposeda_k(nstrain);
    //  for (int j = 0; j < nstrain; ++j) {
    //    proposeda_k[j] = R::rnorm(currenta_k[j], 0.07);
    //  }

  //  double likelihoodcurrenta_k = multGeneralLoglikelihood_cpp2(y, ndept, time, nstrain, currenta_k, acceptedR, acceptedS, acceptedU,Gmat, e_it, acceptedB, Model, Bits, independentChains);
  //  double likelihoodproposeda_k = multGeneralLoglikelihood_cpp2(y, ndept, time, nstrain, proposeda_k, acceptedR, acceptedS, acceptedU,Gmat, e_it, acceptedB, Model, Bits, independentChains);

//    double mh_ratioA = exp(likelihoodproposeda_k
//                             - likelihoodcurrenta_k);

  //  NumericVector accepteda_k(nstrain);
  //  if(R::runif(0, 1) < mh_ratioA) {
  //    for (int a = 0; a < nstrain; ++a) {
  //      MC_chain(i, 5 + time + 12 + ndept + nstrain + a) = proposeda_k[a];
  //      accepteda_k[a] = proposeda_k[a];
  //    }
  //  } else {
  //    for (int a = 0; a < nstrain; ++a) {
  //      MC_chain(i, 5 + time + 12 + ndept + nstrain + a) = currenta_k[a];
  //      accepteda_k[a] = currenta_k[a];
  //    }
  //  }
    //Rcpp::Rcout << "a_k update" << i << std::endl;

    //update Gamma's
    NumericVector tempGs(2);
    NumericVector proposedGs(2);
    NumericVector currentGs(2);
    double priorcurrentGs = 0.0;
    double priorproposedGs = 0.0;
    for(int g = 0; g < 2; ++g){
    currentGs[g] = MC_chain(i-1, g);
    tempGs[g] = std::abs(R::rnorm(currentGs[g], 0.1));
    if(tempGs[0]>1){
      proposedGs[0]=2-tempGs[0];
    }else{
      proposedGs[0]=tempGs[0];
    }
    if(tempGs[1]>1){
      proposedGs[1]=2-tempGs[1];
    }else{
      proposedGs[1]=tempGs[1];
      }
    priorcurrentGs += R::dbeta(currentGs[g], 2.0, 2.0, true);
    priorproposedGs += R::dbeta(proposedGs[g], 2.0, 2.0, true);
    }
    NumericMatrix proposedGmat = makematrix_cpp(proposedGs[0], proposedGs[1]);

    double likelihoodcurrentGs = multGeneralLoglikelihood_cpp2(y, ndept, time, nstrain, currenta_k, acceptedR, acceptedS, acceptedU,Gmat, e_it, acceptedB, Model, Bits, independentChains);
    double likelihoodproposedGs = multGeneralLoglikelihood_cpp2(y, ndept, time, nstrain, currenta_k, acceptedR, acceptedS, acceptedU,proposedGmat, e_it, acceptedB, Model, Bits, independentChains);

    double mh_ratioG = exp(likelihoodproposedGs + priorproposedGs
                         - likelihoodcurrentGs - priorcurrentGs);

    if(R::runif(0,1) < mh_ratioG){
      MC_chain(i, 0) = proposedGs[0];
      MC_chain(i, 1) = proposedGs[1];
    }
    else{
      MC_chain(i, 0) = currentGs[0];
      MC_chain(i, 1) = currentGs[1];
    }
    MC_chain(i, 5+time+12+ndept+nstrain+nstrain) = state_dist_cpp(MC_chain(i, 0), MC_chain(i, 1))[1];
    //Rcpp::Rcout << "gamma update" << i << std::endl;

 //Gibbs update for a_k's
    NumericVector stationDist = state_dist_cpp(MC_chain(i, 0), MC_chain(i, 1));
      double poisMean = 0.0;
      for(int a = 0; a < ndept; ++a){
        for(int b = 0; b < time; ++b){
          int month_index = b % 12;
          poisMean +=  e_it(a, b) * exp(acceptedR[b] + acceptedS[month_index] + acceptedU[a] + dotproduct_cpp(acceptedB, stationDist));
        }
      }

      double scale = 1.0 / (poisMean + 0.01/exp(-15));
      for(int a = 0; a < nstrain; ++a){
      MC_chain(i, 5+time+12+ndept+nstrain+a) = log(R::rgamma(0.01+SumYk_vec[a],  scale));
      accepteda_k[a] = MC_chain(i, 5+time+12+ndept+nstrain+a);
      }

    //Adapting zigmaR
    if(i==5){
      XnR = MC_chain(Range(0,i-1), Range(5, 5+time-1));
      XnbarR = ColMeans_cpp(XnR);
      tempzigmaR = add2matrices_cpp(cov_cpp(XnR), makeDiagMat_cpp(epsilonR, time));
      currentzigmaR = optconstantR * tempzigmaR;
    } else if (i > 5){
      //trend components random-walk update (R)
      NumericVector RWproposedRcomps = rcpp_rmvn(1, acceptedR, currentzigmaR);

      double RWpriorcurrentRcomps = randomwalk2_cpp(acceptedR, MC_chain(i, 2));
      double RWpriorproposedRcomps = randomwalk2_cpp(RWproposedRcomps, MC_chain(i, 2));

      double RWproposalproposedcompsR = rcpp_dmvn(RWproposedRcomps, acceptedR, currentzigmaR);
      double RWproposalcurrentcompsR = rcpp_dmvn(acceptedR, RWproposedRcomps, currentzigmaR);

      double RWlikelihoodcurrentR = multGeneralLoglikelihood_cpp2(y,ndept,time,nstrain,accepteda_k, acceptedR,acceptedS,acceptedU,makematrix_cpp(MC_chain(i, 0), MC_chain(i, 1)),e_it, acceptedB, Model, Bits, independentChains);
      double RWlikelihoodproposedR = multGeneralLoglikelihood_cpp2(y,ndept,time,nstrain,accepteda_k,RWproposedRcomps,acceptedS,acceptedU,makematrix_cpp(MC_chain(i, 0), MC_chain(i, 1)),e_it, acceptedB, Model, Bits, independentChains);

      double RWmh_ratioR = exp(RWlikelihoodproposedR + RWpriorproposedRcomps + RWproposalcurrentcompsR - RWlikelihoodcurrentR - RWpriorcurrentRcomps - RWproposalproposedcompsR);

      NumericVector RWacceptedR(time);
      if(R::runif(0, 1) < RWmh_ratioR){
        for (int t = 0; t < time; ++t) {
          MC_chain(i, 5 + t) = RWproposedRcomps[t];
          RWacceptedR[t] = RWproposedRcomps[t];
        }
      }else{
        for (int t = 0; t < time; ++t) {
          MC_chain(i, 5 + t) = acceptedR[t];
          RWacceptedR[t] = acceptedR[t];
        }
      }
      XnbarPrevR = XnbarR;
      NumericVector tempXnR = i * XnbarR;
      XnbarR = (tempXnR + RWacceptedR)/(i+1);
      NumericMatrix tempzigmaR = rcpp_updateCov(currentzigmaR, RWacceptedR, i-1, XnbarPrevR);
      tempzigmaR = add2matrices_cpp(tempzigmaR, makeDiagMat_cpp(epsilonR, time));
      // Robbins–Monro tuning
      lambdaR = lambdaR * std::exp((2.0 / std::max(1, i - 5)) *
        (std::min(RWmh_ratioR, 1.0) - 0.234));
      currentzigmaR = lambdaR * optconstantR * tempzigmaR;
    }

//Adapting zigmaS
    if(i==5){
      XnS = MC_chain(Range(0,i-1), Range(5+time, 5+time+10));
      XnbarS = ColMeans_cpp(XnS);
      tempzigmaS = add2matrices_cpp(cov_cpp(XnS), makeDiagMat_cpp(epsilonS, 11));
      currentzigmaS = optconstantS * tempzigmaS;
    } else if (i > 5){
      XnbarPrevS = XnbarS;
      NumericVector tempXnS = i * XnbarS;
      NumericVector tempacceptedS = acceptedS;
      tempacceptedS.erase(tempacceptedS.begin()+11);
      XnbarS = (tempXnS + tempacceptedS)/(i+1);
      NumericMatrix tempzigmaS = rcpp_updateCov(currentzigmaS, tempacceptedS, i-1, XnbarPrevS);
      tempzigmaS = add2matrices_cpp(tempzigmaS, makeDiagMat_cpp(epsilonS, 11));
      // Robbins–Monro tuning
      lambdaS = lambdaS * std::exp((2.0 / std::max(1, i - 5)) *
        (std::min(mh_ratioS, 1.0) - 0.234));
      currentzigmaS = lambdaS * optconstantS * tempzigmaS;
    }

    //Adapting zigmaU
    if(i==5){
      XnU = MC_chain(Range(0,i-1), Range(5+time+12, 5+time+12+(ndept-2)));
      XnbarU = ColMeans_cpp(XnU);
      tempzigmaU = add2matrices_cpp(cov_cpp(XnU), makeDiagMat_cpp(epsilonU, ndept-1));
      currentzigmaU = optconstantU * tempzigmaU;
    } else if (i > 5){
      XnbarPrevU = XnbarU;
      NumericVector tempXnU = i * XnbarU;
      NumericVector tempacceptedU = acceptedU;
      tempacceptedU.erase(tempacceptedU.begin()+ndept-1);
      XnbarU = (tempXnU + tempacceptedU)/(i+1);
      NumericMatrix tempzigmaU = rcpp_updateCov(currentzigmaU, tempacceptedU, i-1, XnbarPrevU);
      tempzigmaU = add2matrices_cpp(tempzigmaU, makeDiagMat_cpp(epsilonU, ndept-1));
      // Robbins–Monro tuning
      lambdaU = lambdaU * std::exp((2.0 / std::max(1, i - 5)) *
        (std::min(mh_ratioU, 1.0) - 0.234));
      currentzigmaU = lambdaU * optconstantU * tempzigmaU;
    }
    if(i % 100 == 0) Rcpp::Rcout << "Iteration: " << i << std::endl;
  }
  return MC_chain;
}
