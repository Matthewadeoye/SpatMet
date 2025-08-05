#include <Rcpp.h>
#include <RcppEigen.h>
#include <cmath>
#include <algorithm>
using namespace Rcpp;

// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]
double logSumExp_cpp(NumericVector x) {
  double max_val = max(x);
  return log(sum(exp(x - max_val))) + max_val;
}

// [[Rcpp::export]]
double dot_product(NumericVector v1, NumericVector v2) {
  int lengthV1 = v1.size();
  double result = 0;
  for (int i = 0; i < lengthV1 ; ++i) {
    result += v1[i] * v2[i];
  }
  return result;
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
            full_log_likelihood += R::dpois(y[k * ndept * time + i * time + t], e_it(i, t) * exp(a_k[k] + r[t] + s[month_index] + u[i]), true);
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
            prodEmission(0, n) = prodEmission(0, n) + R::dpois(y[k * ndept * time + i * time + 0], e_it(i, 0) * exp(a_k[k] + r[0] + s[0] + u[i] + dot_product(B, Bits(n, _))), true);
         // Rcpp::Rcout << "DotProd" << dot_product(B, Bits(n, _)) << std::endl;
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
              prodEmission(t, n) = prodEmission(t, n) + R::dpois(y[k * ndept * time + i * time + t], e_it(i, t) * exp(a_k[k] + r[t] + s[month_index] + u[i] + dot_product(B, Bits(n, _))), true);
           // Rcpp::Rcout << "DotProd" << dot_product(B, Bits(n, _)) << std::endl;
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
              prodEmission(0, n) = prodEmission(0, n) + R::dpois(y[k * ndept * time + i * time + 0], e_it(i, 0) * exp(a_k[k] + r[0] + s[0] + u[i] + B[k] * n), true);
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
                prodEmission(t, n) = prodEmission(t, n) + R::dpois(y[k * ndept * time + i * time + t], e_it(i, t) * exp(a_k[k] + r[t] + s[month_index] + u[i] + B[k] * n), true);
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
