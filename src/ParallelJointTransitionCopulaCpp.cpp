#include <RcppArmadillo.h>
#include <RcppThread.h>
#include <omp.h>
#include <mvtnormAPI.h>
#include <cmath>
#include <algorithm>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppThread)]]
// [[Rcpp::depends(mvtnorm)]]
// [[Rcpp::plugins(openmp)]]

using namespace Rcpp;

// [[Rcpp::export]]
int intPower2(int a, int b){
  int res = 1;
  for (int i = 0; i < b; ++i) {
    res = a * res;
  }
  return res;
}


// [[Rcpp::export]]
arma::mat build_corr_from_params_cpp(int d, const arma::vec& params) {
  arma::mat R = arma::eye(d, d);
  int idx = 0;

  for (int i = 0; i < d; i++) {
    for (int j = i + 1; j < d; j++) {
      R(i, j) = params(idx);
      R(j, i) = params(idx);
      idx++;
    }
  }
  return R;
}


void generate_subsets(const std::vector<int>& elements,
                      int idx,
                      std::vector<int>& current,
                      std::vector< std::vector<int> >& out){
  if (idx == (int)elements.size()) {
    out.push_back(current);
    return;
  }

  // Not include element
  generate_subsets(elements, idx + 1, current, out);

  // Include element
  current.push_back(elements[idx]);
  generate_subsets(elements, idx + 1, current, out);
  current.pop_back();
}


// [[Rcpp::export]]
double gaussian_copula_cdf_cpp(const arma::vec& u,
                               const arma::mat& corrMat){
  int d = u.n_elem;

  arma::vec u_clip = u;
  double eps = 1e-12;
  for (int i = 0; i < d; i++) {
    if (u_clip(i) < eps)     u_clip(i) = eps;
    if (u_clip(i) > 1 - eps) u_clip(i) = 1 - eps;
  }

  // Probit transform
  NumericVector upper(d);
  for (int i = 0; i < d; i++) {
    upper[i] = R::qnorm(u_clip(i), 0.0, 1.0, true, false);
  }

  // Convert lower to numeric vector
  NumericVector lower(d, -std::numeric_limits<double>::infinity());

  NumericMatrix sigma = wrap(corrMat);
  NumericVector mean(d, 0.0);

  Environment mvtnorm = Environment::namespace_env("mvtnorm");
  Function pmvnorm = mvtnorm["pmvnorm"];

  SEXP result = pmvnorm(
    _["lower"] = lower,
    _["upper"] = upper,
    _["mean"]  = mean,
    _["sigma"] = sigma
  );

  return as<double>(result);
}

// [[Rcpp::export]]
double Memoisegaussian_copula_cdf_cpp(const arma::vec& u,
                                      const arma::mat& corrMat){
  // Get memoised function from namespace
  Environment env = Environment::namespace_env("SpatMet");
  Function f = env["Memoisegaussian_copula_cdf"];

  // Call memoised function
  return as<double>(
    f(_["u"] = u, _["corrMat"] = corrMat)
  );
}

// [[Rcpp::export]]
arma::mat JointTransitionMatrix_arma_cpp2(arma::mat gamma, int K) {
  int S = intPower2(2, K);
  arma::mat jointGamma(S, S, arma::fill::zeros);
  for (int a = 0; a < S; ++a) {
    for (int b = 0; b < S; ++b) {
      double prob = 1.0;
      for (int k = 0; k < K; ++k) {
        int from_k = (a / intPower2(2, k)) % 2;
        int to_k   = (b / intPower2(2, k)) % 2;
        prob *= gamma(from_k, to_k);
      }
      jointGamma(a, b) = prob;
    }
  }
  return jointGamma;
}

// [[Rcpp::export]]
arma::mat JointTransitionMatrix_per_strain_cpp2(List gamma_list, int K) {
  int S = intPower2(2, K);
  arma::mat jointGamma(S, S, arma::fill::zeros);
  for (int a = 0; a < S; ++a) {
    for (int b = 0; b < S; ++b) {
      double prob = 1.0;
      for (int k = 0; k < K; ++k) {
        int from_k = (a / intPower2(2, k)) % 2;
        int to_k   = (b / intPower2(2, k)) % 2;
        arma::mat currentGamma = gamma_list[k];
        prob *= currentGamma(from_k, to_k);
      }
      jointGamma(a, b) = prob;
    }
  }
  return jointGamma;
}

// [[Rcpp::export]]
arma::mat JointTransitionMatrix_copula_perstrain_cpp(List gamma_list,
                                           int K,
                                           const arma::vec& copParams){
  arma::mat corrMat = build_corr_from_params_cpp(K, copParams);
  int S = intPower2(2, K);
  arma::mat GammaMat(S, S, arma::fill::zeros);



  for (int a = 0; a < S; a++) {
    for (int b = 0; b < S; b++) {

      std::vector<int> Ones;
      std::vector<int> Zeros;
      arma::vec prob(K);

      for (int k = 0; k < K; k++) {
        arma::mat gamma = gamma_list[k];
        arma::mat gamma2 = gamma;
        gamma2(0,0) = gamma(0, 1);
        gamma2(0,1) = gamma(0, 0);
        int from_k = (a >> k) & 1;
        int to_k   = (b >> k) & 1;

        if (from_k == 1)
          Ones.push_back(k);
        else
          Zeros.push_back(k);

        prob(k) = gamma2(from_k, to_k);
      }

      // power set of zeros
      std::vector< std::vector<int> > subsets;
      std::vector<int> cur;
      generate_subsets(Zeros, 0, cur, subsets);

      double total = 0.0;

      for (auto& Tset : subsets) {
        int sign = (Tset.size() % 2 == 0 ? 1 : -1);

        arma::vec u(K, arma::fill::ones);
        for (int idx : Ones) u(idx) = prob(idx);
        for (int idx : Tset) u(idx) = prob(idx);

        total += sign *  Memoisegaussian_copula_cdf_cpp(u, corrMat);
      }
      GammaMat(a,b) = total;
    }
  }

  for (int i = 0; i < S; i++) {
    double s = arma::accu(GammaMat.row(i));
    GammaMat.row(i) /= s;
  }
  return GammaMat;
}

// [[Rcpp::export]]
arma::mat JointTransitionMatrix_copula_cpp(const arma::mat& gamma,
                                           int K,
                                           const arma::vec& copParams){
  arma::mat corrMat = build_corr_from_params_cpp(K, copParams);
  int S = intPower2(2, K);
  arma::mat GammaMat(S, S, arma::fill::zeros);

  arma::mat gamma2 = gamma;
  gamma2(0,0) = gamma(0, 1);
  gamma2(0,1) = gamma(0, 0);

  for (int a = 0; a < S; a++) {
    for (int b = 0; b < S; b++) {

      std::vector<int> Ones;
      std::vector<int> Zeros;
      arma::vec prob(K);

      for (int k = 0; k < K; k++) {
        int from_k = (a >> k) & 1;
        int to_k   = (b >> k) & 1;

        if (from_k == 1)
          Ones.push_back(k);
        else
          Zeros.push_back(k);

        prob(k) = gamma2(from_k, to_k);
      }

      // power set of zeros
      std::vector< std::vector<int> > subsets;
      std::vector<int> cur;
      generate_subsets(Zeros, 0, cur, subsets);

      double total = 0.0;

      for (auto& Tset : subsets) {
        int sign = (Tset.size() % 2 == 0 ? 1 : -1);

        arma::vec u(K, arma::fill::ones);
        for (int idx : Ones) u(idx) = prob(idx);
        for (int idx : Tset) u(idx) = prob(idx);

        total += sign *  Memoisegaussian_copula_cdf_cpp(u, corrMat);
      }
      GammaMat(a,b) = total;
    }
  }

  for (int i = 0; i < S; i++) {
    double s = arma::accu(GammaMat.row(i));
    GammaMat.row(i) /= s;
  }
  return GammaMat;
}


// [[Rcpp::export]]
arma::mat makematrix_arma_cpp2(double g12, double g21){
  arma::mat Gmat(2, 2, arma::fill::zeros);
  Gmat(0, 0) = 1 - g12;
  Gmat(0, 1) = g12;
  Gmat(1, 0) = g21;
  Gmat(1, 1) = 1 - g21;
  return Gmat;
}

// [[Rcpp::export]]
Rcpp::List BuildGamma_list_cpp(const arma::vec& Gs) {
  int n_mat = Gs.n_elem / 2;
  Rcpp::List out(n_mat);

  int index = 0;
  for (int i = 0; i < n_mat; i++) {
    double a = Gs[i + index];
    double b = Gs[i + index + 1];

    out[i] = makematrix_arma_cpp2(a, b);
    index += 1;
  }
  return out;
}

// [[Rcpp::export]]
arma::mat Multipurpose_JointTransitionMatrix_cpp(const arma::vec& gammas,
                                           int K,
                                           const arma::vec& copParams, int Modeltype){
  int nstate = intPower2(2, K);
  arma::mat JointTPM(nstate, nstate, arma::fill::zeros);
  if(Modeltype == 1){
    arma::mat Gmat = makematrix_arma_cpp2(gammas(0), gammas(1));
    JointTPM = JointTransitionMatrix_arma_cpp2(Gmat, K);
  }else if(Modeltype == 2){
    List Glist = BuildGamma_list_cpp(gammas);
    JointTPM = JointTransitionMatrix_per_strain_cpp2(Glist, K);
  }else if(Modeltype == 3){
    arma::mat Gmat = makematrix_arma_cpp2(gammas(0), gammas(1));
    JointTPM = JointTransitionMatrix_copula_cpp(Gmat, K, copParams);
  }else if(Modeltype == 4){
    List Glist = BuildGamma_list_cpp(gammas);
    JointTPM = JointTransitionMatrix_copula_perstrain_cpp(Glist, K, copParams);
  }
  return JointTPM;
}


// [[Rcpp::export]]
double pmvnorm_cpp(
    arma::vec const& x, arma::vec const& mean, arma::mat const& Sigma,
    double abseps = 1e-3){
  int n = x.n_elem;
  arma::vec sds = arma::sqrt(Sigma.diag());
  arma::mat corrmat = arma::diagmat(1.0 / sds) * Sigma * arma::diagmat(1.0 / sds);
  arma::vec bound = (x - mean) / sds;
  arma::vec lowertrivec(n * (n - 1) / 2);
  int k = 0;
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < i; ++j) {
      lowertrivec(k++) = corrmat(i, j);
    }
  }
  int nu = 0;
  int maxpts = 25000;
  double releps = 0;
  int rnd = 1;
  double* bound_ = bound.memptr();
  double* correlationMatrix = lowertrivec.memptr();
  double* lower = new double[n];
  int* infin = new int[n];
  double* delta = new double[n];
  for (int i = 0; i < n; ++i) {
    infin[i] = 0;
    lower[i] = 0.0;
    delta[i] = 0.0;
  }
  double error;
  double value;
  int inform;
  mvtnorm_C_mvtdst(
    &n, &nu, lower, bound_, infin, correlationMatrix, delta, &maxpts, &abseps,
    &releps, &error, &value, &inform, &rnd
  );
  delete[] (lower);
  delete[] (infin);
  delete[] (delta);
  return value;
}

// [[Rcpp::export]]
double gaussian_copula_cdf_cpp2(const arma::vec& u,
                               const arma::mat& corrMat){
  int d = u.n_elem;

  arma::vec u2 = u;
  double eps = 1e-12;
  for (int i = 0; i < d; i++) {
    if (u2(i) < eps) u2(i) = eps;
    if (u2(i) > 1 - eps) u2(i) = 1 - eps;
  }

  arma::vec x(d);
  for (int i = 0; i < d; i++)
    x(i) = R::qnorm(u2(i), 0, 1, 1, 0);
  arma::vec mean(d, arma::fill::zeros);
  return pmvnorm_cpp(x, mean, corrMat);
}


// [[Rcpp::export]]
arma::mat ParallelJointTransitionMatrix_copula_cpp(const arma::mat& gamma,
                                                   int K,
                                                   const arma::vec& copParams){
  arma::mat corrMat = build_corr_from_params_cpp(K, copParams);
  int S = intPower2(2, K);
  arma::mat GammaMat(S, S, arma::fill::zeros);

  arma::mat gamma2 = gamma;
  gamma2(0,0) = gamma(0, 1);
  gamma2(0,1) = gamma(0, 0);

  // parallelize outer loop
  RcppThread::parallelFor(0, S, [&](int a){
        for (int b = 0; b < S; b++) {
          std::vector<int> Ones;
          std::vector<int> Zeros;
          arma::vec prob(K);

          for (int k = 0; k < K; k++) {
            int from_k = (a >> k) & 1;
            int to_k   = (b >> k) & 1;
            if (from_k == 1)
              Ones.push_back(k);
            else
              Zeros.push_back(k);
            prob(k) = gamma2(from_k, to_k);
          }
    // Generate power set of Zeros
          std::vector<std::vector<int>> subsets;
          std::vector<int> cur;
          generate_subsets(Zeros, 0, cur, subsets);

          double total = 0.0;
          for (auto &Tset : subsets) {
            int sign = (Tset.size() % 2 == 0 ? 1 : -1);
            arma::vec u(K, arma::fill::ones);
            for (int idx : Ones) u(idx) = prob(idx);
            for (int idx : Tset) u(idx) = prob(idx);
            total += sign * gaussian_copula_cdf_cpp2(u, corrMat);
          }
          GammaMat(a, b) = total;
        }
  });

    for (int i = 0; i < S; i++) {
      double s = arma::accu(GammaMat.row(i));
      GammaMat.row(i) /= s;
    }
    return GammaMat;
}

// [[Rcpp::export]]
arma::mat ParallelJointTransitionMatrix_copula_cpp2(const arma::mat& gamma,
                                           int K,
                                           const arma::vec& copParams){
  arma::mat corrMat = build_corr_from_params_cpp(K, copParams);
  int S = intPower2(2, K);
  arma::mat GammaMat(S, S, arma::fill::zeros);

  arma::mat gamma2 = gamma;
  gamma2(0,0) = gamma(0, 1);
  gamma2(0,1) = gamma(0, 0);

  // parallelize
#pragma omp parallel for schedule(dynamic)
  for (int a = 0; a < S; a++) {
    for (int b = 0; b < S; b++) {

      std::vector<int> Ones;
      std::vector<int> Zeros;
      arma::vec prob(K);

      for (int k = 0; k < K; k++) {
        int from_k = (a >> k) & 1;
        int to_k   = (b >> k) & 1;

        if (from_k == 1)
          Ones.push_back(k);
        else
          Zeros.push_back(k);

        prob(k) = gamma2(from_k, to_k);
      }

      // power set of zeros
      std::vector<std::vector<int>> subsets;
      std::vector<int> cur;
      generate_subsets(Zeros, 0, cur, subsets);

      double total = 0.0;

      for (auto& Tset : subsets) {
        int sign = (Tset.size() % 2 == 0 ? 1 : -1);

        arma::vec u(K, arma::fill::ones);
        for (int idx : Ones) u(idx) = prob(idx);
        for (int idx : Tset) u(idx) = prob(idx);

        total += sign * gaussian_copula_cdf_cpp2(u, corrMat);
      }

      GammaMat(a, b) = total;
    }
  }

  for (int i = 0; i < S; i++) {
    double s = arma::accu(GammaMat.row(i));
    GammaMat.row(i) /= s;
  }
  return GammaMat;
}


// [[Rcpp::export]]
double frank_cdf_cpp(arma::vec u, double Psi){
  double cdf = 0;
  if(Psi == 0.0){
    cdf = arma::prod(u);
  }else{
  int d = u.n_elem;
  arma::vec tmp = 1.0 - arma::exp(-Psi * u);
  double num = arma::prod(tmp);
  double den = std::pow(1.0 - std::exp(-Psi), d - 1);
  cdf = -(1.0 / Psi) * std::log(1.0 - num / den);
  }
  return cdf;
}

// [[Rcpp::export]]
double frank_cdf_cpp2(const arma::vec &u, double Psi) {
  int d = u.n_elem;

  if (std::abs(Psi) < 1e-12) {
    double prod = 1.0;
    for (int i = 0; i < d; i++) {
      double ui = u[i];
      if (ui <= 0.0) return 0.0;
      if (ui >= 1.0) ui = 1.0;
      prod *= ui;
    }
    return prod;
  }

  for (int i = 0; i < d; i++) {
    if (!std::isfinite(u[i])) return NA_REAL;
    if (u[i] < 0.0 || u[i] > 1.0) return NA_REAL;
  }

  double den = std::expm1(-Psi);
  if (den == 0.0) return NA_REAL;

  long double prod_term = 1.0L;

  for (int i = 0; i < d; i++) {
    double ui = u[i];
    double term = std::expm1(-Psi * ui);
    prod_term *= (long double)term;
  }

  long double ratio = prod_term / std::pow(den, d - 1);

  double result = -(1.0 / Psi) * std::log1p((double)ratio);

  if (!std::isfinite(result)) return NA_REAL;
  if (result < 0.0) result = 0.0;
  if (result > 1.0) result = 1.0;

  return result;
}

// [[Rcpp::export]]
arma::mat ParallelFrankJointTransitionMatrix_copula_cpp(const arma::mat& gamma,
                                                    int K,
                                                    double copParam){
  int S = intPower2(2, K);
  arma::mat GammaMat(S, S, arma::fill::zeros);

  arma::mat gamma2 = gamma;
  gamma2(0,0) = gamma(0, 1);
  gamma2(0,1) = gamma(0, 0);

  // parallelize
#pragma omp parallel for schedule(dynamic)
  for (int a = 0; a < S; a++) {
    for (int b = 0; b < S; b++) {

      std::vector<int> Ones;
      std::vector<int> Zeros;
      arma::vec prob(K);

      for (int k = 0; k < K; k++) {
        int from_k = (a >> k) & 1;
        int to_k   = (b >> k) & 1;

        if (from_k == 1)
          Ones.push_back(k);
        else
          Zeros.push_back(k);

        prob(k) = gamma2(from_k, to_k);
      }

      // power set of zeros
      std::vector<std::vector<int>> subsets;
      std::vector<int> cur;
      generate_subsets(Zeros, 0, cur, subsets);

      double total = 0.0;

      for (auto& Tset : subsets){
        int sign = (Tset.size() % 2 == 0 ? 1 : -1);

        arma::vec u(K, arma::fill::ones);
        for (int idx : Ones) u(idx) = prob(idx);
        for (int idx : Tset) u(idx) = prob(idx);

        total += sign * frank_cdf_cpp2(u, copParam);
      }

      GammaMat(a, b) = total;
    }
  }

  for (int i = 0; i < S; i++) {
    double s = arma::accu(GammaMat.row(i));
    GammaMat.row(i) /= s;
  }
  return GammaMat;
}


// [[Rcpp::export]]
arma::mat ParallelFrankJointTransitionMatrix_copula_perstrain_cpp(List gamma_list,
                                                        int K,
                                                        double copParam){
  int S = intPower2(2, K);
  arma::mat GammaMat(S, S, arma::fill::zeros);

  std::vector<arma::mat> gamma_cpp(K);
  for (int k = 0; k < K; k++)
    gamma_cpp[k] = as<arma::mat>(gamma_list[k]);

  // parallelize
#pragma omp parallel for schedule(dynamic)
  for (int a = 0; a < S; a++) {
    for (int b = 0; b < S; b++) {

      std::vector<int> Ones;
      std::vector<int> Zeros;
      arma::vec prob(K);

      for (int k = 0; k < K; k++) {
        arma::mat gamma = gamma_cpp[k];
        arma::mat gamma2 = gamma;
        gamma2(0,0) = gamma(0, 1);
        gamma2(0,1) = gamma(0, 0);
        int from_k = (a >> k) & 1;
        int to_k   = (b >> k) & 1;

        if (from_k == 1)
          Ones.push_back(k);
        else
          Zeros.push_back(k);

        prob(k) = gamma2(from_k, to_k);
      }

      // power set of zeros
      std::vector<std::vector<int>> subsets;
      std::vector<int> cur;
      generate_subsets(Zeros, 0, cur, subsets);

      double total = 0.0;

      for (auto& Tset : subsets) {
        int sign = (Tset.size() % 2 == 0 ? 1 : -1);

        arma::vec u(K, arma::fill::ones);
        for (int idx : Ones) u(idx) = prob(idx);
        for (int idx : Tset) u(idx) = prob(idx);

        total += sign * frank_cdf_cpp2(u, copParam);
      }

      GammaMat(a, b) = total;
    }
  }

  for (int i = 0; i < S; i++) {
    double s = arma::accu(GammaMat.row(i));
    GammaMat.row(i) /= s;
  }
  return GammaMat;
}

// [[Rcpp::export]]
arma::mat Multipurpose_JointTransitionMatrix_cpp2(const arma::vec& gammas,
                                                 int K,
                                                 double copParam, int Modeltype){
  int nstate = intPower2(2, K);
  arma::mat JointTPM(nstate, nstate, arma::fill::zeros);
  if(Modeltype == 1){
    arma::mat Gmat = makematrix_arma_cpp2(gammas(0), gammas(1));
    JointTPM = JointTransitionMatrix_arma_cpp2(Gmat, K);
  }else if(Modeltype == 2){
    List Glist = BuildGamma_list_cpp(gammas);
    JointTPM = JointTransitionMatrix_per_strain_cpp2(Glist, K);
  }else if(Modeltype == 3){
    arma::mat Gmat = makematrix_arma_cpp2(gammas(0), gammas(1));
    JointTPM = ParallelFrankJointTransitionMatrix_copula_cpp(Gmat, K, copParam);
  }else if(Modeltype == 4){
    List Glist = BuildGamma_list_cpp(gammas);
    JointTPM = ParallelFrankJointTransitionMatrix_copula_perstrain_cpp(Glist, K, copParam);
  }
  return JointTPM;
}

