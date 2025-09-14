#include <RcppArmadillo.h>
#include <random>
#include <vector>
#include <cmath>
#include <stdexcept>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::interfaces(r, cpp)]]

using Rcpp::List;

// Helper function to generate a binary vector with probability P_i
static std::vector<bool> generate_binary_vector(int n_x, double P_i, std::default_random_engine &generator) {
  std::bernoulli_distribution distribution(P_i);
  std::vector<bool> binary_vector(n_x);
  for (int i = 0; i < n_x; i++) {
    binary_vector[i] = distribution(generator);
  }
  return binary_vector;
}

//' Calibrate
//' This function implements a calibration algorithm based on the Direct Search method (DDS).
//'
//' It attempts to find the optimal parameter values that minimize the given objective (fitness) function.
//'
//' @name calibrate
//' @param fitness An Rcpp-exported C++ function (callable from R and C++) that accepts (x, other_data) and returns arma::vec (first element is objective).
//' @param lst_OtherData A field of additional data passed to `fitness`.
//' @param x_Min A numeric vector of lower bounds for each parameter.
//' @param x_Max A numeric vector of upper bounds for each parameter.
//' @param x_Init An optional initial solution (numeric vector). If empty, midpoint of x_Min and x_Max is used.
//' @param max_iter Maximum iterations (default 100).
//' @param r Perturbation factor (default 0.2).
//'
//' @return arma::vec: best parameter set found.
//'
arma::vec cali_DDS(arma::vec (*fitness)(const arma::vec&, const arma::field<arma::vec>&),
                   const arma::field<arma::vec>& lst_OtherData,
                   const arma::vec& x_Min,
                   const arma::vec& x_Max,
                   arma::vec x_Init = arma::vec(),
                   int max_iter = 100,
                   double r = 0.2) {
  if (x_Min.n_elem != x_Max.n_elem) {
    throw std::invalid_argument("x_Min and x_Max must have the same length");
  }
  if (max_iter <= 0) {
    throw std::invalid_argument("max_iter must be > 0");
  }
  if (r <= 0.0) {
    throw std::invalid_argument("r must be > 0");
  }

  int n_x = x_Min.n_elem;

  // Initialize x_Best
  arma::vec x_Best = x_Init.is_empty() ? (x_Min + x_Max) / 2.0 : x_Init;
  if ((int)x_Best.n_elem != n_x) {
    throw std::invalid_argument("x_Init must have the same length as x_Min/x_Max");
  }

  // Evaluate initial fitness
  arma::vec fitness_result0 = fitness(x_Best, lst_OtherData);
  if (fitness_result0.n_elem < 1) {
    throw std::runtime_error("fitness must return a vector with at least one element (objective)");
  }
  double y_Best = fitness_result0[0];

  // Precompute perturbation probabilities P_i
  arma::vec P_i(max_iter);
  for (int i = 0; i < max_iter; i++) {
    P_i[i] = 1.0 - std::log(i + 1.0) / std::log(max_iter);
    if (P_i[i] < 0.0) P_i[i] = 0.0;
    if (P_i[i] > 1.0) P_i[i] = 1.0;
  }

  std::random_device rd;
  std::default_random_engine generator(rd());
  std::normal_distribution<double> normal_dist(0.0, 1.0);

  std::vector<std::vector<int>> lst_Cali_x(max_iter);
  for (int i = 0; i < max_iter; i++) {
    std::vector<bool> randm_Para = generate_binary_vector(n_x, P_i[i], generator);
    std::vector<int> selected_indices;
    for (int j = 0; j < n_x; j++) {
      if (randm_Para[j]) selected_indices.push_back(j);
    }
    if (selected_indices.empty()) {
      std::uniform_int_distribution<int> uniform_dist(0, n_x - 1);
      selected_indices.push_back(uniform_dist(generator));
    }
    lst_Cali_x[i] = std::move(selected_indices);
  }

  arma::vec sigma_ = x_Max - x_Min;

  Rcpp::Rcout << "Calibration in progress...\n";
  int progress_step = std::max(1, max_iter / 10);

  for (int iter = 1; iter < max_iter; iter++) {
    arma::vec x_New = x_Best;

    arma::vec N_01(n_x);
    for (int j = 0; j < n_x; j++) N_01[j] = normal_dist(generator);

    arma::vec x_New0 = x_Best + r * N_01 % sigma_;

    arma::vec x_New1 = arma::max(2.0 * x_Min - x_New0, x_Min);
    arma::vec x_New2 = arma::min(2.0 * x_Max - x_New0, x_Max);

    for (int j = 0; j < n_x; j++) {
      if (x_New0[j] < x_Min[j]) x_New0[j] = x_New1[j];
      if (x_New0[j] > x_Max[j]) x_New0[j] = x_New2[j];
    }

    for (int idx : lst_Cali_x[iter]) {
      x_New[idx] = x_New0[idx];
    }

    arma::vec fitness_result_new = fitness(x_New, lst_OtherData);
    if (fitness_result_new.n_elem < 1) {
      throw std::runtime_error("fitness must return a vector with at least one element (objective)");
    }
    double y_New = fitness_result_new[0];

    if (y_New < y_Best) {
      x_Best = x_New;
      y_Best = y_New;
    }

    if ((iter % progress_step) == 0) {
      Rcpp::Rcout << "Iteration " << iter << " / " << max_iter << " completed. Current best = " << y_Best << "\n";
    }
  }

  return x_Best;
}
