#include "utils.h"
#include <stdexcept>

//' Evaluate metrics
//' @name evaluate
//' @param num_Sim A numeric vector of simulated values.
//' @param num_Obs A numeric vector of observed values. NA values are removed along with corresponding values in num_Sim.
//' @return A double representing the Evaluate metrics.
//' @export
double evalute_NSE(arma::vec num_Sim, arma::vec num_Obs) {
  // Ensure both vectors have the same length
  if (num_Sim.n_elem != num_Obs.n_elem) {
    throw std::invalid_argument("Simulated and observed vectors must have the same length.");
  }
  
  // Find non-NA indices in observed data
  arma::uvec non_na_indices = find_finite(num_Obs);
  
  // Filter out NA values
  arma::vec filtered_Obs = num_Obs.elem(non_na_indices);
  arma::vec filtered_Sim = num_Sim.elem(non_na_indices);
  
  // Ensure there are remaining values to compute NSE
  if (filtered_Obs.n_elem == 0) {
    throw std::invalid_argument("All observed values are NA; cannot calculate NSE.");
  }
  
  // Calculate the mean of observed values
  double obs_mean = arma::mean(filtered_Obs);
  
  // Compute numerator and denominator for NSE
  double numerator = arma::accu(arma::square(filtered_Obs - filtered_Sim));
  double denominator = arma::accu(arma::square(filtered_Obs - obs_mean));
  
  // Return NSE
  return 1.0 - (numerator / denominator);
}

//' @rdname evaluate
//' @param factor_r,factor_alpha,factor_beta A double specifying the weight for the correlation term (r - 1), (alpha - 1) and (beta - 1). Default is 1.0.
//' @export
double evalute_KGE(arma::vec num_Sim, arma::vec num_Obs,
                   double factor_r = 1.0, double factor_alpha = 1.0, double factor_beta = 1.0) {
  // Ensure both vectors have the same length
  if (num_Sim.n_elem != num_Obs.n_elem) {
    throw std::invalid_argument("Simulated and observed vectors must have the same length.");
  }
  
  // Find non-NA indices in observed data
  arma::uvec non_na_indices = find_finite(num_Obs);
  
  // Filter out NA values
  arma::vec filtered_Obs = num_Obs.elem(non_na_indices);
  arma::vec filtered_Sim = num_Sim.elem(non_na_indices);
  
  // Ensure there are remaining values to compute KGE
  if (filtered_Obs.n_elem == 0) {
    throw std::invalid_argument("All observed values are NA; cannot calculate KGE.");
  }
  
  // Calculate means of observed and simulated values
  double obs_mean = arma::mean(filtered_Obs);
  double sim_mean = arma::mean(filtered_Sim);
  
  // Calculate standard deviations
  double obs_sd = arma::stddev(filtered_Obs);
  double sim_sd = arma::stddev(filtered_Sim);
  
  // Calculate correlation coefficient - manual calculation
  double numerator_corr = 0.0;
  for (size_t i = 0; i < filtered_Sim.n_elem; i++) {
    numerator_corr += (filtered_Obs[i] - obs_mean) * (filtered_Sim[i] - sim_mean);
  }
  double correlation = numerator_corr / ((filtered_Sim.n_elem - 1) * obs_sd * sim_sd);
  
  // Alternative: using Armadillo's cov function and extracting scalar value
  // arma::mat corr_matrix = arma::cor(filtered_Obs, filtered_Sim);
  // double correlation = corr_matrix(0, 0);
  
  // Calculate KGE components
  double alpha = sim_sd / obs_sd;
  double beta = sim_mean / obs_mean;
  
  // Apply factors to each component
  double sum_Factor = factor_r + factor_alpha + factor_beta;
  factor_r = factor_r / sum_Factor;
  factor_alpha = factor_alpha / sum_Factor;
  factor_beta = factor_beta / sum_Factor;
  
  double term_r = pow(factor_r * (correlation - 1.0), 2);
  double term_alpha = pow(factor_alpha * (alpha - 1.0), 2);
  double term_beta = pow(factor_beta * (beta - 1.0), 2);
  
  // Calculate and return KGE
  return 1.0 - sqrt(term_r + term_alpha + term_beta);
}