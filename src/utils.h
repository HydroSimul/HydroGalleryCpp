#ifndef UTILS
#define UTILS

// #include <RcppArmadillo.h>
#define ARMA_ALIEN_MEM_ALLOC_FUNCTION   // define as empty → disables hijack
#define ARMA_ALIEN_MEM_FREE_FUNCTION    // same
#define ARMA_DONT_USE_BLAS
#define ARMA_DONT_USE_LAPACK
#include <armadillo>                    // now uses malloc/free
#include <Rcpp.h>                       // R package glue still works



#include <string>
#include <set>
#include <algorithm>

arma::mat load_mat(const std::string& path);
arma::vec load_vec(const std::string& path);

#endif // UTILS
