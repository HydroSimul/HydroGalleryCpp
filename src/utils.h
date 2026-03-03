#ifndef UTILS
#define UTILS
#define ARMA_DONT_USE_BLAS
#define ARMA_DONT_USE_LAPACK
#include <armadillo>                    // now uses malloc/free



#include <string>
#include <set>
#include <algorithm>

arma::mat load_mat(const std::string& path);
arma::vec load_vec(const std::string& path);

#endif // UTILS
