#ifndef CALC_ALO_HH
#define CALC_ALO_HH

#include <map>

#include <Eigen/Dense>

Eigen::SparseMatrix<double> calc_ALO(const std::size_t wl_value_hash);

#endif
