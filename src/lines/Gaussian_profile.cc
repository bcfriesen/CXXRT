#include <cmath>

#include "../constants.hh"

double gauss_profile(const double lambda, const double lambda_0, const double Delta_lambda) {
    return (1.0 / (Delta_lambda * std::sqrt(2.0 * pi))) * std::exp(-std::pow(lambda - lambda_0, 2) / (2.0 * std::pow(Delta_lambda, 2)));
}
