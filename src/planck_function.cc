#include <cmath>

#include "constants.hh"

double planck_function(const double lambda, const double temperature) {
    return (2.0 * h_planck * std::pow(c_light,2) / std::pow(lambda,5)) / (std::exp((h_planck * c_light) / (lambda * k_boltzmann * temperature)) - 1.0);
}
