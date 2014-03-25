#include <cmath>

#include "constants.hh"

double planck_function(const double lambda, const double temperature) {
  return (2.0 * h_planck * pow(c_light,2) / pow(lambda,5)) / (exp((h_planck * c_light) / (lambda * k_boltzmann * temperature)) - 1.0);
}
