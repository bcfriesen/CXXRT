#include <cmath>

#include "planck_function.hh"

// Use centered difference formula since dB_lambda/dT is not analytic.
double planck_function_temperature_derivative(const double lambda, const double temperature) {
    const double delta_T = 1.0e-8;
    const double forward_delta_B = planck_function(lambda, temperature + delta_T);
    const double backward_delta_B = planck_function(lambda, temperature - delta_T);

    return (forward_delta_B - backward_delta_B) / (2.0 * delta_T);
}
