#include <cmath>

#include "../constants.hh"

double lorentz_profile(const double lambda, const double lambda_0, const double Delta_lambda) {
    return (1.0 / pi) * (0.5 * Delta_lambda) / (std::pow(lambda - lambda_0, 2) + std::pow(0.5 * Delta_lambda, 2))
}
