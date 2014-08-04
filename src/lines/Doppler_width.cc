#include <cmath>

#include "../constants.hh"

double Doppler_width(const double lambda_0, const double temperature) {
    return std::sqrt((8.0 * k_boltzmann * temperature * std::log(2.0)) / (m_electron * std::pow(c_light, 2))) * (c_light / lambda_0);
}
