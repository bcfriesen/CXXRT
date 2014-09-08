#include <cmath>

#include "../constants.hh"
#include "../EOS/atoms.hh"

double Doppler_width(const double lambda_0, const double temperature, const Atom &atom) {
    return std::sqrt((8.0 * k_boltzmann * temperature * std::log(2.0)) / (atom.atomic_weight * std::pow(c_light, 2))) * lambda_0;
}
