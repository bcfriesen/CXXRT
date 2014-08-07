#include <cmath>

#include "atoms.hh"
#include "../constants.hh"

// These functions are useful for calculating LTE populations of particular
// levels once the LTE equation of state has solved for n_e.

// Mihalas defines the Phi function in his Eq (5-14)
double Phi(const AtomicLevel level, const Ion ion, const Atom atom, const double temperature) {
    const double C_I = 0.5 * std::pow(std::pow(h_planck, 2) / (2.0 * pi * m_electron * k_boltzmann), 3.0/2.0); // Mihalas defines this constant in his Eq (5-14)

    return (level.g / ion.continuum_state->g) * C_I * std::pow(temperature, -3.0/2.0) * std::exp((ion.ionization_potential - level.energy) / (k_boltzmann * temperature));
}


// Mihalas defines the Phi_tilde function in his Eq (5-15)
double Phi_tilde(const AtomicLevel level, const Ion ion, const Atom atom, const double temperature) {
    return (ion.continuum_state->g / ion.partition_function) * Phi(level, ion, atom, temperature);
}
