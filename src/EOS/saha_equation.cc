#include <iostream>
#include <cmath>
#include <limits>

#include "atoms.hh"
#include "../constants.hh"

double saha_equation(Atom* atom, const unsigned int lower_ionization_stage, const double n_e, const double temperature) {

    double result = 0.0;
    const double g_e = 2.0;

    std::vector<Ion>::const_iterator lower_ion_it;
    std::vector<Ion>::const_iterator upper_ion_it;

    // This is a very slight (and totally worth it) price to pay for exchanging
    // array indices for iterators: we can no longer assume that the i-th
    // ionization stage has index i in the "ions" array, and instead we must
    // search for the ionization stage itself, and save the iterator which
    // points to it. Ultimately this works in our favor because this method is
    // immune to indexing errors which can arise if elements of the ions array
    // are added or removed.

    for (lower_ion_it = atom->ions.begin(); lower_ion_it != atom->ions.end(); ++lower_ion_it) {
        if (lower_ion_it->ionization_stage == lower_ionization_stage) break;
    }
    for (upper_ion_it = atom->ions.begin(); upper_ion_it != atom->ions.end(); ++upper_ion_it) {
        if (upper_ion_it->ionization_stage == lower_ionization_stage+1) break;
    }

    double lower_partition_function = lower_ion_it->partition_function;
    double upper_partition_function;

    // If the upper ionization stage is the fully ionized atom, then it will have no levels so skip looking for the partition function because it's just 1.
    if (lower_ionization_stage == atom->atomic_number-1) {
        upper_partition_function = 1.0;
    } else {
        upper_partition_function = upper_ion_it->partition_function;
    }

    result = (1.0 / n_e) * (std::pow(2.0 * pi * m_electron * k_boltzmann * temperature, 1.5) / std::pow(h_planck, 3)) * g_e * (upper_partition_function / lower_partition_function) * std::exp(-lower_ion_it->ionization_potential / (k_boltzmann * temperature));

    return result;
}
