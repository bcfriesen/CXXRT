#include <cstdio>
#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <limits>

#include <gsl/gsl_sf.h>

#include "atoms.hh"
#include "../globals.hh"
#include "../constants.hh"
#include "../lines/Gaussian_profile.hh"
#include "../lines/Doppler_width.hh"
#include "../wavelength_grid.hh"

Atom::Atom(const unsigned int atomic_number_in)
    : atomic_number(atomic_number_in) {
    const std::string atomic_data_file_name = config["atomic_data_root_folder"].as<std::string>() + "/" + "atoms.dat";
    std::ifstream atomic_data_file;
    std::string one_line;
    atomic_data_file.open(atomic_data_file_name.c_str()); // in C++11 the argument of open() can be a string
    unsigned int atomic_number_from_file;
    while (std::getline(atomic_data_file, one_line)) {
        std::istringstream iss(one_line);
        iss >> atomic_number_from_file;
        if (atomic_number_from_file == atomic_number)
            iss >> atomic_weight >> atomic_symbol;
    }
    atomic_data_file.close();

    for (unsigned int i = 0; i <= atomic_number_in; ++i) {
        Ion ion(atomic_number, i);
        ions.push_back(ion);
    }
}


// We assume all continuum processes promote a bound state to the ground state
// of the next ionization stage (as opposed to promotion from the bound state
// to an excited state of the next ionization stage).
void Atom::set_continuum_pointers() {
    for (auto &ion: ions) {

        // A fully ionized atom has no higher continuum state.
        if (ion.ionization_stage == atomic_number)
            ion.continuum_state = nullptr;

        for (auto &next_ion: ions) {
            if (next_ion.ionization_stage == ion.ionization_stage+1) {
                ion.next_ion = &next_ion;
                for (auto &level: next_ion.levels) {
                    if (level.energy < std::numeric_limits<double>::epsilon()) {
                        ion.continuum_state = &level;
                        break;
                    }
                }
            }
        }
    }
}


Ion::Ion(const unsigned int atomic_number_in, const unsigned int ionization_stage_in)
    : atomic_number(atomic_number_in),
      ionization_stage(ionization_stage_in) {
}


void Ion::read_atomic_data() {
    if (ionization_stage > atomic_number) {
        std::cerr << "ERROR! in atom " << atomic_number << " you tried to add an ion with ionization stage " << ionization_stage << std::endl;
        exit(1);
    } else if (ionization_stage == atomic_number) { // A fully ionized atom is just a nucleus which has no bound states, only a ground state;
        class AtomicLevel ground_state;
        ground_state.energy = 0.0;
        ground_state.J = 0.0;
        ground_state.g = 1;
        levels.push_back(ground_state);
    } else {
        std::ostringstream convert; // for compilers without std::to_string (a C++11 feature)
        convert << atomic_number*100 + ionization_stage;
        const std::string atomic_data_file_name = config["atomic_data_root_folder"].as<std::string>() + "/" + convert.str() + ".dat";
        std::ifstream atomic_data_file;
        atomic_data_file.open(atomic_data_file_name.c_str()); // in C++11 the argument of open() can be a string
        std::string one_line;
        double wavelength, log_gf, first_energy_level, J_first, second_energy_level, J_second;
        bool duplicate;
        const double tolerance = 1.0e-30;

        // The first line is the ionization potential.
        std::getline(atomic_data_file, one_line);
        std::istringstream iss(one_line);
        iss >> ionization_potential;
        ionization_potential *= h_planck * c_light;

        while (std::getline(atomic_data_file, one_line)) {
            std::istringstream iss(one_line);
            iss >> wavelength >> log_gf >> first_energy_level >> J_first >> second_energy_level >> J_second;

            class AtomicLevel upper_level, lower_level;
            lower_level.energy = first_energy_level * h_planck * c_light;  // convert cm^-1 to erg
            lower_level.J = J_first;
            lower_level.g = int(2.0*lower_level.J + 1.0);
            upper_level.energy = second_energy_level * h_planck * c_light; // convert cm^-1 to erg
            upper_level.J = J_second;
            upper_level.g = int(2.0*upper_level.J + 1.0);

            // Add the lower level of the line to the model atom only if it hasn't been
            // added already. We detect duplicates by checking energy differences; if
            // two levels are similar to within some tolerance, they are defined as
            // duplicates.

            for (auto level: levels) {
                double difference;
                if (level.energy < std::numeric_limits<double>::epsilon()) {
                    difference = std::abs(lower_level.energy - level.energy);
                } else {
                    difference = std::abs(lower_level.energy - level.energy) / level.energy;
                }
                if (difference > tolerance) {
                    duplicate = false;
                } else {
                    duplicate = true;
                    break;
                }
            }
            if (!duplicate) {
                levels.push_back(lower_level);
            }
            // Now add the upper level if necessary.
            for (auto level: levels) {
                const double difference = std::abs(upper_level.energy - level.energy) / level.energy;
                if (difference > tolerance) {
                    duplicate = false;
                } else {
                    duplicate = true;
                    break;
                }
            }
            if (!duplicate) {
                levels.push_back(upper_level);
            }
        }

        // Set up the pointer for the ground state. By definition the energy of
        // the ground state is zero.
        for (auto &level: levels) {
            if (level.energy < std::numeric_limits<double>::epsilon())
                ground_state = &level;
        }

        atomic_data_file.close();

        // Now read the file again to get the lines.

        atomic_data_file.open(atomic_data_file_name.c_str()); // in C++ the argument of open() can be a string

        // The first line is the ionization potential so skip it this time.
        std::getline(atomic_data_file, one_line);

        while (std::getline(atomic_data_file, one_line)) {
            std::istringstream iss(one_line);
            iss >> wavelength >> log_gf >> first_energy_level >> J_first >> second_energy_level >> J_second;

            class AtomicLine atomic_line;
            atomic_line.wavelength = wavelength*1.0e-7; // Kurucz wavelengths are in nanometers
            const int g = int(2.0*J_first + 1.0);
            atomic_line.oscillator_strength = std::pow(10.0, log_gf) / double(g);

            // Now search through the model atom energy levels to find the lower and upper levels for this line.
            for (auto &level: levels) {
                if (level.energy < std::numeric_limits<double>::epsilon()) {  // if we're comparing to the ground state (which has energy 0), don't divide by it
                    if (std::abs(first_energy_level*h_planck*c_light - level.energy) < tolerance)
                        atomic_line.lower_level = &level;
                } else if ((std::abs(first_energy_level*h_planck*c_light - level.energy) / level.energy) < tolerance) {
                    atomic_line.lower_level = &level;
                }
            }
            for (auto &level: levels) {
                if (level.energy < std::numeric_limits<double>::epsilon()) {  // if we're comparing to the ground state (which has energy 0), don't divide by it
                    if (std::abs(second_energy_level*h_planck*c_light - level.energy) < tolerance)
                        atomic_line.upper_level = &level;
                } else if ((std::abs(second_energy_level*h_planck*c_light - level.energy) / level.energy) < tolerance) {
                    atomic_line.upper_level = &level;
                }
            }

            duplicate = false;
            // We assume the line lists don't have duplicate lines, but it can't hurt to check anyway.
            for (auto line: lines) {
                if ((std::abs(atomic_line.wavelength - line.wavelength) / line.wavelength) > tolerance) {
                    duplicate = false;
                } else {
                    duplicate = true;
                    break;
                }
            }

            // TODO: allow user to choose which line profile to use. For now the default is Gaussian.
            atomic_line.line_profile = gauss_profile;

            if (!duplicate) {
                lines.push_back(atomic_line);
            }
        }
        atomic_data_file.close();
    }
}


void Ion::calc_partition_function(const double temperature) {
    const double beta = 1.0 / (k_boltzmann * temperature);
    double result = 0.0;
    for (auto level: levels) {
        result += level.g * std::exp(-beta * level.energy);
    }
    partition_function = result;
}

std::ostream& operator<<(std::ostream& os, const Ion& ion) {
    os << std::setw(15) << "ION DATA: ";
    os << std::setw(20) << "atomic symbol" << std::setw(20) << "atomic number" << std::setw(20) << "ionization stage" << std::setw(20) << "atomic weight" << std::endl;
    os << std::setw(35) << ion.atomic_symbol << std::setw(20) << ion.atomic_number << std::setw(20) << ion.ionization_stage << std::setw(20) << ion.atomic_weight << std::endl;
    return os;
}


double AtomicLine::Einstein_B() {
    return (4.0 * std::pow(pi, 2) * std::pow(e_charge, 2) * wavelength) / (h_planck * m_electron * std::pow(c_light, 2)) * oscillator_strength;
}


double AtomicLine::Einstein_A() {
    return ((2.0 * h_planck * c_light) / std::pow(wavelength, 3)) * Einstein_B();
}


double AtomicLine::alpha(const double lambda) {
    return ((h_planck * c_light) / (4.0 * pi)) * (wavelength / c_light) * Einstein_B() * line_profile(lambda, wavelength, Delta_lambda);
}


void AtomicLine::set_line_width(const double temperature) {
    // TODO: let user choose among Doppler broadening and ... other types of broadening.
    Delta_lambda = Doppler_width(wavelength, temperature);
}


double AtomicLine::radiative_rate_absorption(const std::vector<GridWavelengthPoint> wavelength_grid) {
    double result = 0.0;
    for (auto it = wavelength_grid.begin(); it != wavelength_grid.end()-1; ++it) {
        auto it_next = it;
        std::advance(it_next, 1); // use std::next when more C++ compilers are C++11-compliant
        result += 0.5 * (*(it_next->lambda) - *(it->lambda)) * (alpha(*(it->lambda)) * it->J + alpha(*(it_next->lambda)) * it_next->J);
    }
    result *= (4.0 * pi) / (h_planck * c_light);
    return result;
}


double AtomicLine::radiative_rate_emission(const std::vector<GridWavelengthPoint> wavelength_grid, const double temperature) {
    double result = 0.0;
    for (auto it = wavelength_grid.begin(); it != wavelength_grid.end()-1; ++it) {
        auto it_next = it;
        std::advance(it_next, 1); // use std::next when more C++ compilers are C++11-compliant
        double f1, f2;
        f1 = alpha(*(it->lambda)) * (((2.0 * h_planck * std::pow(c_light, 2)) / std::pow(*(it->lambda), 5)) + it->J) * std::exp(-(h_planck * c_light) / (k_boltzmann **(it->lambda) * temperature));
        f2 = alpha(*(it_next->lambda)) * (((2.0 * h_planck * std::pow(c_light, 2)) / std::pow(*(it_next->lambda), 5)) + it_next->J) * std::exp(-(h_planck * c_light) / (k_boltzmann **(it_next->lambda) * temperature));
        result += 0.5 * (*(it_next->lambda) - *(it->lambda)) * (f1 + f2);
    }
    result *= (4.0 * pi) / (h_planck * c_light);
    return result;
}


// TODO: allow user to choose other rates besides this (the Van Regemorter (1962) rate)
double AtomicLine::collisional_rate_absorption(const double n_e, const double temperature) {
    double result = 0.0;
    const double y = (upper_level->energy - lower_level->energy) / (k_boltzmann * temperature); // defined in van Regemorter (1962)
    const double P_y = (3.0 * std::sqrt(3)) / (2.0 * pi) * (1.0 - y * std::exp(y) * gsl_sf_expint_E1(y));; // the integral in Eq 18 of van Regemorter (1962)

    // TODO: implement van Regemorter's more sophisticated values of P(y)
    result = 20.6 * std::pow(wavelength, 3) * n_e * std::pow(temperature, -0.5) * Einstein_A() * P_y;

    return result;
}


double Ion::photo_xs(const double lambda) const {

    if (ionization_stage == atomic_number) {
        std::cerr << "ERROR: A fully ionized atom has no photoionization cross-section ..." << std::endl;
        exit(1);
    }

    // Convert photon from wavelength in cm to energy in eV
    const double E = (h_planck * c_light / lambda) / ev2erg;

    // Fitting parameters from Verner et al.
    double sigma_0;
    double E_0;
    double y_w;
    double y_a;
    double P;
    double y_0;
    double y_1;

    switch (atomic_number) {
    case 1:
        switch (ionization_stage) {
        case 0:
            sigma_0 = 5.475e+4;
            E_0     = 4.298e-1;
            y_w     = 0.0e+0;
            y_a     = 3.288e+1;
            P       = 2.963e+0;
            y_0     = 0.0e+0;
            y_1     = 0.0e+0;

        }
    }

    const double x = (E / E_0) - y_0;
    const double y = std::sqrt(std::pow(x, 2) + std::pow(y_1, 2));

    const double F_y = (std::pow(x - 1.0, 2) + std::pow(y_w, 2)) * std::pow(y, 0.5 * P - 5.5) * std::pow(1.0 + std::sqrt(y / y_a), -P);

    return (sigma_0 * F_y * 1.0e-18);
}
