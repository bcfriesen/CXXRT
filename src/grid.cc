#include <algorithm>
#include <cmath>
#include <limits>

#include "globals.hh"
#include "grid.hh"
#include "constants.hh"
#include "EOS/Phi.hh"
#include "EOS/LTE_EOS.hh"
#include "misc/planck_function.hh"
#include "misc/planck_function_temperature_derivative.hh"

void GridVoxel::calc_J(const std::size_t wl_value_hash) {
    double result = 0.0;

    for (std::vector<RayIntersectionData>::const_iterator it = ray_intersection_data.begin(); it != ray_intersection_data.end()-1; ++it) {
        const std::vector<RayData>::iterator real_it = it->ray->raydata.begin() + it->intersection_point;

        std::vector<RayIntersectionData>::const_iterator it_next = it;
        std::advance(it_next, 1); // use std::next when more C++ compilers are C++11-compliant
        const std::vector<RayData>::iterator real_it_next = it_next->ray->raydata.begin() + it_next->intersection_point;

        result += 0.5 * (real_it->wavelength_grid[wl_value_hash].I + real_it_next->wavelength_grid[wl_value_hash].I) * (real_it_next->mu - real_it->mu);
    }
    wavelength_grid[wl_value_hash].J = 0.5 * result;
}


void GridVoxel::calc_H(const std::size_t wl_value_hash) {
    double result = 0.0;

    for (std::vector<RayIntersectionData>::const_iterator it = ray_intersection_data.begin(); it != ray_intersection_data.end()-1; ++it) {
        const std::vector<RayData>::iterator real_it = it->ray->raydata.begin() + it->intersection_point;

        std::vector<RayIntersectionData>::const_iterator it_next = it;
        std::advance(it_next, 1); // use std::next when more C++ compilers are C++11-compliant
        const std::vector<RayData>::iterator real_it_next = it_next->ray->raydata.begin() + it_next->intersection_point;

        result += 0.5 * (real_it->mu * real_it->wavelength_grid[wl_value_hash].I + real_it_next->mu * real_it_next->wavelength_grid[wl_value_hash].I) * (real_it_next->mu - real_it->mu);
    }
    wavelength_grid[wl_value_hash].H = 0.5 * result;
}


void GridVoxel::calc_K(const std::size_t wl_value_hash) {
    double result = 0.0;

    for (std::vector<RayIntersectionData>::const_iterator it = ray_intersection_data.begin(); it != ray_intersection_data.end()-1; ++it) {
        const std::vector<RayData>::iterator real_it = it->ray->raydata.begin() + it->intersection_point;

        std::vector<RayIntersectionData>::const_iterator it_next = it;
        std::advance(it_next, 1); // use std::next when more C++ compilers are C++11-compliant
        const std::vector<RayData>::iterator real_it_next = it_next->ray->raydata.begin() + it_next->intersection_point;

        result += 0.5 * (std::pow(real_it->mu, 2) * real_it->wavelength_grid[wl_value_hash].I + std::pow(real_it_next->mu, 2) * real_it_next->wavelength_grid[wl_value_hash].I) * (real_it_next->mu - real_it->mu);
    }
    wavelength_grid[wl_value_hash].K = 0.5 * result;
}


void GridVoxel::calc_LTE_populations() {
    for (std::vector<Atom>::iterator atom = atoms.begin(); atom != atoms.end(); ++atom) {
        for (std::vector<Ion>::iterator ion = atom->ions.begin(); ion != atom->ions.end(); ++ion) {
            for (std::vector<AtomicLevel>::iterator level = ion->levels.begin(); level != ion->levels.end(); ++level) {
                if (ion->ionization_stage <= atom->max_ionization_stage) {
                    level->number_density = atom->number_fraction * (n_g - n_e) * n_e * f_ij(*atom, *(ion->next_ion), n_e, temperature) * Phi_tilde(*level, *ion, *atom, temperature);
                } else {
                    // Treat fully ionized atom specially.
                    level->number_density = atom->number_fraction * (n_g - n_e) * f_ij(*atom, *ion, n_e, temperature);
                }
            }
        }
    }
}


void GridVoxel::calculate_emissivity_and_opacity(const std::size_t wl_value_hash) {
    double eta_tot = 0.0;
    double kappa_tot = 0.0;
    for (std::vector<Atom>::const_iterator atom = atoms.begin(); atom != atoms.end(); ++atom) {
        for (std::vector<Ion>::const_iterator ion = atom->ions.begin(); ion != atom->ions.end()-1; ++ion) {
            // Add up contributions from lines.
            for (std::vector<AtomicLine>::const_iterator line = ion->lines.begin(); line != ion->lines.end(); ++line) {
                eta_tot += line->eta(wavelength_values[wl_value_hash]);
                kappa_tot += line->kappa(wavelength_values[wl_value_hash]);
            }
            if (ion->ionization_stage <= atom->max_ionization_stage) {
                // Add up contributions from continua.
                for (std::vector<AtomicLevel>::const_iterator level = ion->levels.begin(); level != ion->levels.end(); ++level) {
                    eta_tot += ion->eta(wavelength_values[wl_value_hash], *level, n_e, temperature);
                    kappa_tot += ion->kappa(wavelength_values[wl_value_hash], *level, n_e, temperature);
                }
            }
        }
    }

    wavelength_grid[wl_value_hash].eta = eta_tot;
    wavelength_grid[wl_value_hash].kappa = kappa_tot;
    sigma = n_e * sigma_T;
    wavelength_grid[wl_value_hash].chi = wavelength_grid[wl_value_hash].kappa * (1.0 - std::exp((-h_planck * c_light) / (wavelength_values[wl_value_hash] * k_boltzmann * temperature))) + sigma;
    wavelength_grid[wl_value_hash].epsilon = wavelength_grid[wl_value_hash].kappa / (wavelength_grid[wl_value_hash].kappa + sigma);
}


void GridVoxel::calc_J_wl_integral() {
    J_wl_integral = 0.0;
    bool first_time = true; // I can't iterate over part of a map like I can a vector (i.e., "it != wavelength_grid.end()-1") so I have to use this bool hack.
    for (std::map<std::size_t, GridWavelengthPoint>::const_iterator it = wavelength_grid.begin(); it != wavelength_grid.end(); ++it) {
        if (!first_time) {
            std::map<std::size_t, GridWavelengthPoint>::const_iterator it_prev = it;
            std::advance(it_prev, -1);
            J_wl_integral += 0.5 * (it->second.lambda - it_prev->second.lambda) * (it->second.J + it_prev->second.J);
        }
        first_time = false;
    }
}


void GridVoxel::calc_H_wl_integral() {
    H_wl_integral = 0.0;
    bool first_time = true; // I can't iterate over part of a map like I can a vector (i.e., "it != wavelength_grid.end()-1") so I have to use this bool hack.
    for (std::map<std::size_t, GridWavelengthPoint>::const_iterator it = wavelength_grid.begin(); it != wavelength_grid.end(); ++it) {
        if (!first_time) {
            std::map<std::size_t, GridWavelengthPoint>::const_iterator it_prev = it;
            std::advance(it_prev, -1);
            H_wl_integral += 0.5 * (it->second.lambda - it_prev->second.lambda) * (it->second.H + it_prev->second.H);
        }
        first_time = false;
    }
}


void GridVoxel::calc_K_wl_integral() {
    K_wl_integral = 0.0;
    bool first_time = true; // I can't iterate over part of a map like I can a vector (i.e., "it != wavelength_grid.end()-1") so I have to use this bool hack.
    for (std::map<std::size_t, GridWavelengthPoint>::const_iterator it = wavelength_grid.begin(); it != wavelength_grid.end(); ++it) {
        if (!first_time) {
            std::map<std::size_t, GridWavelengthPoint>::const_iterator it_prev = it;
            std::advance(it_prev, -1);
            K_wl_integral += 0.5 * (it->second.lambda - it_prev->second.lambda) * (it->second.K + it_prev->second.K);
        }
        first_time = false;
    }
}


void GridVoxel::calc_chi_H() {
    chi_H = 0.0;

    bool first_time = true; // I can't iterate over part of a map like I can a vector (i.e., "it != wavelength_grid.end()-1") so I have to use this bool hack.
    for (std::map<std::size_t, double>::const_iterator it = wavelength_values.begin(); it != wavelength_values.end(); ++it) {
        if (!first_time) {
            std::map<std::size_t, double>::const_iterator it_prev = it;
            std::advance(it_prev, -1);
            const double lambda = wavelength_grid[it->first].lambda;
            const double lambda_prev = wavelength_grid[it_prev->first].lambda;
            const double H = wavelength_grid[it->first].H;
            const double H_prev = wavelength_grid[it_prev->first].H;
            const double chi = wavelength_grid[it->first].chi;
            const double chi_prev = wavelength_grid[it_prev->first].chi;
            chi_H += 0.5 * (lambda - lambda_prev) * (chi * H + chi_prev * H_prev);
        }
        first_time = false;
    }
    // Make sure we call calc_H_wl_integral before this!
    chi_H /= H_wl_integral;
}


void GridVoxel::calc_kappa_J() {
    kappa_J = 0.0;

    bool first_time = true; // I can't iterate over part of a map like I can a vector (i.e., "it != wavelength_grid.end()-1") so I have to use this bool hack.
    for (std::map<std::size_t, double>::const_iterator it = wavelength_values.begin(); it != wavelength_values.end(); ++it) {
        if (!first_time) {
            std::map<std::size_t, double>::const_iterator it_prev = it;
            std::advance(it_prev, -1);
            const double lambda = wavelength_grid[it->first].lambda;
            const double lambda_prev = wavelength_grid[it_prev->first].lambda;
            const double J = wavelength_grid[it->first].J;
            const double J_prev = wavelength_grid[it_prev->first].J;
            const double kappa = wavelength_grid[it->first].kappa;
            const double kappa_prev = wavelength_grid[it_prev->first].kappa;
            kappa_J += 0.5 * (lambda - lambda_prev) * (kappa * J + kappa_prev * J_prev);
        }
        first_time = false;
    }
    // Make sure we call calc_J_wl_integral before this!
    kappa_J /= J_wl_integral;
}


void GridVoxel::calc_kappa_B() {
    kappa_B = 0.0;

    bool first_time = true; // I can't iterate over part of a map like I can a vector (i.e., "it != wavelength_grid.end()-1") so I have to use this bool hack.
    for (std::map<std::size_t, double>::const_iterator it = wavelength_values.begin(); it != wavelength_values.end(); ++it) {
        if (!first_time) {
            std::map<std::size_t, double>::const_iterator it_prev = it;
            std::advance(it_prev, -1);
            const double lambda = wavelength_grid[it->first].lambda;
            const double lambda_prev = wavelength_grid[it_prev->first].lambda;
            const double B = planck_function(wavelength_values[it->first], temperature);
            const double B_prev = planck_function(wavelength_values[it_prev->first], temperature);
            const double kappa = wavelength_grid[it->first].kappa;
            const double kappa_prev = wavelength_grid[it_prev->first].kappa;
            kappa_B += 0.5 * (lambda - lambda_prev) * (kappa * B + kappa_prev * B_prev);
        }
        first_time = false;
    }
    kappa_B /= (sigma_stefan * std::pow(temperature, 4) / pi);
}


void GridVoxel::calc_Eddington_factor_f() {
    // Make sure we call calc_K_wl_integral() and calc_J_wl_integral before this!
    Eddington_factor_f = K_wl_integral / J_wl_integral;
}


void GridVoxel::calc_eta_minus_chi_J_wl_integral() {
    eta_minus_chi_J_wl_integral = 0.0;

    bool first_time = true; // I can't iterate over part of a map like I can a vector (i.e., "it != wavelength_grid.end()-1") so I have to use this bool hack.
    for (std::map<std::size_t, double>::const_iterator it = wavelength_values.begin(); it != wavelength_values.end(); ++it) {
        if (!first_time) {
            std::map<std::size_t, double>::const_iterator it_prev = it;
            std::advance(it_prev, -1);
            const double lambda = wavelength_grid[it->first].lambda;
            const double lambda_prev = wavelength_grid[it_prev->first].lambda;
            const double J = wavelength_grid[it->first].J;
            const double J_prev = wavelength_grid[it_prev->first].J;
            const double eta = wavelength_grid[it->first].eta;
            const double eta_prev = wavelength_grid[it_prev->first].eta;
            const double chi = wavelength_grid[it->first].chi;
            const double chi_prev = wavelength_grid[it_prev->first].chi;
            eta_minus_chi_J_wl_integral += 0.5 * (lambda - lambda_prev) * ((eta - chi * J) + (eta_prev - chi_prev * J_prev));
        }
        first_time = false;
    }
}


void GridVoxel::calc_Rosseland_mean_opacity() {
    Rosseland_mean_opacity = 0.0;

    bool first_time = true; // I can't iterate over part of a map like I can a vector (i.e., "it != wavelength_grid.end()-1") so I have to use this bool hack.
    for (std::map<std::size_t, double>::const_iterator wlv_it = wavelength_values.begin(); wlv_it != wavelength_values.end(); ++wlv_it) {
        if (!first_time) {
            std::map<std::size_t, double>::const_iterator wlv_it_prev = wlv_it;
            std::advance(wlv_it_prev, -1);

            // Just so it's easier to see what we're doing.
            const double lambda = wlv_it->second;
            const double lambda_prev = wlv_it_prev->second;

            const double chi = wavelength_grid[wlv_it->first].chi;
            const double chi_prev = wavelength_grid[wlv_it_prev->first].chi;

            Rosseland_mean_opacity += 0.5 * ((1.0 / chi) * planck_function_temperature_derivative(lambda, temperature) + (1.0 / chi_prev) * planck_function_temperature_derivative(lambda_prev, temperature)) * (lambda - lambda_prev);
        }
        first_time = false;
    }
    Rosseland_mean_opacity *= pi / (a_rad * c_light * std::pow(temperature, 3));
    Rosseland_mean_opacity = 1.0 / Rosseland_mean_opacity;
}


void GridVoxel::calc_source_fn(const std::size_t wl_value_hash) {
    wavelength_grid[wl_value_hash].source_fn = wavelength_grid[wl_value_hash].epsilon * planck_function(wavelength_grid[wl_value_hash].lambda, temperature) + (1.0 - wavelength_grid[wl_value_hash].epsilon) * wavelength_grid[wl_value_hash].J;
}
