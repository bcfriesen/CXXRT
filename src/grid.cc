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

void GridVoxel::calc_J(const unsigned int wl_index) {
    double result = 0.0;

    for (std::vector<RayIntersectionData>::const_iterator it = ray_intersection_data.begin(); it != ray_intersection_data.end()-1; ++it) {
        const std::vector<RayData>::iterator real_it = it->ray->raydata.begin() + it->intersection_point;

        std::vector<RayIntersectionData>::const_iterator it_next = it;
        std::advance(it_next, 1); // use std::next when more C++ compilers are C++11-compliant
        const std::vector<RayData>::iterator real_it_next = it_next->ray->raydata.begin() + it_next->intersection_point;

        result += 0.5 * (real_it->wavelength_grid.at(wl_index).I + real_it_next->wavelength_grid.at(wl_index).I) * (real_it_next->mu - real_it->mu);
    }
    wavelength_grid.at(wl_index).J = 0.5 * result;
}


void GridVoxel::calc_H(const unsigned int wl_index) {
    double result = 0.0;

    for (std::vector<RayIntersectionData>::const_iterator it = ray_intersection_data.begin(); it != ray_intersection_data.end()-1; ++it) {
        const std::vector<RayData>::iterator real_it = it->ray->raydata.begin() + it->intersection_point;

        std::vector<RayIntersectionData>::const_iterator it_next = it;
        std::advance(it_next, 1); // use std::next when more C++ compilers are C++11-compliant
        const std::vector<RayData>::iterator real_it_next = it_next->ray->raydata.begin() + it_next->intersection_point;

        result += 0.5 * (real_it->mu * real_it->wavelength_grid.at(wl_index).I + real_it_next->mu * real_it_next->wavelength_grid.at(wl_index).I) * (real_it_next->mu - real_it->mu);
    }
    wavelength_grid.at(wl_index).H = 0.5 * result;
}


void GridVoxel::calc_K(const unsigned int wl_index) {
    double result = 0.0;

    for (std::vector<RayIntersectionData>::const_iterator it = ray_intersection_data.begin(); it != ray_intersection_data.end()-1; ++it) {
        const std::vector<RayData>::iterator real_it = it->ray->raydata.begin() + it->intersection_point;

        std::vector<RayIntersectionData>::const_iterator it_next = it;
        std::advance(it_next, 1); // use std::next when more C++ compilers are C++11-compliant
        const std::vector<RayData>::iterator real_it_next = it_next->ray->raydata.begin() + it_next->intersection_point;

        result += 0.5 * (std::pow(real_it->mu, 2) * real_it->wavelength_grid.at(wl_index).I + std::pow(real_it_next->mu, 2) * real_it_next->wavelength_grid.at(wl_index).I) * (real_it_next->mu - real_it->mu);
    }
    wavelength_grid.at(wl_index).K = 0.5 * result;
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


void GridVoxel::calculate_emissivity_and_opacity(const unsigned int wl_index) {
    double eta_tot = 0.0;
    double kappa_tot = 0.0;
    for (std::vector<Atom>::const_iterator atom = atoms.begin(); atom != atoms.end(); ++atom) {
        for (std::vector<Ion>::const_iterator ion = atom->ions.begin(); ion != atom->ions.end()-1; ++ion) {
            // Add up contributions from lines.
            for (std::vector<AtomicLine>::const_iterator line = ion->lines.begin(); line != ion->lines.end(); ++line) {
                eta_tot += line->eta(wavelength_values.at(wl_index));
                kappa_tot += line->kappa(wavelength_values.at(wl_index));
            }
            if (ion->ionization_stage <= atom->max_ionization_stage) {
                // Add up contributions from continua.
                for (std::vector<AtomicLevel>::const_iterator level = ion->levels.begin(); level != ion->levels.end(); ++level) {
                    eta_tot += ion->eta(wavelength_values.at(wl_index), *level, n_e, temperature);
                    kappa_tot += ion->kappa(wavelength_values.at(wl_index), *level, n_e, temperature);
                }
            }
        }
    }

    wavelength_grid.at(wl_index).eta = eta_tot;
    wavelength_grid.at(wl_index).kappa = kappa_tot;
    sigma = n_e * sigma_T;
    wavelength_grid.at(wl_index).chi = wavelength_grid.at(wl_index).kappa * (1.0 - std::exp((-h_planck * c_light) / (wavelength_values.at(wl_index) * k_boltzmann * temperature))) + sigma;
    wavelength_grid.at(wl_index).epsilon = wavelength_grid.at(wl_index).kappa / (wavelength_grid.at(wl_index).kappa + sigma);
}


void GridVoxel::calc_J_wl_integral() {
    J_wl_integral = 0.0;
    for (std::vector<GridWavelengthPoint>::const_iterator it = wavelength_grid.begin()+1; it != wavelength_grid.end(); ++it) {
        std::vector<GridWavelengthPoint>::const_iterator it_prev = it;
        std::advance(it_prev, -1);
        J_wl_integral += 0.5 * (it->lambda - it_prev->lambda) * (it->J + it_prev->J);
    }
}


void GridVoxel::calc_H_wl_integral() {
    H_wl_integral = 0.0;
    for (std::vector<GridWavelengthPoint>::const_iterator it = wavelength_grid.begin()+1; it != wavelength_grid.end(); ++it) {
        std::vector<GridWavelengthPoint>::const_iterator it_prev = it;
        std::advance(it_prev, -1);
        H_wl_integral += 0.5 * (it->lambda - it_prev->lambda) * (it->H + it_prev->H);
    }
}


void GridVoxel::calc_K_wl_integral() {
    K_wl_integral = 0.0;
    for (std::vector<GridWavelengthPoint>::const_iterator it = wavelength_grid.begin()+1; it != wavelength_grid.end(); ++it) {
        std::vector<GridWavelengthPoint>::const_iterator it_prev = it;
        std::advance(it_prev, -1);
        K_wl_integral += 0.5 * (it->lambda - it_prev->lambda) * (it->K + it_prev->K);
    }
}


void GridVoxel::calc_chi_H() {
    chi_H = 0.0;

    for (unsigned int i = 1; i < wavelength_values.size(); ++i) {
        const double lambda = wavelength_grid.at(i).lambda;
        const double lambda_prev = wavelength_grid.at(i-1).lambda;
        const double H = wavelength_grid.at(i).H;
        const double H_prev = wavelength_grid.at(i-1).H;
        const double chi = wavelength_grid.at(i).chi;
        const double chi_prev = wavelength_grid.at(i-1).chi;
        chi_H += 0.5 * (lambda - lambda_prev) * (chi * H + chi_prev * H_prev);
    }
    // Make sure we call calc_H_wl_integral before this!
    chi_H /= H_wl_integral;
}


void GridVoxel::calc_kappa_J() {
    kappa_J = 0.0;

    for (unsigned int i = 1; i < wavelength_values.size(); ++i) {
        const double lambda = wavelength_grid.at(i).lambda;
        const double lambda_prev = wavelength_grid.at(i-1).lambda;
        const double J = wavelength_grid.at(i).J;
        const double J_prev = wavelength_grid.at(i-1).J;
        const double kappa = wavelength_grid.at(i).kappa;
        const double kappa_prev = wavelength_grid.at(i-1).kappa;
        kappa_J += 0.5 * (lambda - lambda_prev) * (kappa * J + kappa_prev * J_prev);
    }
    // Make sure we call calc_J_wl_integral before this!
    kappa_J /= J_wl_integral;
}


void GridVoxel::calc_kappa_B() {
    kappa_B = 0.0;

    for (unsigned int i = 1; i < wavelength_values.size(); ++i) {
        const double lambda = wavelength_grid.at(i).lambda;
        const double lambda_prev = wavelength_grid.at(i-1).lambda;
        const double B = planck_function(wavelength_values.at(i), temperature);
        const double B_prev = planck_function(wavelength_values.at(i-1), temperature);
        const double kappa = wavelength_grid.at(i).kappa;
        const double kappa_prev = wavelength_grid.at(i-1).kappa;
        kappa_B += 0.5 * (lambda - lambda_prev) * (kappa * B + kappa_prev * B_prev);
    }
    kappa_B /= (sigma_stefan * std::pow(temperature, 4) / pi);
}


void GridVoxel::calc_Eddington_factor_f() {
    // Make sure we call calc_K_wl_integral() and calc_J_wl_integral before this!
    Eddington_factor_f = K_wl_integral / J_wl_integral;
}


void GridVoxel::calc_eta_minus_chi_J_wl_integral() {
    eta_minus_chi_J_wl_integral = 0.0;

    for (unsigned int i = 1; i < wavelength_values.size(); ++i) {
        const double lambda = wavelength_grid.at(i).lambda;
        const double lambda_prev = wavelength_grid.at(i-1).lambda;
        const double J = wavelength_grid.at(i).J;
        const double J_prev = wavelength_grid.at(i-1).J;
        const double eta = wavelength_grid.at(i).eta;
        const double eta_prev = wavelength_grid.at(i-1).eta;
        const double chi = wavelength_grid.at(i).chi;
        const double chi_prev = wavelength_grid.at(i-1).chi;
        eta_minus_chi_J_wl_integral += 0.5 * (lambda - lambda_prev) * ((eta - chi * J) + (eta_prev - chi_prev * J_prev));
    }
}


void GridVoxel::calc_Rosseland_mean_opacity() {
    Rosseland_mean_opacity = 0.0;

    for (unsigned int i = 1; i < wavelength_values.size(); ++i) {

        // Just so it's easier to see what we're doing.
        const double lambda = wavelength_values.at(i);
        const double lambda_prev = wavelength_values.at(i-1);

        const double chi = wavelength_grid.at(i).chi;
        const double chi_prev = wavelength_grid.at(i-1).chi;

        Rosseland_mean_opacity += 0.5 * ((1.0 / chi) * planck_function_temperature_derivative(lambda, temperature) + (1.0 / chi_prev) * planck_function_temperature_derivative(lambda_prev, temperature)) * (lambda - lambda_prev);
    }
    Rosseland_mean_opacity *= pi / (a_rad * c_light * std::pow(temperature, 3));
    Rosseland_mean_opacity = 1.0 / Rosseland_mean_opacity;
}


void GridVoxel::calc_source_fn(const unsigned int wl_index) {
    wavelength_grid.at(wl_index).source_fn = wavelength_grid.at(wl_index).epsilon * planck_function(wavelength_grid.at(wl_index).lambda, temperature) + (1.0 - wavelength_grid.at(wl_index).epsilon) * wavelength_grid.at(wl_index).J;
}
