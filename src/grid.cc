#include <algorithm>
#include <cmath>
#include <limits>

#include "globals.hh"
#include "grid.hh"
#include "constants.hh"
#include "EOS/Phi.hh"
#include "EOS/LTE_EOS.hh"
#include "planck_function.hh"

void GridVoxel::calc_J(const std::size_t wl_value_hash) {
    double result = 0.0;

    for (auto it = ray_intersection_data.begin(); it != ray_intersection_data.end()-1; ++it) {
        const auto real_it = it->ray->raydata.begin() + it->intersection_point;

        auto it_next = it;
        std::advance(it_next, 1); // use std::next when more C++ compilers are C++11-compliant
        const auto real_it_next = it_next->ray->raydata.begin() + it_next->intersection_point;

        result += 0.5 * (real_it->wavelength_grid[wl_value_hash].I + real_it_next->wavelength_grid[wl_value_hash].I) * (real_it_next->mu - real_it->mu);
    }
    wavelength_grid[wl_value_hash].J = 0.5 * result;
}


void GridVoxel::calc_H(const std::size_t wl_value_hash) {
    double result = 0.0;

    for (auto it = ray_intersection_data.begin(); it != ray_intersection_data.end()-1; ++it) {
        const auto real_it = it->ray->raydata.begin() + it->intersection_point;

        auto it_next = it;
        std::advance(it_next, 1); // use std::next when more C++ compilers are C++11-compliant
        const auto real_it_next = it_next->ray->raydata.begin() + it_next->intersection_point;

        result += 0.5 * (real_it->mu * real_it->wavelength_grid[wl_value_hash].I + real_it_next->mu * real_it_next->wavelength_grid[wl_value_hash].I) * (real_it_next->mu - real_it->mu);
    }
    wavelength_grid[wl_value_hash].H = 0.5 * result;
}


void GridVoxel::calc_K(const std::size_t wl_value_hash) {
    double result = 0.0;

    for (auto it = ray_intersection_data.begin(); it != ray_intersection_data.end()-1; ++it) {
        const auto real_it = it->ray->raydata.begin() + it->intersection_point;

        auto it_next = it;
        std::advance(it_next, 1); // use std::next when more C++ compilers are C++11-compliant
        const auto real_it_next = it_next->ray->raydata.begin() + it_next->intersection_point;

        result += 0.5 * (std::pow(real_it->mu, 2) * real_it->wavelength_grid[wl_value_hash].I + std::pow(real_it_next->mu, 2) * real_it_next->wavelength_grid[wl_value_hash].I) * (real_it_next->mu - real_it->mu);
    }
    wavelength_grid[wl_value_hash].K = 0.5 * result;
}


void GridVoxel::calc_LTE_populations() {
    for (auto atom = atoms.begin(); atom != atoms.end(); ++atom) {
        for (auto ion = atom->ions.begin(); ion != atom->ions.end(); ++ion) {
            for (auto level = ion->levels.begin(); level != ion->levels.end(); ++level) {
                if (ion->ionization_stage < ion->atomic_number) {
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
    for (auto atom = atoms.begin(); atom != atoms.end(); ++atom) {
        for (auto ion = atom->ions.begin(); ion != atom->ions.end(); ++ion) {
            // Add up contributions from lines.
            for (auto line = ion->lines.begin(); line != ion->lines.end(); ++line) {
                eta_tot += line->eta(wavelength_values[wl_value_hash]);
                kappa_tot += line->kappa(wavelength_values[wl_value_hash]);
            }
            if (ion->ionization_stage < ion->atomic_number) {
                // Add up contributions from continua.
                for (auto level: ion->levels) {
                    eta_tot += ion->eta(wavelength_values[wl_value_hash], level, n_e, temperature);
                    kappa_tot += ion->kappa(wavelength_values[wl_value_hash], level, n_e, temperature);
                }
            }
        }
    }

    // Since we don't treat anisotropic sources or sinks, all intersecting rays at this voxel get the same value for eta and kappa.
    for (auto rid_it = ray_intersection_data.begin(); rid_it != ray_intersection_data.end(); ++rid_it) {
        const auto ray_it = rid_it->ray->raydata.begin() + rid_it->intersection_point;

        ray_it->wavelength_grid[wl_value_hash].eta = eta_tot;
        ray_it->wavelength_grid[wl_value_hash].kappa = kappa_tot;
        ray_it->wavelength_grid[wl_value_hash].sigma = n_e * sigma_T;
        ray_it->wavelength_grid[wl_value_hash].chi = ray_it->wavelength_grid[wl_value_hash].kappa + ray_it->wavelength_grid[wl_value_hash].sigma;
        ray_it->wavelength_grid[wl_value_hash].epsilon = ray_it->wavelength_grid[wl_value_hash].kappa / (ray_it->wavelength_grid[wl_value_hash].kappa + ray_it->wavelength_grid[wl_value_hash].sigma);

        // In general the thermalization parameter is a ray-dependent quantity
        // since it depends on kappa (which is ray-dependent). However in the
        // equivalent-two-level-atom formalism it is ray-independent, so just
        // set the grid scalar value to the ray value.
        wavelength_grid[wl_value_hash].epsilon = ray_it->wavelength_grid[wl_value_hash].kappa / (ray_it->wavelength_grid[wl_value_hash].kappa + ray_it->wavelength_grid[wl_value_hash].sigma);
    }
}


double GridVoxel::calc_J_wl_integral() const {
    double result = 0.0;
    bool first_time = true; // I can't iterate over part of a map like I can a vector (i.e., "it != wavelength_grid.end()-1") so I have to use this bool hack.
    for (auto it = wavelength_grid.begin(); it != wavelength_grid.end(); ++it) {
        if (!first_time) {
            auto it_prev = it;
            std::advance(it_prev, -1);
            result += 0.5 * (it->second.lambda - it_prev->second.lambda) * (it->second.J - it_prev->second.J);
        }
        first_time = false;
    }

    return result;
}


double GridVoxel::calc_H_wl_integral() const {
    double result = 0.0;
    bool first_time = true; // I can't iterate over part of a map like I can a vector (i.e., "it != wavelength_grid.end()-1") so I have to use this bool hack.
    for (auto it = wavelength_grid.begin(); it != wavelength_grid.end(); ++it) {
        if (!first_time) {
            auto it_prev = it;
            std::advance(it_prev, -1);
            result += 0.5 * (it->second.lambda - it_prev->second.lambda) * (it->second.H - it_prev->second.H);
        }
        first_time = false;
    }

    return result;
}


double GridVoxel::calc_K_wl_integral() const {
    double result = 0.0;
    bool first_time = true; // I can't iterate over part of a map like I can a vector (i.e., "it != wavelength_grid.end()-1") so I have to use this bool hack.
    for (auto it = wavelength_grid.begin(); it != wavelength_grid.end(); ++it) {
        if (!first_time) {
            auto it_prev = it;
            std::advance(it_prev, -1);
            result += 0.5 * (it->second.lambda - it_prev->second.lambda) * (it->second.K - it_prev->second.K);
        }
        first_time = false;
    }

    return result;
}


double GridVoxel::calc_chi_H() {
    double result = 0.0;

    // chi is in principle a ray-dependent quantity, but all rays should have
    // the same value of chi in a given voxel. So just grab it from the first
    // ray we can find.
    const auto first_ray = ray_intersection_data.begin();
    const auto first_ray_raydata = first_ray->ray->raydata.begin() + first_ray->intersection_point;

    bool first_time = true; // I can't iterate over part of a map like I can a vector (i.e., "it != wavelength_grid.end()-1") so I have to use this bool hack.
    for (auto it = wavelength_values.begin(); it != wavelength_values.end(); ++it) {
        if (!first_time) {
            auto it_prev = it;
            std::advance(it_prev, -1);
            const double lambda = wavelength_grid[it->first].lambda;
            const double lambda_prev = wavelength_grid[it_prev->first].lambda;
            const double H = wavelength_grid[it->first].H;
            const double H_prev = wavelength_grid[it_prev->first].H;
            const double chi = first_ray_raydata->wavelength_grid[it->first].chi;
            const double chi_prev = first_ray_raydata->wavelength_grid[it_prev->first].chi;
            result += 0.5 * (lambda - lambda_prev) * (chi * H - chi_prev * H_prev);
        }
        first_time = false;
    }
    result /= calc_H_wl_integral();

    return result;
}


double GridVoxel::calc_kappa_J() {
    double result = 0.0;

    // kappa is in principle a ray-dependent quantity, but all rays should have
    // the same value of kappa in a given voxel. So just grab it from the first
    // ray we can find.
    const auto first_ray = ray_intersection_data.begin();
    const auto first_ray_raydata = first_ray->ray->raydata.begin() + first_ray->intersection_point;

    bool first_time = true; // I can't iterate over part of a map like I can a vector (i.e., "it != wavelength_grid.end()-1") so I have to use this bool hack.
    for (auto it = wavelength_values.begin(); it != wavelength_values.end(); ++it) {
        if (!first_time) {
            auto it_prev = it;
            std::advance(it_prev, -1);
            const double lambda = wavelength_grid[it->first].lambda;
            const double lambda_prev = wavelength_grid[it_prev->first].lambda;
            const double J = wavelength_grid[it->first].J;
            const double J_prev = wavelength_grid[it_prev->first].J;
            const double kappa = first_ray_raydata->wavelength_grid[it->first].kappa;
            const double kappa_prev = first_ray_raydata->wavelength_grid[it_prev->first].kappa;
            result += 0.5 * (lambda - lambda_prev) * (kappa * J - kappa_prev * J_prev);
        }
        first_time = false;
    }
    result /= calc_J_wl_integral();

    return result;
}


double GridVoxel::calc_kappa_B() {
    double result = 0.0;

    // kappa is in principle a ray-dependent quantity, but all rays should have
    // the same value of kappa in a given voxel. So just grab it from the first
    // ray we can find.
    const auto first_ray = ray_intersection_data.begin();
    const auto first_ray_raydata = first_ray->ray->raydata.begin() + first_ray->intersection_point;

    bool first_time = true; // I can't iterate over part of a map like I can a vector (i.e., "it != wavelength_grid.end()-1") so I have to use this bool hack.
    for (auto it = wavelength_values.begin(); it != wavelength_values.end(); ++it) {
        if (!first_time) {
            auto it_prev = it;
            std::advance(it_prev, -1);
            const double lambda = wavelength_grid[it->first].lambda;
            const double lambda_prev = wavelength_grid[it_prev->first].lambda;
            const double B = planck_function(wavelength_values[it->first], temperature);
            const double B_prev = planck_function(wavelength_values[it_prev->first], temperature);
            const double kappa = first_ray_raydata->wavelength_grid[it->first].kappa;
            const double kappa_prev = first_ray_raydata->wavelength_grid[it_prev->first].kappa;
            result += 0.5 * (lambda - lambda_prev) * (kappa * B - kappa_prev * B_prev);
        }
        first_time = false;
    }
    result /= (sigma_stefan * std::pow(temperature, 4) / pi);

    return result;
}


double GridVoxel::calc_Eddington_factor_f() const {
    return calc_K_wl_integral() / calc_J_wl_integral();
}


double GridVoxel::calc_eta_minus_chi_J_wl_integral() {
    double result = 0.0;

    // chi and eta are in principle ray-dependent quantities, but all rays
    // should have the same values of ch and eta in a given voxel. So just grab
    // them from the first ray we can find.
    const auto first_ray = ray_intersection_data.begin();
    const auto first_ray_raydata = first_ray->ray->raydata.begin() + first_ray->intersection_point;

    bool first_time = true; // I can't iterate over part of a map like I can a vector (i.e., "it != wavelength_grid.end()-1") so I have to use this bool hack.
    for (auto it = wavelength_values.begin(); it != wavelength_values.end(); ++it) {
        if (!first_time) {
            auto it_prev = it;
            std::advance(it_prev, -1);
            const double lambda = wavelength_grid[it->first].lambda;
            const double lambda_prev = wavelength_grid[it_prev->first].lambda;
            const double J = wavelength_grid[it->first].J;
            const double J_prev = wavelength_grid[it_prev->first].J;
            const double eta = first_ray_raydata->wavelength_grid[it->first].eta;
            const double eta_prev = first_ray_raydata->wavelength_grid[it_prev->first].eta;
            const double chi = first_ray_raydata->wavelength_grid[it->first].chi;
            const double chi_prev = first_ray_raydata->wavelength_grid[it_prev->first].chi;
            result += 0.5 * (lambda - lambda_prev) * ((eta - chi * J) - (eta_prev - chi_prev * J_prev));
        }
        first_time = false;
    }
    return result;
}
