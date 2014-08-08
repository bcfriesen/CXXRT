#include <algorithm>
#include <cmath>
#include <limits>

#include "grid.hh"
#include "constants.hh"
#include "EOS/Phi.hh"
#include "EOS/LTE_EOS.hh"

bool ray_angle_sort_function(const struct RayIntersectionData rd1, const struct RayIntersectionData rd2) {
    const auto it1 = rd1.ray->raydata.begin() + rd1.intersection_point;
    const auto it2 = rd2.ray->raydata.begin() + rd2.intersection_point;
    return (it1->mu < it2->mu);
}


void GridVoxel::calc_J(const double lambda) {
    std::sort(ray_intersection_data.begin(), ray_intersection_data.end(), ray_angle_sort_function);
    double result = 0.0;

    // Find the requested wavelength point on the grid voxel.
    // TODO: make this faster than a crude linear search.
    std::vector<GridWavelengthPoint>::iterator grid_wlp;
    for (grid_wlp = wavelength_grid.begin(); grid_wlp != wavelength_grid.end(); ++grid_wlp) {
        if (std::abs(*(grid_wlp->lambda) - lambda) < std::numeric_limits<double>::epsilon())
            break;
    }

    for (auto it = ray_intersection_data.begin(); it != ray_intersection_data.end()-1; ++it) {
        const auto real_it = it->ray->raydata.begin() + it->intersection_point;

        // Find the requested wavelength point on the rays.
        std::vector<RayWavelengthPoint>::const_iterator wlp_it;
        for (wlp_it = real_it->wavelength_grid.begin(); wlp_it != real_it->wavelength_grid.end(); ++wlp_it) {
            if (std::abs(*(wlp_it->lambda) - lambda) < std::numeric_limits<double>::epsilon())
                break;
        }

        auto it_next = it;
        std::advance(it_next, 1); // use std::next when more C++ compilers are C++11-compliant
        const auto real_it_next = it_next->ray->raydata.begin() + it_next->intersection_point;

        std::vector<RayWavelengthPoint>::const_iterator wlp_it_next;
        for (wlp_it_next = real_it_next->wavelength_grid.begin(); wlp_it_next != real_it_next->wavelength_grid.end(); ++wlp_it_next) {
            if (std::abs(*(wlp_it_next->lambda) - lambda) < std::numeric_limits<double>::epsilon())
                break;
        }

        result += 0.5 * (wlp_it->I + wlp_it_next->I) * (real_it_next->mu - real_it->mu);
    }
    grid_wlp->J = 0.5 * result;
}


void GridVoxel::calc_H(const double lambda) {
    std::sort(ray_intersection_data.begin(), ray_intersection_data.end(), ray_angle_sort_function);
    double result = 0.0;

    // Find the requested wavelength point on the grid voxel.
    // TODO: make this faster than a crude linear search.
    std::vector<GridWavelengthPoint>::iterator grid_wlp;
    for (grid_wlp = wavelength_grid.begin(); grid_wlp != wavelength_grid.end(); ++grid_wlp) {
        if (std::abs(*(grid_wlp->lambda) - lambda) < std::numeric_limits<double>::epsilon())
            break;
    }

    for (auto it = ray_intersection_data.begin(); it != ray_intersection_data.end()-1; ++it) {
        const auto real_it = it->ray->raydata.begin() + it->intersection_point;

        // Find the requested wavelength point on the rays.
        std::vector<RayWavelengthPoint>::const_iterator wlp_it;
        for (wlp_it = real_it->wavelength_grid.begin(); wlp_it != real_it->wavelength_grid.end(); ++wlp_it) {
            if (std::abs(*(wlp_it->lambda) - lambda) < std::numeric_limits<double>::epsilon())
                break;
        }

        auto it_next = it;
        std::advance(it_next, 1); // use std::next when more C++ compilers are C++11-compliant
        const auto real_it_next = it_next->ray->raydata.begin() + it_next->intersection_point;

        std::vector<RayWavelengthPoint>::const_iterator wlp_it_next;
        for (wlp_it_next = real_it_next->wavelength_grid.begin(); wlp_it_next != real_it_next->wavelength_grid.end(); ++wlp_it_next) {
            if (std::abs(*(wlp_it_next->lambda) - lambda) < std::numeric_limits<double>::epsilon())
                break;
        }

        result += 0.5 * (real_it->mu * wlp_it->I + real_it_next->mu * wlp_it_next->I) * (real_it_next->mu - real_it->mu);
    }
    grid_wlp->H = 0.5 * result;
}


void GridVoxel::calc_K(const double lambda) {
    std::sort(ray_intersection_data.begin(), ray_intersection_data.end(), ray_angle_sort_function);
    double result = 0.0;

    // Find the requested wavelength point on the grid voxel.
    // TODO: make this faster than a crude linear search.
    std::vector<GridWavelengthPoint>::iterator grid_wlp;
    for (grid_wlp = wavelength_grid.begin(); grid_wlp != wavelength_grid.end(); ++grid_wlp) {
        if (std::abs(*(grid_wlp->lambda) - lambda) < std::numeric_limits<double>::epsilon())
            break;
    }

    for (auto it = ray_intersection_data.begin(); it != ray_intersection_data.end()-1; ++it) {
        const auto real_it = it->ray->raydata.begin() + it->intersection_point;

        // Find the requested wavelength point on the rays.
        std::vector<RayWavelengthPoint>::const_iterator wlp_it;
        for (wlp_it = real_it->wavelength_grid.begin(); wlp_it != real_it->wavelength_grid.end(); ++wlp_it) {
            if (std::abs(*(wlp_it->lambda) - lambda) < std::numeric_limits<double>::epsilon())
                break;
        }

        auto it_next = it;
        std::advance(it_next, 1); // use std::next when more C++ compilers are C++11-compliant
        const auto real_it_next = it_next->ray->raydata.begin() + it_next->intersection_point;

        std::vector<RayWavelengthPoint>::const_iterator wlp_it_next;
        for (wlp_it_next = real_it_next->wavelength_grid.begin(); wlp_it_next != real_it_next->wavelength_grid.end(); ++wlp_it_next) {
            if (std::abs(*(wlp_it_next->lambda) - lambda) < std::numeric_limits<double>::epsilon())
                break;
        }

        result += 0.5 * (std::pow(real_it->mu, 2) * wlp_it->I + std::pow(real_it_next->mu, 2) * wlp_it_next->I) * (real_it_next->mu - real_it->mu);
    }
    grid_wlp->K = 0.5 * result;
}


void GridVoxel::calc_LTE_populations() {
    for (auto &atom: atoms) {
        for (auto &ion: atom.ions) {
            for (auto &level: ion.levels) {
                if (ion.ionization_stage < ion.atomic_number) {
                    level.number_density = atom.number_fraction * (n_g - n_e) * n_e * f_ij(atom, *(ion.next_ion), n_e, temperature) * Phi_tilde(level, ion, atom, temperature);
                } else {
                    // Treat fully ionized atom specially.
                    level.number_density = atom.number_fraction * (n_g - n_e) * f_ij(atom, ion, n_e, temperature);
                }
            }
        }
    }
}


void GridVoxel::calculate_emissivity_and_opacity(const double lambda) {
    double eta_tot = 0.0;
    double kappa_tot = 0.0;
    for (auto &atom: atoms) {
        for (auto &ion: atom.ions) {
            for (auto &line: ion.lines) {
                eta_tot += line.eta(lambda);
                kappa_tot += line.kappa(lambda);
            }
        }
    }

    // Since we don't treat anisotropic sources or sinks, all intersecting rays at this voxel get the same value for eta and kappa.
    for (auto rid_it = ray_intersection_data.begin(); rid_it != ray_intersection_data.end(); ++rid_it) {
        const auto ray_it = rid_it->ray->raydata.begin() + rid_it->intersection_point;

        // Find the requested wavelength point on the rays.
        // TODO: make this faster than a crude linear search.
        std::vector<RayWavelengthPoint>::iterator wlp_it;
        for (wlp_it = ray_it->wavelength_grid.begin(); wlp_it != ray_it->wavelength_grid.end(); ++wlp_it) {
            if (std::abs(*(wlp_it->lambda) - lambda) < std::numeric_limits<double>::epsilon())
                break;
        }
        wlp_it->eta = eta_tot;
        wlp_it->kappa = kappa_tot;
        wlp_it->sigma = n_e * sigma_T;
        wlp_it->chi = wlp_it->kappa + wlp_it->sigma;

        // Find the requested wavelength point on the grid voxel.
        // TODO: make this faster than a crude linear search.
        std::vector<GridWavelengthPoint>::iterator grid_wlp;
        for (grid_wlp = wavelength_grid.begin(); grid_wlp != wavelength_grid.end(); ++grid_wlp) {
            if (std::abs(*(grid_wlp->lambda) - lambda) < std::numeric_limits<double>::epsilon())
                break;
        }
        // In general the thermalization parameter is a ray-dependent quantity
        // since it depends on kappa (which is ray-dependent). However in the
        // equivalent-two-level-atom formalism it is ray-independent, so just
        // set the grid scalar value to the ray value.
        grid_wlp->epsilon = wlp_it->kappa / (wlp_it->kappa + wlp_it->sigma);
    }
}
