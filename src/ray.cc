#include <algorithm>
#include <iomanip>
#include <iterator>
#include <cmath>
#include <fstream>

#include "planck_function.hh"
#include "ray.hh"
#include "globals.hh"
#include "grid.hh"

void Ray::bind_to_grid(const double mu) {
    const unsigned int n_depth_pts = grid.size();
    raydata.resize(grid.size());
    for (std::vector<RayData>::iterator rd = raydata.begin(); rd != raydata.end(); ++rd) {
        rd->mu = mu;
    }
    int i = 0;
    for (std::vector<RayData>::iterator it = raydata.begin(); it != raydata.end(); ++it) {
        struct RayIntersectionData ray_intersection_data;
        if (raydata.front().mu < 0.0) {
            it->gridvoxel = &grid[i];
            ray_intersection_data.intersection_point = it - raydata.begin();
        } else {
            it->gridvoxel = &grid[(n_depth_pts-1)-i];
            ray_intersection_data.intersection_point = (raydata.end()-1) - it;
        }
        ray_intersection_data.ray = this;
        grid.at(i).ray_intersection_data.push_back(ray_intersection_data);
        ++i;
    }
}

void Ray::print_ray_data(const std::size_t wl_value_hash) {
    log_file << std::endl;
    log_file << std::setw(15) << "# RAY DATA:";
    log_file << std::setw(15) << "z";
    log_file << std::setw(15) << "mu";
    log_file << std::setw(15) << "lambda";
    log_file << std::setw(15) << "I_lam";
    log_file << std::setw(15) << "sigma";
    log_file << std::setw(15) << "kappa";
    log_file << std::setw(15) << "eta";
    log_file << std::setw(15) << "tau";
    log_file << std::setw(15) << "chi";
    log_file << std::setw(15) << "source_fn";
    log_file << std::setw(15) << "Delta_tau";
    log_file << std::setw(15) << "alpha";
    log_file << std::setw(15) << "beta";
    log_file << std::setw(15) << "gamma";
    log_file << std::endl;

    for (std::vector<RayData>::iterator rd = raydata.begin(); rd != raydata.end(); ++rd) {
        log_file << std::setw(30) << std::scientific << rd->gridvoxel->z;
        log_file << std::setw(15) << std::scientific << rd->mu;
        log_file << std::setw(15) << std::scientific << rd->wavelength_grid[wl_value_hash].lambda * 1.0e+8;
        log_file << std::setw(15) << std::scientific << rd->wavelength_grid[wl_value_hash].I;
        log_file << std::setw(15) << std::scientific << rd->wavelength_grid[wl_value_hash].sigma;
        log_file << std::setw(15) << std::scientific << rd->wavelength_grid[wl_value_hash].kappa;
        log_file << std::setw(15) << std::scientific << rd->wavelength_grid[wl_value_hash].eta;
        log_file << std::setw(15) << std::scientific << rd->wavelength_grid[wl_value_hash].tau;
        log_file << std::setw(15) << std::scientific << rd->wavelength_grid[wl_value_hash].chi;
        log_file << std::setw(15) << std::scientific << rd->wavelength_grid[wl_value_hash].source_fn;
        log_file << std::setw(15) << std::scientific << rd->wavelength_grid[wl_value_hash].Delta_tau;
        log_file << std::setw(15) << std::scientific << rd->wavelength_grid[wl_value_hash].alpha;
        log_file << std::setw(15) << std::scientific << rd->wavelength_grid[wl_value_hash].beta;
        log_file << std::setw(15) << std::scientific << rd->wavelength_grid[wl_value_hash].gamma;
        log_file << std::endl;
    }
    log_file << std::endl;
}


RayData::RayData() {
    gridvoxel = NULL;

    for (std::map<std::size_t, double>::const_iterator wlv = wavelength_values.begin(); wlv != wavelength_values.end(); ++wlv) {
        RayWavelengthPoint wlp_tmp;
        wlp_tmp.lambda = wlv->second;
        wavelength_grid[wlv->first] = wlp_tmp;
    }

    for (std::map<std::size_t, RayWavelengthPoint>::iterator wlp = wavelength_grid.begin(); wlp != wavelength_grid.end(); ++wlp) {
        wlp->second.I = 0.0;
        wlp->second.tau = 0.0;
        wlp->second.chi = 0.0;
        wlp->second.source_fn = 0.0;
    }
}

void Ray::calc_tau(const std::size_t wl_value_hash) {
    for (std::vector<RayData>::iterator it = raydata.begin(); it != raydata.end(); ++it) {
        if (it == raydata.begin()) {
            it->wavelength_grid[wl_value_hash].tau = 0.0;
        } else {
            std::vector<RayData>::iterator it_prev = it;
            std::advance(it_prev, -1); // use std::prev when more C++ compilers are C++11-compliant
            it->wavelength_grid[wl_value_hash].tau = it_prev->wavelength_grid[wl_value_hash].tau + (0.5 * (it_prev->wavelength_grid[wl_value_hash].chi + it->wavelength_grid[wl_value_hash].chi) * std::abs(it->gridvoxel->z - it_prev->gridvoxel->z) / std::abs(it->mu));
        }
    }
}

void Ray::calc_SC_coeffs(const std::size_t wl_value_hash) {
    for (std::vector<RayData>::iterator it = raydata.begin(); it != raydata.end(); ++it) {

        if (it == raydata.begin()) {
            it->wavelength_grid[wl_value_hash].alpha = 0.0;
            it->wavelength_grid[wl_value_hash].beta = 0.0;
            it->wavelength_grid[wl_value_hash].gamma = 0.0;
            it->wavelength_grid[wl_value_hash].Delta_tau = 0.0;
        } else {
            std::vector<RayData>::iterator it_prev = it;
            std::advance(it_prev, -1); // use std::prev when more C++ compilers are C++11-compliant
            it->wavelength_grid[wl_value_hash].Delta_tau = it->wavelength_grid[wl_value_hash].tau - it_prev->wavelength_grid[wl_value_hash].tau;
            // Sometimes optical depths get so large that the difference
            // between two of them amounts to the difference between two
            // extremely large numbers, and we can't keep enough significant
            // digits, so the result is numerically 0. Without using a
            // multiprecision library, which will be slow as shit, the only way
            // to fix this is by lowering the maximum optical depth. This
            // should not pose any problems though because there is no
            // interesting physics which happens at an optical depth of 10^8
            // which doesn't already happen at 10^2.
            if (std::fabs(it->wavelength_grid[wl_value_hash].Delta_tau) < std::numeric_limits<double>::epsilon()) {
                std::cerr << "ERROR: Delta_tau = 0! Your maximum optical depth is too high for floating-point precision!" << std::endl;
                std::cerr << "The first unresolvable depth point is at z = " << it->gridvoxel->z << std::endl;
                exit(1);
            }
            it->wavelength_grid[wl_value_hash].alpha = 1.0 - std::exp(-it->wavelength_grid[wl_value_hash].Delta_tau) - ((it->wavelength_grid[wl_value_hash].Delta_tau - 1.0 + std::exp(-it->wavelength_grid[wl_value_hash].Delta_tau)) / it->wavelength_grid[wl_value_hash].Delta_tau);
            it->wavelength_grid[wl_value_hash].beta = (it->wavelength_grid[wl_value_hash].Delta_tau - 1.0 + std::exp(-it->wavelength_grid[wl_value_hash].Delta_tau)) / it->wavelength_grid[wl_value_hash].Delta_tau;
            // TODO: fill in gamma for parabolic interpolation
            it->wavelength_grid[wl_value_hash].gamma = 0.0;
        }
    }
}

void Ray::formal_soln(const std::size_t wl_value_hash) {
    for (std::vector<RayData>::iterator it = raydata.begin(); it != raydata.end(); ++it) {

        if (it == raydata.begin()) {
            if (it->mu > 0.0) {
                it->wavelength_grid[wl_value_hash].I = planck_function(it->wavelength_grid[wl_value_hash].lambda, it->gridvoxel->temperature);
            } else {
                it->wavelength_grid[wl_value_hash].I = 0.0;
            }
        } else {
            std::vector<RayData>::iterator it_prev = it;
            std::advance(it_prev, -1); // use std::prev when more C++ compilers are C++11-compliant
            const double Delta_I = (it->wavelength_grid[wl_value_hash].alpha * it_prev->wavelength_grid[wl_value_hash].source_fn) + (it->wavelength_grid[wl_value_hash].beta * it->wavelength_grid[wl_value_hash].source_fn);
            it->wavelength_grid[wl_value_hash].I = it_prev->wavelength_grid[wl_value_hash].I * std::exp(-it->wavelength_grid[wl_value_hash].Delta_tau) + Delta_I;
        }
    }
}

void RayData::calc_source_fn(const std::size_t wl_value_hash) {
    wavelength_grid[wl_value_hash].source_fn = wavelength_grid[wl_value_hash].epsilon * planck_function(wavelength_grid[wl_value_hash].lambda, gridvoxel->temperature) + (1.0 - wavelength_grid[wl_value_hash].epsilon) * gridvoxel->wavelength_grid[wl_value_hash].J;
}
