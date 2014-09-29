#include <algorithm>
#include <iomanip>
#include <iterator>
#include <cmath>
#include <fstream>

#include "misc/planck_function.hh"
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

void Ray::print_ray_data(const unsigned int wl_index) {
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
    log_file << std::setw(15) << "tau_Rosseland";
    log_file << std::setw(15) << "chi";
    log_file << std::setw(15) << "temperature";
    log_file << std::setw(15) << "source_fn";
    log_file << std::setw(15) << "planck_fn";
    log_file << std::setw(15) << "J_lam";
    log_file << std::setw(15) << "H_lam";
    log_file << std::setw(15) << "K_lam";
    log_file << std::setw(15) << "Delta_tau";
    log_file << std::setw(15) << "alpha";
    log_file << std::setw(15) << "beta";
    log_file << std::setw(15) << "gamma";
    log_file << std::endl;

    for (std::vector<RayData>::iterator rd = raydata.begin(); rd != raydata.end(); ++rd) {
        log_file << std::setw(30) << std::scientific << rd->gridvoxel->z;
        log_file << std::setw(15) << std::scientific << rd->mu;
        log_file << std::setw(15) << std::scientific << rd->wavelength_grid.at(wl_index).lambda * 1.0e+8;
        log_file << std::setw(15) << std::scientific << rd->wavelength_grid.at(wl_index).I;
        log_file << std::setw(15) << std::scientific << rd->gridvoxel->sigma;
        log_file << std::setw(15) << std::scientific << rd->gridvoxel->wavelength_grid.at(wl_index).kappa;
        log_file << std::setw(15) << std::scientific << rd->gridvoxel->wavelength_grid.at(wl_index).eta;
        log_file << std::setw(15) << std::scientific << rd->wavelength_grid.at(wl_index).tau;
        log_file << std::setw(15) << std::scientific << rd->tau_Rosseland;
        log_file << std::setw(15) << std::scientific << rd->gridvoxel->wavelength_grid.at(wl_index).chi;
        log_file << std::setw(15) << std::scientific << rd->gridvoxel->temperature;
        log_file << std::setw(15) << std::scientific << rd->gridvoxel->wavelength_grid.at(wl_index).source_fn;
        log_file << std::setw(15) << std::scientific << planck_function(rd->wavelength_grid.at(wl_index).lambda, rd->gridvoxel->temperature);
        log_file << std::setw(15) << std::scientific << rd->gridvoxel->wavelength_grid.at(wl_index).J;
        log_file << std::setw(15) << std::scientific << rd->gridvoxel->wavelength_grid.at(wl_index).H;
        log_file << std::setw(15) << std::scientific << rd->gridvoxel->wavelength_grid.at(wl_index).K;
        log_file << std::setw(15) << std::scientific << rd->wavelength_grid.at(wl_index).Delta_tau;
        log_file << std::setw(15) << std::scientific << rd->wavelength_grid.at(wl_index).alpha;
        log_file << std::setw(15) << std::scientific << rd->wavelength_grid.at(wl_index).beta;
        log_file << std::setw(15) << std::scientific << rd->wavelength_grid.at(wl_index).gamma;
        log_file << std::endl;
    }
    log_file << std::endl;
}


RayData::RayData() {
    gridvoxel = NULL;

    wavelength_grid.resize(wavelength_values.size());
    for (unsigned int i = 0; i < wavelength_values.size(); ++i) {
        RayWavelengthPoint wlp_tmp;
        wlp_tmp.lambda = wavelength_values.at(i);
        wavelength_grid.at(i) = wlp_tmp;
    }

    for (std::vector<RayWavelengthPoint>::iterator wlp = wavelength_grid.begin(); wlp != wavelength_grid.end(); ++wlp) {
        wlp->I = 0.0;
        wlp->tau = 0.0;
    }
}

void Ray::calc_tau(const unsigned int wl_index) {
    for (std::vector<RayData>::iterator it = raydata.begin(); it != raydata.end(); ++it) {
        if (it == raydata.begin()) {
            it->wavelength_grid.at(wl_index).tau = 0.0;
        } else {
            std::vector<RayData>::iterator it_prev = it;
            std::advance(it_prev, -1); // use std::prev when more C++ compilers are C++11-compliant
            it->wavelength_grid.at(wl_index).tau = it_prev->wavelength_grid.at(wl_index).tau + (0.5 * (it_prev->gridvoxel->wavelength_grid.at(wl_index).chi + it->gridvoxel->wavelength_grid.at(wl_index).chi) * std::abs(it->gridvoxel->z - it_prev->gridvoxel->z) / std::abs(it->mu));
        }
    }
}


void Ray::calc_tau_Rosseland() {
    for (std::vector<RayData>::iterator it = raydata.begin(); it != raydata.end(); ++it) {
        if (it == raydata.begin()) {
            it->tau_Rosseland = 0.0;
        } else {
            std::vector<RayData>::iterator it_prev = it;
            std::advance(it_prev, -1); // use std::prev when more C++ compilers are C++11-compliant
            it->tau_Rosseland = it_prev->tau_Rosseland + (0.5 * (it_prev->gridvoxel->Rosseland_mean_opacity + it->gridvoxel->Rosseland_mean_opacity) * std::abs(it->gridvoxel->z - it_prev->gridvoxel->z) / std::abs(it->mu));
        }
    }
}


void Ray::calc_SC_coeffs(const unsigned int wl_index) {
    for (std::vector<RayData>::iterator it = raydata.begin(); it != raydata.end(); ++it) {

        if (it == raydata.begin()) {
            it->wavelength_grid.at(wl_index).alpha = 0.0;
            it->wavelength_grid.at(wl_index).beta = 0.0;
            it->wavelength_grid.at(wl_index).gamma = 0.0;
            it->wavelength_grid.at(wl_index).Delta_tau = 0.0;
        } else {
            std::vector<RayData>::iterator it_prev = it;
            std::advance(it_prev, -1); // use std::prev when more C++ compilers are C++11-compliant
            it->wavelength_grid.at(wl_index).Delta_tau = it->wavelength_grid.at(wl_index).tau - it_prev->wavelength_grid.at(wl_index).tau;
            // Sometimes optical depths get so large that the difference
            // between two of them amounts to the difference between two
            // extremely large numbers, and we can't keep enough significant
            // digits, so the result is numerically 0. Without using a
            // multiprecision library, which will be slow as shit, the only way
            // to fix this is by lowering the maximum optical depth. This
            // should not pose any problems though because there is no
            // interesting physics which happens at an optical depth of 10^8
            // which doesn't already happen at 10^2.
            if (std::fabs(it->wavelength_grid.at(wl_index).Delta_tau) < std::numeric_limits<double>::epsilon()) {
                std::cerr << "ERROR: Delta_tau = 0! Your maximum optical depth is too high for floating-point precision!" << std::endl;
                std::cerr << "The first unresolvable depth point is at z = " << it->gridvoxel->z << std::endl;
                exit(1);
            }
            it->wavelength_grid.at(wl_index).alpha = 1.0 - std::exp(-it->wavelength_grid.at(wl_index).Delta_tau) - ((it->wavelength_grid.at(wl_index).Delta_tau - 1.0 + std::exp(-it->wavelength_grid.at(wl_index).Delta_tau)) / it->wavelength_grid.at(wl_index).Delta_tau);
            it->wavelength_grid.at(wl_index).beta = (it->wavelength_grid.at(wl_index).Delta_tau - 1.0 + std::exp(-it->wavelength_grid.at(wl_index).Delta_tau)) / it->wavelength_grid.at(wl_index).Delta_tau;
            // TODO: fill in gamma for parabolic interpolation
            it->wavelength_grid.at(wl_index).gamma = 0.0;
        }
    }
}

void Ray::formal_soln(const unsigned int wl_index) {
    for (std::vector<RayData>::iterator it = raydata.begin(); it != raydata.end(); ++it) {

        if (it == raydata.begin()) {
            if (it->mu > 0.0) {
                std::vector<RayData>::const_iterator it_next = it;
                std::advance(it_next, +1); // use std::next when more C++ compilers are C++11-compliant
                const double B = planck_function(it->wavelength_grid.at(wl_index).lambda, it->gridvoxel->temperature);
                const double B_next = planck_function(it_next->wavelength_grid.at(wl_index).lambda, it_next->gridvoxel->temperature);
                const double tau = it->wavelength_grid.at(wl_index).tau;
                const double tau_next = it_next->wavelength_grid.at(wl_index).tau;
                // Derivative of Planck function w.r.t. optical depth. See Mihalas Eq (2-88)
                const double dB_dtau = (B - B_next) / std::abs(tau - tau_next);
                // Diffusion condition for rays coming up from depth.
                it->wavelength_grid.at(wl_index).I = planck_function(it->wavelength_grid.at(wl_index).lambda, it->gridvoxel->temperature) + (it->mu * dB_dtau);
            } else {
                it->wavelength_grid.at(wl_index).I = 0.0;
            }
        } else {
            std::vector<RayData>::const_iterator it_prev = it;
            std::advance(it_prev, -1); // use std::prev when more C++ compilers are C++11-compliant
            const double Delta_I = (it->wavelength_grid.at(wl_index).alpha * it_prev->gridvoxel->wavelength_grid.at(wl_index).source_fn) + (it->wavelength_grid.at(wl_index).beta * it->gridvoxel->wavelength_grid.at(wl_index).source_fn);
            it->wavelength_grid.at(wl_index).I = it_prev->wavelength_grid.at(wl_index).I * std::exp(-it->wavelength_grid.at(wl_index).Delta_tau) + Delta_I;
        }
    }
}
