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
    for (RayData& rd: raydata) {
        rd.mu = mu;
    }
    int i = 0;
    for (auto it = raydata.begin(); it != raydata.end(); ++it) {
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

void Ray::print_ray_data(const double lambda) {
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

    for (auto rd: raydata) {
        std::vector<RayWavelengthPoint>::iterator wlp_it;
        for (wlp_it = rd.wavelength_grid.begin(); wlp_it != rd.wavelength_grid.end(); ++wlp_it) {
            if (std::abs(*(wlp_it->lambda) - lambda) < std::numeric_limits<double>::epsilon())
                break;
        }

        log_file << std::setw(30) << std::scientific << rd.gridvoxel->z;
        log_file << std::setw(15) << std::scientific << rd.mu;
        log_file << std::setw(15) << std::scientific << *(wlp_it->lambda) * 1.0e+8;
        log_file << std::setw(15) << std::scientific << wlp_it->I;
        log_file << std::setw(15) << std::scientific << wlp_it->sigma;
        log_file << std::setw(15) << std::scientific << wlp_it->kappa;
        log_file << std::setw(15) << std::scientific << wlp_it->eta;
        log_file << std::setw(15) << std::scientific << wlp_it->tau;
        log_file << std::setw(15) << std::scientific << wlp_it->chi;
        log_file << std::setw(15) << std::scientific << wlp_it->source_fn;
        log_file << std::setw(15) << std::scientific << wlp_it->Delta_tau;
        log_file << std::setw(15) << std::scientific << wlp_it->alpha;
        log_file << std::setw(15) << std::scientific << wlp_it->beta;
        log_file << std::setw(15) << std::scientific << wlp_it->gamma;
        log_file << std::endl;
    }
    log_file << std::endl;
}


RayData::RayData() {
    gridvoxel = nullptr;

    for (auto wlv = wavelength_values.begin(); wlv != wavelength_values.end(); ++wlv) {
        RayWavelengthPoint wlp_tmp;
        wlp_tmp.lambda = &(*wlv);
        wavelength_grid.push_back(wlp_tmp);
    }

    for (RayWavelengthPoint& wlp: wavelength_grid) {
        wlp.I = 0.0;
        wlp.tau = 0.0;
        wlp.chi = 0.0;
        wlp.source_fn = 0.0;
    }
}

void Ray::calc_tau(const double lambda) {
    for (auto it = raydata.begin(); it != raydata.end(); ++it) {

        // Find the requested wavelength point.
        // TODO: make this faster than a crude linear search.
        std::vector<RayWavelengthPoint>::iterator wlp_it;
        for (wlp_it = it->wavelength_grid.begin(); wlp_it != it->wavelength_grid.end(); ++wlp_it) {
            if (std::abs(*(wlp_it->lambda) - lambda) < std::numeric_limits<double>::epsilon())
                break;
        }

        if (it == raydata.begin()) {
            wlp_it->tau = 0.0;
        } else {
            auto it_prev = it;
            std::advance(it_prev, -1); // use std::prev when more C++ compilers are C++11-compliant

            std::vector<RayWavelengthPoint>::iterator wlp_it_prev;
            for (wlp_it_prev = it_prev->wavelength_grid.begin(); wlp_it_prev != it_prev->wavelength_grid.end(); ++wlp_it_prev) {
                if (std::abs(*(wlp_it_prev->lambda) - lambda) < std::numeric_limits<double>::epsilon())
                    break;
            }

            wlp_it->tau = wlp_it_prev->tau + (0.5 * (wlp_it_prev->chi + wlp_it->chi) * std::abs(it->gridvoxel->z - it_prev->gridvoxel->z) / std::abs(it->mu));
        }
    }
}

void Ray::calc_SC_coeffs(const double lambda) {
    for (auto it = raydata.begin(); it != raydata.end(); ++it) {

        // Find the requested wavelength point.
        // TODO: make this faster than a crude linear search.
        std::vector<RayWavelengthPoint>::iterator wlp_it;
        for (wlp_it = it->wavelength_grid.begin(); wlp_it != it->wavelength_grid.end(); ++wlp_it) {
            if (std::abs(*(wlp_it->lambda) - lambda) < std::numeric_limits<double>::epsilon())
                break;
        }

        if (it == raydata.begin()) {
            wlp_it->alpha = 0.0;
            wlp_it->beta = 0.0;
            wlp_it->gamma = 0.0;
            wlp_it->Delta_tau = 0.0;
        } else {
            auto it_prev = it;
            std::advance(it_prev, -1); // use std::prev when more C++ compilers are C++11-compliant

            std::vector<RayWavelengthPoint>::iterator wlp_it_prev;
            for (wlp_it_prev = it_prev->wavelength_grid.begin(); wlp_it_prev != it_prev->wavelength_grid.end(); ++wlp_it_prev) {
                if (std::abs(*(wlp_it_prev->lambda) - lambda) < std::numeric_limits<double>::epsilon())
                    break;
            }

            wlp_it->Delta_tau = wlp_it->tau - wlp_it_prev->tau;
            // Sometimes optical depths get so large that the difference
            // between two of them amounts to the difference between two
            // extremely large numbers, and we can't keep enough significant
            // digits, so the result is numerically 0. Without using a
            // multiprecision library, which will be slow as shit, the only way
            // to fix this is by lowering the maximum optical depth. This
            // should not pose any problems though because there is no
            // interesting physics which happens at an optical depth of 10^8
            // which doesn't already happen at 10^2.
            if (std::fabs(wlp_it->Delta_tau) < std::numeric_limits<double>::epsilon()) {
                std::cerr << "ERROR: Delta_tau = 0! Your maximum optical depth is too high for floating-point precision!" << std::endl;
                std::cerr << "The first unresolvable depth point is at z = " << it->gridvoxel->z << std::endl;
                exit(1);
            }
            wlp_it->alpha = 1.0 - std::exp(-wlp_it->Delta_tau) - ((wlp_it->Delta_tau - 1.0 + std::exp(-wlp_it->Delta_tau)) / wlp_it->Delta_tau);
            wlp_it->beta = (wlp_it->Delta_tau - 1.0 + std::exp(-wlp_it->Delta_tau)) / wlp_it->Delta_tau;
            // TODO: fill in gamma for parabolic interpolation
            wlp_it->gamma = 0.0;
        }
    }
}

void Ray::formal_soln(const double lambda) {
    for (auto it = raydata.begin(); it != raydata.end(); ++it) {

        // Find the requested wavelength point.
        // TODO: make this faster than a crude linear search.
        std::vector<RayWavelengthPoint>::iterator wlp_it;
        for (wlp_it = it->wavelength_grid.begin(); wlp_it != it->wavelength_grid.end(); ++wlp_it) {
            if (std::abs(*(wlp_it->lambda) - lambda) < std::numeric_limits<double>::epsilon())
                break;
        }

        if (it == raydata.begin()) {
            if (it->mu > 0.0) {
                wlp_it->I = planck_function(*(wlp_it->lambda), it->gridvoxel->temperature);
            } else {
                wlp_it->I = 0.0;
            }
        } else {
            auto it_prev = it;
            std::advance(it_prev, -1); // use std::prev when more C++ compilers are C++11-compliant

            std::vector<RayWavelengthPoint>::iterator wlp_it_prev;
            for (wlp_it_prev = it_prev->wavelength_grid.begin(); wlp_it_prev != it_prev->wavelength_grid.end(); ++wlp_it_prev) {
                if (std::abs(*(wlp_it_prev->lambda) - lambda) < std::numeric_limits<double>::epsilon())
                    break;
            }

            const double Delta_I = (wlp_it->alpha * wlp_it_prev->source_fn) + (wlp_it->beta * wlp_it->source_fn);
            wlp_it->I = wlp_it_prev->I * std::exp(-wlp_it->Delta_tau) + Delta_I;
        }
    }
}

void RayData::calc_source_fn(const double lambda) {

    // Find the requested wavelength point.
    // TODO: make this faster than a crude linear search.
    std::vector<RayWavelengthPoint>::iterator wlp;
    for (wlp = wavelength_grid.begin(); wlp != wavelength_grid.end(); ++wlp) {
        if (std::abs(*(wlp->lambda) - lambda) < std::numeric_limits<double>::epsilon())
            break;
    }

    std::vector<GridWavelengthPoint>::const_iterator grid_wlp;
    for (grid_wlp = gridvoxel->wavelength_grid.begin(); grid_wlp != gridvoxel->wavelength_grid.end(); ++grid_wlp) {
        if (std::abs(*(grid_wlp->lambda) - lambda) < std::numeric_limits<double>::epsilon())
            break;
    }

    wlp->source_fn = wlp->epsilon * planck_function(*(wlp->lambda), gridvoxel->temperature) + (1.0 - wlp->epsilon) * grid_wlp->J;
}
