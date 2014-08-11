#include <cmath>

#include <Eigen/Sparse>

#include "globals.hh"
#include "ray.hh"
#include "wavelength_grid.hh"

Eigen::SparseMatrix<double> calc_ALO (const double lambda) {
    const unsigned int n_depth_pts = config["n_depth_pts"].as<int>();
    Eigen::SparseMatrix<double> Lambda_star(n_depth_pts, n_depth_pts);

    std::vector< Eigen::Triplet<double> > tripletList;
    tripletList.reserve(3*n_depth_pts-2);

    std::vector<double>  Lambda_star_contrib;
    for (unsigned int i = 0; i < n_depth_pts; ++i) {
        for (struct RayIntersectionData& rid: grid.at(i).ray_intersection_data) {
            auto it = rid.ray->raydata.begin() + rid.intersection_point;

            // Find the requested wavelength point.
            // TODO: make this faster than a crude linear search.
            std::vector<RayWavelengthPoint>::const_iterator wlp_it;
            for (wlp_it = it->wavelength_grid.begin(); wlp_it != it->wavelength_grid.end(); ++wlp_it) {
                if (std::abs(*(wlp_it->lambda) - lambda) < std::numeric_limits<double>::epsilon())
                    break;
            }

            auto it_prev = it;
            std::vector<RayWavelengthPoint>::const_iterator wlp_it_prev;
            if (it != rid.ray->raydata.begin()) {
                std::advance(it_prev, -1); // use std::prev when more C++ compilers are C++11-compliant
                for (wlp_it_prev = it_prev->wavelength_grid.begin(); wlp_it_prev != it_prev->wavelength_grid.end(); ++wlp_it_prev) {
                    if (std::abs(*(wlp_it_prev->lambda) - lambda) < std::numeric_limits<double>::epsilon())
                        break;
                }
            }

            auto it_next = it;
            std::vector<RayWavelengthPoint>::const_iterator wlp_it_next;
            if (it != rid.ray->raydata.end()-1) {
                std::advance(it_next, +1); // use std::next when more C++ compilers are C++11-compliant
                for (wlp_it_next = it_next->wavelength_grid.begin(); wlp_it_next != it_next->wavelength_grid.end(); ++wlp_it_next) {
                    if (std::abs(*(wlp_it_next->lambda) - lambda) < std::numeric_limits<double>::epsilon())
                        break;
                }
            }

            if (it == rid.ray->raydata.begin()) {
                it_prev->gridvoxel->I_hat_im1.push_back(0.0);
                it->gridvoxel->I_hat_i.push_back(wlp_it->beta);
                it_next->gridvoxel->I_hat_ip1.push_back(wlp_it->beta * std::exp(-wlp_it->Delta_tau) + wlp_it_next->alpha);
                it_prev->gridvoxel->mu_im1.push_back(it->mu);
                it->gridvoxel->mu_i.push_back(it->mu);
                it_next->gridvoxel->mu_ip1.push_back(it->mu);
            } else {
                it_prev->gridvoxel->I_hat_im1.push_back(wlp_it_prev->gamma);
                it->gridvoxel->I_hat_i.push_back(wlp_it_prev->gamma * std::exp(-wlp_it_prev->Delta_tau) + wlp_it->beta);
                if (it == rid.ray->raydata.end()-1) {
                    it_next->gridvoxel->I_hat_ip1.push_back((wlp_it_prev->gamma * std::exp(-wlp_it_prev->Delta_tau) + wlp_it->beta) * std::exp(-wlp_it->Delta_tau));
                } else {
                    it_next->gridvoxel->I_hat_ip1.push_back((wlp_it_prev->gamma * std::exp(-wlp_it_prev->Delta_tau) + wlp_it->beta) * std::exp(-wlp_it->Delta_tau) + wlp_it_next->alpha);
                }
                it_prev->gridvoxel->mu_im1.push_back(it->mu);
                it->gridvoxel->mu_i.push_back(it->mu);
                it_next->gridvoxel->mu_ip1.push_back(it->mu);
            }
        }

        double result = 0.0;

        if (i > 0) {
            for (unsigned int j = 0; j < grid.at(i).I_hat_im1.size()-1; ++j) {
                result += 0.5 * (grid.at(i).I_hat_im1.at(j) + grid.at(i).I_hat_im1.at(j+1)) * (grid.at(i).mu_im1.at(j+1) - grid.at(i).mu_im1.at(j));
            }
            result *= 0.5;
            tripletList.push_back(Eigen::Triplet<double> (i-1, i, result));
            grid.at(i).I_hat_im1.erase(grid.at(i).I_hat_im1.begin(), grid.at(i).I_hat_im1.end());
        }

        result = 0.0;
        for (unsigned int j = 0; j < grid.at(i).I_hat_i.size()-1; ++j) {
            result += 0.5 * (grid.at(i).I_hat_i.at(j) + grid.at(i).I_hat_i.at(j+1)) * (grid.at(i).mu_i.at(j+1) - grid.at(i).mu_i.at(j));
        }
        result *= 0.5;
        tripletList.push_back(Eigen::Triplet<double> (i, i, result));
        grid.at(i).I_hat_i.erase(grid.at(i).I_hat_i.begin(), grid.at(i).I_hat_i.end());

        if (i < n_depth_pts-1) {
            result = 0.0;
            for (unsigned int j = 0; j < grid.at(i).I_hat_ip1.size()-1; ++j) {
                result += 0.5 * (grid.at(i).I_hat_ip1.at(j) + grid.at(i).I_hat_ip1.at(j+1)) * (grid.at(i).mu_ip1.at(j+1) - grid.at(i).mu_ip1.at(j));
            }
            result *= 0.5;
            tripletList.push_back(Eigen::Triplet<double> (i+1, i, result));
            grid.at(i).I_hat_ip1.erase(grid.at(i).I_hat_ip1.begin(), grid.at(i).I_hat_ip1.end());
        }
    }
    Lambda_star.setFromTriplets(tripletList.begin(), tripletList.end());
    std::cout << Lambda_star << std::endl;
    return (Lambda_star);
}
