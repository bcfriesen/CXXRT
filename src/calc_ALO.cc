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
        std::vector<double> I_hat_i;
        std::vector<double> I_hat_im1;
        std::vector<double> I_hat_ip1;
        std::vector<double> mu_i;
        std::vector<double> mu_im1;
        std::vector<double> mu_ip1;
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
                I_hat_im1.push_back(0.0);
                I_hat_i.push_back(wlp_it->beta);
                I_hat_ip1.push_back(wlp_it->beta * std::exp(-wlp_it->Delta_tau) + wlp_it_next->alpha);
                mu_im1.push_back(0.0);
                mu_i.push_back(it->mu);
                mu_ip1.push_back(it_next->mu);
            } else {
                I_hat_im1.push_back(wlp_it_prev->gamma);
                I_hat_i.push_back(wlp_it_prev->gamma * std::exp(-wlp_it_prev->Delta_tau) + wlp_it->beta);
                if (it == rid.ray->raydata.end()-1) {
                    I_hat_ip1.push_back((wlp_it_prev->gamma * std::exp(-wlp_it_prev->Delta_tau) + wlp_it->beta) * std::exp(-wlp_it->Delta_tau));
                    mu_ip1.push_back(0.0);
                } else {
                    I_hat_ip1.push_back((wlp_it_prev->gamma * std::exp(-wlp_it_prev->Delta_tau) + wlp_it->beta) * std::exp(-wlp_it->Delta_tau) + wlp_it_next->alpha);
                    mu_ip1.push_back(it_next->mu);
                }
                mu_im1.push_back(it_prev->mu);
                mu_i.push_back(it->mu);
            }
        }

        double result = 0.0;

        if (i > 0) {
            for (unsigned int j = 0; j < I_hat_im1.size()-1; ++j) {
                result += 0.5 * (I_hat_im1.at(j) + I_hat_im1.at(j+1)) * (mu_im1.at(j+1) - mu_im1.at(j));
            }
            result *= 0.5;
            tripletList.push_back(Eigen::Triplet<double> (i-1, i, result));
        }

        result = 0.0;
        for (unsigned int j = 0; j < I_hat_i.size()-1; ++j) {
            result += 0.5 * (I_hat_i.at(j) + I_hat_i.at(j+1)) * (mu_i.at(j+1) - mu_i.at(j));
        }
        result *= 0.5;
        tripletList.push_back(Eigen::Triplet<double> (i, i, result));

        if (i < n_depth_pts-1) {
            result = 0.0;
            for (unsigned int j = 0; j < I_hat_ip1.size()-1; ++j) {
                result += 0.5 * (I_hat_ip1.at(j) + I_hat_ip1.at(j+1)) * (mu_ip1.at(j+1) - mu_ip1.at(j));
            }
            result *= 0.5;
            tripletList.push_back(Eigen::Triplet<double> (i+1, i, result));
        }
    }
    Lambda_star.setFromTriplets(tripletList.begin(), tripletList.end());
    std::cout << Lambda_star << std::endl;
    return (Lambda_star);
}
