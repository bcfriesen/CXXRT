#include <cmath>

#include <Eigen/Sparse>

#include "globals.hh"
#include "ray.hh"
#include "wavelength_grid.hh"

Eigen::SparseMatrix<double> calc_ALO (const unsigned int wl_index) {
    const unsigned int n_depth_pts = grid.size();
    Eigen::SparseMatrix<double> Lambda_star(n_depth_pts, n_depth_pts);

    std::vector< Eigen::Triplet<double> > tripletList;
    tripletList.reserve(n_depth_pts);

    for (unsigned int i = 0; i < n_depth_pts; ++i) {
        std::vector<double> I_hat;
        std::vector<double> mu;
        for (std::vector<RayIntersectionData>::const_iterator rid = grid.at(i).ray_intersection_data.begin(); rid != grid.at(i).ray_intersection_data.end(); ++rid) {
            std::vector<RayData>::iterator it = rid->ray->raydata.begin() + rid->intersection_point;

            if (it == rid->ray->raydata.begin()) {
                I_hat.push_back(it->wavelength_grid.at(wl_index).beta);
            } else {
                std::vector<RayData>::iterator it_prev = it;
                std::advance(it_prev, -1); // use std::prev when more C++ compilers are C++11-compliant
                I_hat.push_back(it_prev->wavelength_grid.at(wl_index).gamma * std::exp(-it_prev->wavelength_grid.at(wl_index).Delta_tau) + it->wavelength_grid.at(wl_index).beta);
            }
            mu.push_back(it->mu);
        }

        double result = 0.0;
        for (unsigned int j = 0; j < I_hat.size()-1; ++j) {
            result += 0.5 * (I_hat.at(j) + I_hat.at(j+1)) * (mu.at(j+1) - mu.at(j));
        }
        result *= 0.5;
        tripletList.push_back(Eigen::Triplet<double> (i, i, result));
    }
    Lambda_star.setFromTriplets(tripletList.begin(), tripletList.end());
    return (Lambda_star);
}
