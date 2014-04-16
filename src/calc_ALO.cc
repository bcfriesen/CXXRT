#include <cmath>

#include <Eigen/Dense>

#include "globals.hh"
#include "ray.hh"

Eigen::MatrixXd calc_ALO () {
    Eigen::MatrixXd Lambda_star(grid.size(), grid.size());

    std::vector<double>  Lambda_star_contrib;
    for (unsigned int i = 0; i < grid.size(); ++i) {
        std::vector<double> I_hat;
        std::vector<double> mu;
        for (struct RayIntersectionData& rid: grid.at(i).ray_intersection_data) {
            auto it = rid.ray->raydata.begin() + rid.intersection_point;
            if (it == rid.ray->raydata.begin()) {
              I_hat.push_back(it->beta);
            } else {
              auto it_prev = it; std::advance(it_prev, -1); // use std::prev when more C++ compilers are C++11-compliant
              I_hat.push_back(it_prev->gamma * std::exp(-it_prev->Delta_tau) + it->beta);
            }
            mu.push_back(it->mu);
        }

        double result = 0.0;
        for (unsigned int j = 0; j < I_hat.size()-1; ++j) {
            result += 0.5 * (I_hat.at(j) + I_hat.at(j+1)) * (mu.at(j+1) - mu.at(j));
        }
        result *= 0.5;
        Lambda_star(i, i) = result;
    }
    return (Lambda_star);
}
