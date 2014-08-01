#include <cmath>

#include <Eigen/Dense>

#include "globals.hh"
#include "ray.hh"
#include "wavelength_grid.hh"

Eigen::MatrixXd calc_ALO (const double lambda) {
    const unsigned int n_depth_pts = config["n_depth_pts"].as<int>();
    Eigen::MatrixXd Lambda_star(n_depth_pts, n_depth_pts);

    std::vector<double>  Lambda_star_contrib;
    for (unsigned int i = 0; i < n_depth_pts; ++i) {
        std::vector<double> I_hat;
        std::vector<double> mu;
        for (struct RayIntersectionData& rid: grid.at(i).ray_intersection_data) {
            auto it = rid.ray->raydata.begin() + rid.intersection_point;

            // Find the requested wavelength point.
            // TODO: make this faster than a crude linear search.
            std::vector<RayWavelengthPoint>::const_iterator wlp_it;
            for (wlp_it = it->wavelength_grid.begin(); wlp_it != it->wavelength_grid.end(); ++wlp_it) {
              if (std::abs(*(wlp_it->lambda) - lambda) < std::numeric_limits<double>::epsilon())
                  break;
            }

            if (it == rid.ray->raydata.begin()) {
              I_hat.push_back(wlp_it->beta);
            } else {
              auto it_prev = it; std::advance(it_prev, -1); // use std::prev when more C++ compilers are C++11-compliant

              std::vector<RayWavelengthPoint>::const_iterator wlp_it_prev;
              for (wlp_it_prev = it_prev->wavelength_grid.begin(); wlp_it_prev != it_prev->wavelength_grid.end(); ++wlp_it_prev) {
                if (std::abs(*(wlp_it_prev->lambda) - lambda) < std::numeric_limits<double>::epsilon())
                    break;
              }

              I_hat.push_back(wlp_it_prev->gamma * std::exp(-wlp_it_prev->Delta_tau) + wlp_it->beta);
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
