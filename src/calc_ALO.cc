#include <cmath>

#include <Eigen/Dense>

#include "globals.hh"
#include "ray.hh"

Eigen::MatrixXd calc_ALO (std::vector<Ray> rays) {
    Eigen::MatrixXd Lambda_star(grid.size(), grid.size());

    std::vector<double>  Lambda_star_contrib;
    for (unsigned int i = 0; i < grid.size(); ++i) {
        std::vector<double> I_hat;
        std::vector<double> mu;
        for (auto it = grid.at(i).ray_intersection_data.begin(); it != grid.at(i).ray_intersection_data.end(); ++it) {
            auto real_it = it->ray->raydata.begin() + it->intersection_point;
            if (real_it == it->ray->raydata.begin()) {
              I_hat.push_back(real_it->beta);
            } else {
              const auto it_prev = std::prev(real_it, 1);
              I_hat.push_back(it_prev->gamma * std::exp(-it_prev->Delta_tau) + real_it->beta);
            }
            mu.push_back(real_it->mu);
        }

        double result = 0.0;
        for (unsigned int i = 0; i < I_hat.size()-1; ++i) {
            result += 0.5 * (I_hat.at(i) + I_hat.at(i+1)) * (mu.at(i+1) - mu.at(i));
        }
        result *= 0.5;
        Lambda_star(i, i) = result;
    }
    return (Lambda_star);
}
