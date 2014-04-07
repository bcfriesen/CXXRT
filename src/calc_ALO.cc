#include <Eigen/Dense>

#include "globals.hh"
#include "ray.hh"

Eigen::MatrixXd calc_ALO (std::vector<Ray> rays) {
    Eigen::MatrixXd Lambda_star(grid.size(), grid.size());

    std::vector<double>  Lambda_star_contrib;
    for (unsigned int i = 1; i < grid.size()-1; ++i) {
        std::vector<double> I_hat;
        std::vector<double> mu;
        for (auto it = grid.at(i).ray_intersection_data.begin(); it != grid.at(i).ray_intersection_data.end(); ++it) {
            const auto it_prev = std::prev(it, 1);
            I_hat.push_back(it_prev->data->gamma * std::exp(-it_prev->data->Delta_tau) + it->data->beta);
            mu.push_back(it->data->mu);
        }

        double result = 0.0;
        for (unsigned int i = 0; i < I_hat.size()-1; ++i) {
            result += 0.5 * (I_hat.at(i) + I_hat.at(i+1)) * (mu.at(i+1) - mu.at(i));
        }
        result *= 0.5;
        std::cout << result << std::endl;

    }
    return (Lambda_star);
}
