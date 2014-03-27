#include <iostream>

#include <yaml-cpp/yaml.h>

#include "grid.hh"
#include "ray.hh"

int main(int argc, char *argv[]) {

    if (argc == 1) {
        std::cout << "Usage: <executable name> <YAML control file>" << std::endl;
        exit(0);
    }

    YAML::Node config = YAML::LoadFile(argv[1]);

    const int n_depth_pts = config["n_depth_pts"].as<int>();
    std::vector<struct GridVoxel> grid(n_depth_pts);

    const int n_mu_pts = config["n_mu_pts"].as<int>();
    const double mu_min = -1.0;
    const double mu_max = +1.0;

    std::vector<Ray> rays(n_mu_pts);

    {
      double mu = mu_min;
      for (Ray& r: rays) {
        r.mu = mu;
        mu += (mu_max - mu_min) / double(rays.size());
        r.bind_to_grid(grid);
      }
    }

    return(0);
}
