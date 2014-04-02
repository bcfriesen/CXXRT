#include <iostream>

#include <yaml-cpp/yaml.h>

#include "grid.hh"
#include "ray.hh"
#include "globals.hh"
#include "calc_moments.hh"

std::vector<struct GridVoxel> grid;

int main(int argc, char *argv[]) {

    if (argc == 1) {
        std::cout << "Usage: <executable name> <YAML control file>" << std::endl;
        exit(0);
    }

    YAML::Node config = YAML::LoadFile(argv[1]);

    const int n_depth_pts = config["n_depth_pts"].as<int>();
    grid.resize(n_depth_pts);

    double rho = 1.0e-5;
    for (GridVoxel& gv: grid) {
        gv.rho = rho;
        rho *= 2.0;
        gv.temperature = 5778.0;
    }

    double z_tmp = 1.0;
    for (auto it = grid.begin(); it != grid.end(); ++it) {
      it->z = z_tmp;
      z_tmp += 1.0;
    }

    const int n_mu_pts = config["n_mu_pts"].as<int>();
    const double mu_min = -1.0;
    const double mu_max = +1.0;

    std::vector<Ray> rays(n_mu_pts);

    double mu = mu_min;
    for (Ray& r: rays) {
      r.bind_to_grid(mu);
      for (RayData& rd: r.raydata) {
        rd.lambda = 5.0;
      }
      mu += (mu_max - mu_min) / double(rays.size());
      r.calc_chi();
      r.calc_tau();
      r.formal_soln();
    }

    for (GridVoxel& gv: grid) {
      calc_J(gv);
    }

    return(0);
}
