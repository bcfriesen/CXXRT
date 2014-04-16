#include <iostream>
#include <limits>
#include <cmath>
#include <iomanip>

#include <yaml-cpp/yaml.h>
#include <Eigen/Dense>

#include "calc_ALO.hh"
#include "grid.hh"
#include "ray.hh"
#include "globals.hh"

std::vector<class GridVoxel> grid;
std::vector<Ray> rays;

int main(int argc, char *argv[]) {

    if (argc == 1) {
        std::cerr << "Usage: <executable name> <YAML control file>" << std::endl;
        exit(0);
    }

    YAML::Node config = YAML::LoadFile(argv[1]);

    std::cout << std::scientific;

    const int n_depth_pts = config["n_depth_pts"].as<int>();
    grid.resize(n_depth_pts);

    const double log10_rho_min = config["log10_rho_min"].as<double>();
    const double log10_rho_max = config["log10_rho_max"].as<double>();
    double log10_rho = log10_rho_min;
    const double log10_delta_rho = (log10_rho_max - log10_rho_min) / double(grid.size()-1);
    for (GridVoxel& gv: grid) {
        gv.rho = std::pow(10.0, log10_rho);
        gv.temperature = 5778.0;
        log10_rho += log10_delta_rho;
    }

    double z_tmp = 1.0;
    for (auto it = grid.rbegin(); it != grid.rend(); ++it) {
      it->z = z_tmp;
      z_tmp += 1.0;
    }

    const int n_mu_pts = config["n_mu_pts"].as<int>();
    rays.resize(n_mu_pts);
    const double mu_min = -1.0;
    const double mu_max = +1.0;
    double mu = mu_min;
    for (Ray& r: rays) {
      if (std::fabs(mu) < std::numeric_limits<double>::epsilon()) {
        std::cerr << "ERROR: mu too close to zero! : " << mu << std::endl;
        exit(0);
      }
      r.bind_to_grid(mu);
      for (RayData& rd: r.raydata) {
        rd.lambda = 5.0e-5;
        rd.epsilon = config["epsilon"].as<double>();
      }
      mu += (mu_max - mu_min) / double(rays.size());
      r.calc_chi();
      r.calc_tau();
      r.set_to_LTE();
      r.calc_SC_coeffs();
      r.formal_soln();
    }

    for (GridVoxel& gv: grid) {
      gv.calc_J();
    }

    if (config["print_every_iter"].as<bool>()) {
      for (GridVoxel& gv: grid) {
        std::cout << std::setw(15) << gv.rho << std::setw(15) << gv.J_lam << std::endl;
      }
      std::cout << std::endl;
    }

    Eigen::MatrixXd Lambda_star = calc_ALO();

    Eigen::VectorXd J_old(grid.size());
    Eigen::VectorXd J_new(grid.size());
    Eigen::VectorXd J_fs(grid.size());
    for (unsigned int i = 0; i < grid.size(); ++i) {
        J_old(i) = grid.at(i).J_lam;
    }
    Eigen::VectorXd rhs;
    Eigen::MatrixXd mtx;
    for (int i = 0; i < config["max_iter"].as<int>(); ++i) {
        for (Ray& r: rays) {
          r.calc_source_fn();
          r.formal_soln();
        }
        for (GridVoxel& gv: grid) {
          gv.calc_J();
        }
        for (unsigned int i = 0; i < grid.size(); ++i) {
            J_fs(i) = grid.at(i).J_lam;
        }
        rhs = J_fs - (1.0 - config["epsilon"].as<double>())*Lambda_star*J_old;
        mtx = Eigen::MatrixXd::Identity(grid.size(), grid.size()) - (1.0 - config["epsilon"].as<double>())*Lambda_star;
        J_new = mtx.colPivHouseholderQr().solve(rhs);
        J_old = J_new;
        for (unsigned int i = 0; i < grid.size(); ++i) {
          grid.at(i).J_lam = J_old(i);
        }

        if (config["print_every_iter"].as<bool>()) {
          for (GridVoxel& gv: grid) {
            std::cout << std::setw(15) << gv.rho << std::setw(15) << gv.J_lam << std::endl;
          }
          std::cout << std::endl;
        }
    }

    return(0);
}
