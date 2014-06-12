#include <iostream>
#include <limits>
#include <cmath>
#include <iomanip>
#include <fstream>

#include <yaml-cpp/yaml.h>
#include <Eigen/Dense>

#include "EOS/atoms.hh"
#include "calc_ALO.hh"
#include "grid.hh"
#include "ray.hh"
#include "globals.hh"

std::vector<class GridVoxel> grid;
std::vector<Ray> rays;
YAML::Node config;

int main(int argc, char *argv[]) {

    if (argc == 1) {
        std::cerr << "Usage: <executable name> <YAML control file>" << std::endl;
        exit(0);
    }
    config = YAML::LoadFile(argv[1]);
    const std::string moments_file_name = config["moments_file"].as<std::string>();
    const unsigned int n_depth_pts = config["n_depth_pts"].as<int>();
    const double log10_rho_min = config["log10_rho_min"].as<double>();
    const double log10_rho_max = config["log10_rho_max"].as<double>();

    double log10_rho = log10_rho_min;
    const double log10_delta_rho = (log10_rho_max - log10_rho_min) / double(n_depth_pts-1);
    grid.resize(n_depth_pts);

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

    std::cout << std::scientific;

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
      mu += (mu_max - mu_min) / double(rays.size()-1);
      r.calc_chi();
      r.calc_tau();
      r.calc_SC_coeffs();
      r.set_to_LTE();
      r.formal_soln();
    }

    for (GridVoxel& gv: grid) {
      gv.calc_J();
      gv.calc_H();
      gv.calc_K();
    }

    std::ofstream moments_file;
    moments_file.open(moments_file_name);
    moments_file << std::scientific;
    moments_file << std::setw(15) << "z" << std::setw(15) << "rho" << std::setw(15) << "J_lam" << std::setw(15) << "H_lam" << std::setw(15) << "K_lam" << std::endl;

    if (config["print_every_iter"].as<bool>()) {
      for (GridVoxel& gv: grid) {
        moments_file << std::setw(15) << gv.z << std::setw(15) << gv.rho << std::setw(15) << gv.J_lam << std::setw(15) << gv.H_lam << std::setw(15) << gv.K_lam << std::endl;
      }
      moments_file << std::endl;
    }

    Eigen::MatrixXd Lambda_star = calc_ALO();

    Eigen::VectorXd J_old(n_depth_pts);
    Eigen::VectorXd J_new(n_depth_pts);
    Eigen::VectorXd J_fs(n_depth_pts);
    for (unsigned int i = 0; i < n_depth_pts; ++i) {
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
          gv.calc_H();
          gv.calc_K();
        }
        for (unsigned int i = 0; i < n_depth_pts; ++i) {
            J_fs(i) = grid.at(i).J_lam;
        }
        rhs = J_fs - (1.0 - config["epsilon"].as<double>())*Lambda_star*J_old;
        mtx = Eigen::MatrixXd::Identity(n_depth_pts, n_depth_pts) - (1.0 - config["epsilon"].as<double>())*Lambda_star;
        J_new = mtx.colPivHouseholderQr().solve(rhs);
        J_old = J_new;
        for (unsigned int i = 0; i < n_depth_pts; ++i) {
          grid.at(i).J_lam = J_old(i);
        }

        if (config["print_every_iter"].as<bool>() || i == config["max_iter"].as<int>()-1) {
          for (GridVoxel& gv: grid) {
            moments_file << std::setw(15) << gv.z << std::setw(15) << gv.rho << std::setw(15) << gv.J_lam << std::setw(15) << gv.H_lam << std::setw(15) << gv.K_lam << std::endl;
          }
          moments_file << std::endl;
        }
    }

    Atom H(1);

    moments_file.close();

    return(0);
}
