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
#include "rmsd.hh"

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
    const std::string log_file_name = config["log_file"].as<std::string>();
    const unsigned int n_depth_pts = config["n_depth_pts"].as<int>();
    const double log10_rho_min = config["log10_rho_min"].as<double>();
    const double log10_rho_max = config["log10_rho_max"].as<double>();
    const unsigned int max_iter = config["max_iter"].as<int>();
    const double epsilon = config["epsilon"].as<double>();

    double log10_rho = log10_rho_min;
    const double log10_delta_rho = (log10_rho_max - log10_rho_min) / double(n_depth_pts-1);
    grid.resize(n_depth_pts);

    std::ofstream log_file;
    log_file.open(log_file_name.c_str());
    log_file << std::scientific;

    log_file << "PARAMETERS USED:" << std::endl;
    log_file << config << std::endl << std::endl;

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
        rd.epsilon = epsilon;
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
    moments_file.open(moments_file_name.c_str());
    moments_file << std::scientific;
    moments_file << "#" << std::setw(15) << "z" << std::setw(15) << "rho" << std::setw(15) << "J_lam" << std::setw(15) << "H_lam" << std::setw(15) << "K_lam" << std::endl;

    if (config["print_every_iter"].as<bool>()) {
      for (GridVoxel& gv: grid) {
        moments_file << std::setw(16) << gv.z << std::setw(15) << gv.rho << std::setw(15) << gv.J_lam << std::setw(15) << gv.H_lam << std::setw(15) << gv.K_lam << std::endl;
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

    log_file << "Beginning ALI ..." << std::endl;
    for (unsigned int i = 0; i < max_iter; ++i) {
        for (Ray& r: rays) {
          r.calc_source_fn();
          r.formal_soln();
        }
        for (GridVoxel& gv: grid) {
          gv.calc_J();
          gv.calc_H();
          gv.calc_K();
        }
        for (unsigned int j = 0; j < n_depth_pts; ++j) {
            J_fs(j) = grid.at(j).J_lam;
        }
        rhs = J_fs - (1.0 - epsilon)*Lambda_star*J_old;
        mtx = Eigen::MatrixXd::Identity(n_depth_pts, n_depth_pts) - (1.0 - epsilon)*Lambda_star;
        J_new = mtx.colPivHouseholderQr().solve(rhs);
        double rmsd = calc_rmsd(J_old, J_new);
        log_file << "RMSD of relative change in J: " << rmsd << std::endl;
        J_old = J_new;
        for (unsigned int j = 0; j < n_depth_pts; ++j) {
          grid.at(j).J_lam = J_old(j);
        }

        if (config["print_every_iter"].as<bool>() || i == max_iter-1) {
          for (GridVoxel& gv: grid) {
            moments_file << std::setw(16) << gv.z << std::setw(15) << gv.rho << std::setw(15) << gv.J_lam << std::setw(15) << gv.H_lam << std::setw(15) << gv.K_lam << std::endl;
          }
          moments_file << std::endl;
        }
    }

    Atom H(1);

    moments_file.close();
    log_file.close();

    return(0);
}
