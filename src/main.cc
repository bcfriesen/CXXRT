#include <iostream>
#include <limits>
#include <cmath>
#include <iomanip>
#include <fstream>

#include <yaml-cpp/yaml.h>
#include <Eigen/Dense>

#include "EOS/atoms.hh"
#include "EOS/LTE_EOS.hh"
#include "calc_ALO.hh"
#include "grid.hh"
#include "ray.hh"
#include "globals.hh"
#include "rmsd.hh"
#include "constants.hh"
#include "planck_function.hh"

std::vector<class GridVoxel> grid;
std::vector<Ray> rays;
YAML::Node config;
std::vector<double> wavelength_values;

int main(int argc, char *argv[]) {

    if (argc == 1) {
        std::cerr << "Usage: <executable name> <YAML control file>" << std::endl;
        exit(1);
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

    // set up the wavelength grid
    const int n_wavelength_pts = config["n_wavelength_pts"].as<int>();
    const double wl_min = config["wl_min"].as<double>();
    const double wl_max = config["wl_max"].as<double>();
    for (unsigned int i = 0; i < n_wavelength_pts; ++i) {
        wavelength_values.push_back(wl_min + double(i) * (wl_max - wl_min) / double(n_wavelength_pts-1));
    }

    std::ofstream log_file;
    log_file.open(log_file_name.c_str());
    log_file << std::scientific;

    log_file << "PARAMETERS USED:" << std::endl;
    log_file << config << std::endl << std::endl;

    for (GridVoxel& gv: grid) {
        gv.rho = std::pow(10.0, log10_rho);
        gv.temperature = 5778.0;
        log10_rho += log10_delta_rho;
        for (double &wlv: wavelength_values) {
            GridWavelengthPoint gwlp_tmp;
            gwlp_tmp.lambda = &wlv;
            gv.wavelength_grid.push_back(gwlp_tmp);
        }
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
        exit(1);
      }
      r.bind_to_grid(mu);
      for (RayData& rd: r.raydata) {
        for (RayWavelengthPoint& rwlp: rd.wavelength_grid) {
          rwlp.epsilon = epsilon;
        }
      }
      mu += (mu_max - mu_min) / double(rays.size()-1);

      for (RayData& rd: r.raydata) {
        for (RayWavelengthPoint& rwlp: rd.wavelength_grid) {
          rwlp.calc_chi(rd.gridvoxel->rho, *(rwlp.lambda));
        }
      }

      for (auto wlv: wavelength_values) {
        r.calc_tau(wlv);
        r.calc_SC_coeffs(wlv);
      }

      for (RayData& rd: r.raydata) {
        for (RayWavelengthPoint& rwlp: rd.wavelength_grid) {
          rwlp.set_to_LTE(rd.gridvoxel->temperature);
        }
      }

      for (auto wlv: wavelength_values) {
        r.formal_soln(wlv);
      }
    }

    for (GridVoxel& gv: grid) {
      for (GridWavelengthPoint& wlp: gv.wavelength_grid) {
        gv.calc_J(*(wlp.lambda));
        gv.calc_H(*(wlp.lambda));
        gv.calc_K(*(wlp.lambda));
      }
    }

    std::ofstream moments_file;
    moments_file.open(moments_file_name.c_str());
    moments_file << std::scientific;
    moments_file << "#" << std::setw(15) << "z" << std::setw(15) << "rho" << std::setw(15) << "lambda" << std::setw(15) << "J_lam" << std::setw(15) << "H_lam" << std::setw(15) << "K_lam" << std::setw(15) << "B_lam" << std::endl;

    if (config["print_every_iter"].as<bool>()) {
      for (GridVoxel& gv: grid) {
        for (GridWavelengthPoint& wlp: gv.wavelength_grid) {
          moments_file << std::setw(16) << gv.z << std::setw(15) << gv.rho << std::setw(15) << *(wlp.lambda) << std::setw(15) << wlp.J << std::setw(15) << wlp.H << std::setw(15) << wlp.K << std::setw(15) << planck_function(*(wlp.lambda), gv.temperature) << std::endl;
        }
      }
      moments_file << std::endl;
    }

    for (double wlv: wavelength_values) {
        Eigen::MatrixXd Lambda_star = calc_ALO(wlv);

        Eigen::VectorXd J_old(n_depth_pts);
        Eigen::VectorXd J_new(n_depth_pts);
        Eigen::VectorXd J_fs(n_depth_pts);
        for (unsigned int i = 0; i < n_depth_pts; ++i) {
            // Find the requested wavelength point on the grid voxel.
            // TODO: make this faster than a crude linear search.
            std::vector<GridWavelengthPoint>::iterator grid_wlp;
            for (grid_wlp = grid.at(i).wavelength_grid.begin(); grid_wlp != grid.at(i).wavelength_grid.end(); ++grid_wlp) {
              if (std::abs(*(grid_wlp->lambda) - wlv) < std::numeric_limits<double>::epsilon())
                  break;
            }
            J_old(i) = grid_wlp->J;
        }
        Eigen::VectorXd rhs;
        Eigen::MatrixXd mtx;

        log_file << "Beginning ALI ..." << std::endl;
        for (unsigned int i = 0; i < max_iter; ++i) {
            for (Ray& r: rays) {
              for (RayData &rd: r.raydata) {
                rd.calc_source_fn(wlv);
              }
              r.formal_soln(wlv);
            }
            for (GridVoxel& gv: grid) {
              gv.calc_J(wlv);
              gv.calc_H(wlv);
              gv.calc_K(wlv);
            }
            for (unsigned int j = 0; j < n_depth_pts; ++j) {
                // Find the requested wavelength point on the grid voxel.
                // TODO: make this faster than a crude linear search.
                std::vector<GridWavelengthPoint>::iterator grid_wlp;
                for (grid_wlp = grid.at(j).wavelength_grid.begin(); grid_wlp != grid.at(j).wavelength_grid.end(); ++grid_wlp) {
                  if (std::abs(*(grid_wlp->lambda) - wlv) < std::numeric_limits<double>::epsilon())
                      break;
                }
                J_fs(j) = grid_wlp->J;
            }
            rhs = J_fs - (1.0 - epsilon)*Lambda_star*J_old;
            mtx = Eigen::MatrixXd::Identity(n_depth_pts, n_depth_pts) - (1.0 - epsilon)*Lambda_star;
            J_new = mtx.colPivHouseholderQr().solve(rhs);
            double rmsd = calc_rmsd(J_old, J_new);
            log_file << "RMSD of relative change in J: " << rmsd << std::endl;
            J_old = J_new;
            for (unsigned int j = 0; j < n_depth_pts; ++j) {
              // Find the requested wavelength point on the grid voxel.
              // TODO: make this faster than a crude linear search.
              std::vector<GridWavelengthPoint>::iterator grid_wlp;
              for (grid_wlp = grid.at(j).wavelength_grid.begin(); grid_wlp != grid.at(j).wavelength_grid.end(); ++grid_wlp) {
                if (std::abs(*(grid_wlp->lambda) - wlv) < std::numeric_limits<double>::epsilon())
                    break;
              }
              grid_wlp->J = J_old(j);
            }

            if (config["print_every_iter"].as<bool>() || i == max_iter-1) {
              for (GridVoxel& gv: grid) {
                for (GridWavelengthPoint& wlp: gv.wavelength_grid) {
                  moments_file << std::setw(16) << gv.z << std::setw(15) << gv.rho << std::setw(15) << *(wlp.lambda) << std::setw(15) << wlp.J << std::setw(15) << wlp.H << std::setw(15) << wlp.K << std::setw(15) << planck_function(*(wlp.lambda), gv.temperature) << std::endl;
                }
              }
              moments_file << std::endl;
            }
        }
    }

    for (auto &gv: grid) {
        gv.temperature = 5778.0;
        // TODO: this works only for hydrogen! fix when adding more elements!!
        gv.n_g = gv.rho / H_mass;
    }

    for (auto &gv: grid) {
        Atom H(1);
        gv.atoms.push_back(H);
    }

    for (auto &gv: grid) {
        for (auto &atom: gv.atoms) {
            atom.number_fraction = 1.0;
            for (auto &ion: atom.ions) {
                ion.calc_partition_function(gv.temperature);
            }
        }
    }

    // Normalize number fractions. They must add up to 1.
    for (auto &gv: grid) {
        double tmp = 0.0;
        for (auto &atom: gv.atoms) {
            tmp += atom.number_fraction;
        }
        for (auto &atom: gv.atoms) {
            atom.number_fraction /= tmp;
        }
    }

    log_file << std::endl;
    log_file << "GRID VALUES:" << std::endl;
    log_file << std::setw(15) << "rho" << std::setw(15) << "temperature" << std::setw(15) << "n_e" << std::endl;
    for (auto &gv: grid) {
        calc_n_e_LTE(gv);
        log_file << std::setw(15) << gv.rho << std::setw(15) << gv.temperature << std::setw(15) << gv.n_e << std::endl;
    }

    moments_file.close();
    log_file.close();

    return(0);
}
