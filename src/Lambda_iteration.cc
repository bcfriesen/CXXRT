#include <fstream>

#include <Eigen/Dense>

#include "grid.hh"
#include "ray.hh"
#include "globals.hh"
#include "wavelength_grid.hh"
#include "rmsd.hh"

void Lambda_iteration() {

    const unsigned int n_depth_pts = config["n_depth_pts"].as<int>();

    double rmsd;
    const double max_tol = 1.0e-8;
    Eigen::VectorXd J_old(n_depth_pts);
    Eigen::VectorXd J_new(n_depth_pts);

    for (unsigned int i = 0; i < wavelength_values.size(); ++i) {
        log_file << "Starting Lambda iteration on wavelength point " << wavelength_values.at(i)*1.0e+8 << " A ... ";
        unsigned int iter = 0;

        for (unsigned int j = 0; j < n_depth_pts; ++j) {
            J_old(j) = grid.at(j).wavelength_grid.at(i).J;
        }

        do {
            for (std::vector<Ray>::iterator r = rays.begin(); r != rays.end(); ++r) {
              for (std::vector<GridVoxel>::iterator gv = grid.begin(); gv != grid.end(); ++gv) {
                gv->calc_source_fn(i);
              }
              r->formal_soln(i);
            }

            for (std::vector<GridVoxel>::iterator gv = grid.begin(); gv != grid.end(); ++gv) {
              for (unsigned int j = 0; j < wavelength_values.size(); ++j) {
                  gv->calc_J(j);
              }
            }

            for (unsigned int j = 0; j < n_depth_pts; ++j) {
                J_new(j) = grid.at(j).wavelength_grid.at(i).J;
            }
            rmsd = calc_rmsd(J_old, J_new);
            J_old = J_new;
            std::cout << rmsd << std::endl;
            iter++;
        } while (rmsd > max_tol);
        log_file << " Converged to relative error: " << rmsd << " after " << iter << " iterations." << std::endl;
    }
}
