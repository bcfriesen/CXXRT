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

    for (std::map<std::size_t, double>::const_iterator wlv = wavelength_values.begin(); wlv != wavelength_values.end(); ++wlv) {
        log_file << "Starting Lambda iteration on wavelength point " << wlv->second*1.0e+8 << " A ... ";
        unsigned int iter = 0;

        for (unsigned int i = 0; i < n_depth_pts; ++i) {
            J_old(i) = grid.at(i).wavelength_grid[wlv->first].J;
        }

        do {
            for (std::vector<Ray>::iterator r = rays.begin(); r != rays.end(); ++r) {
              for (std::vector<GridVoxel>::iterator gv = grid.begin(); gv != grid.end(); ++gv) {
                gv->calc_source_fn(wlv->first);
              }
              r->formal_soln(wlv->first);
            }

            for (std::vector<GridVoxel>::iterator gv = grid.begin(); gv != grid.end(); ++gv) {
              for (std::map<std::size_t, double>::const_iterator wlp = wavelength_values.begin(); wlp != wavelength_values.end(); ++wlp) {
                  gv->calc_J(wlp->first);
              }
            }

            for (unsigned int i = 0; i < n_depth_pts; ++i) {
                J_new(i) = grid.at(i).wavelength_grid[wlv->first].J;
            }
            rmsd = calc_rmsd(J_old, J_new);
            J_old = J_new;
            std::cout << rmsd << std::endl;
            iter++;
        } while (rmsd > max_tol);
        log_file << " Converged to relative error: " << rmsd << " after " << iter << " iterations." << std::endl;
    }
}
