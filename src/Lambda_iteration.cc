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

    for (double lambda: wavelength_values) {
        log_file << "Starting Lambda iteration on wavelength point " << lambda*1.0e+8 << " A ... ";
        unsigned int iter = 0;

        for (unsigned int i = 0; i < n_depth_pts; ++i) {
            // Find the requested wavelength point on the grid voxel.
            // TODO: make this faster than a crude linear search.
            std::vector<GridWavelengthPoint>::iterator grid_wlp;
            for (grid_wlp = grid.at(i).wavelength_grid.begin(); grid_wlp != grid.at(i).wavelength_grid.end(); ++grid_wlp) {
                if (std::abs(*(grid_wlp->lambda) - lambda) < std::numeric_limits<double>::epsilon())
                    break;
            }
            J_old(i) = grid_wlp->J;
        }

        do {
            for (Ray& r: rays) {
              for (RayData& rd: r.raydata) {
                rd.calc_source_fn(lambda);
              }
              r.formal_soln(lambda);
            }

            for (auto gv: grid) {
              for (auto wlp: gv.wavelength_grid) {
                  gv.calc_J(*(wlp.lambda));
              }
            }

            for (unsigned int i = 0; i < n_depth_pts; ++i) {
                // Find the requested wavelength point on the grid voxel.
                // TODO: make this faster than a crude linear search.
                std::vector<GridWavelengthPoint>::iterator grid_wlp;
                for (grid_wlp = grid.at(i).wavelength_grid.begin(); grid_wlp != grid.at(i).wavelength_grid.end(); ++grid_wlp) {
                    if (std::abs(*(grid_wlp->lambda) - lambda) < std::numeric_limits<double>::epsilon())
                        break;
                }
                J_new(i) = grid_wlp->J;
            }
            rmsd = calc_rmsd(J_old, J_new);
            J_old = J_new;
            std::cout << rmsd << std::endl;
            iter++;
        } while (rmsd > max_tol);
        log_file << " Converged to relative error: " << rmsd << " after " << iter << " iterations." << std::endl;
    }
}
