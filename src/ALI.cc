#include <fstream>
#include <iomanip>

#include <Eigen/Sparse>

#include "calc_ALO.hh"
#include "globals.hh"
#include "grid.hh"
#include "planck_function.hh"
#include "rmsd.hh"
#include "wavelength_grid.hh"

void do_ALI() {

    const unsigned int n_depth_pts = config["n_depth_pts"].as<int>();
    const unsigned int max_iter = config["max_iter"].as<int>();
    const double epsilon = config["epsilon"].as<double>();

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

            std::vector< Eigen::Triplet<double> > tripletList;
            tripletList.reserve(n_depth_pts);

            for (unsigned int k = 0; k < n_depth_pts; ++k) {
                tripletList.push_back(Eigen::Triplet<double> (k, k, 1.0 - (1.0 - epsilon)*(Lambda_star(k, k))));
            }
            Eigen::SparseMatrix<double> mtx(n_depth_pts, n_depth_pts);
            mtx.setFromTriplets(tripletList.begin(), tripletList.end());
            Eigen::ConjugateGradient< Eigen::SparseMatrix<double> > cg;
            cg.compute(mtx);
            // Iterate the conjugate gradient method on the ALI linear system until it converges.
            unsigned int n_cg_iter = 0;
            const unsigned int max_cg_iter = 20;
            const double max_err_tol = 1.0e-10;
            while (true) {
                J_new = cg.solve(rhs);
                n_cg_iter++;
                if (cg.error() <= max_err_tol) break;
                if (n_cg_iter > max_cg_iter) {
                    std::cerr << "ERROR: could not solve sparse ALI system! Error after " << n_cg_iter << " CG iterations: " << cg.error() << std::endl;
                    exit(1);
                }
            }

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


}