#include <fstream>
#include <iomanip>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <Eigen/Sparse>

#include "calc_ALO.hh"
#include "globals.hh"
#include "grid.hh"
#include "planck_function.hh"
#include "rmsd.hh"
#include "wavelength_grid.hh"

void do_ALI() {

    const unsigned int n_depth_pts = config["n_depth_pts"].as<int>();

    double rmsd;
    const double max_tol = 1.0e-8;

#pragma omp parallel
    {
    for (std::map<std::size_t, double>::const_iterator wlv = wavelength_values.begin(); wlv != wavelength_values.end(); ++wlv) {
#pragma omp single nowait
        {
        // use thread-specific buffers to store output, then dump them to the
        // log file at the end. This way the text isn't garbled by racy output.
        std::stringstream thread_buf;
        thread_buf << std::scientific;

        thread_buf << "Starting ALI on wavelength point " << wlv->second * 1.0e+8 << " A ... ";
        unsigned int iter = 0;
        Eigen::MatrixXd Lambda_star = calc_ALO(wlv->first);

        Eigen::VectorXd J_old(n_depth_pts);
        Eigen::VectorXd J_new(n_depth_pts);
        Eigen::VectorXd J_fs(n_depth_pts);
        for (unsigned int i = 0; i < n_depth_pts; ++i) {
            J_old(i) = grid.at(i).wavelength_grid[wlv->first].J;
        }
        Eigen::VectorXd rhs;
        Eigen::VectorXd epsilon(n_depth_pts);

        do {
            for (Ray& r: rays) {
                for (RayData &rd: r.raydata) {
                    rd.calc_source_fn(wlv->first);
                }
                r.formal_soln(wlv->first);
            }
            for (GridVoxel& gv: grid) {
                gv.calc_J(wlv->first);
                gv.calc_H(wlv->first);
                gv.calc_K(wlv->first);
            }
            for (unsigned int j = 0; j < n_depth_pts; ++j) {
                J_fs(j) = grid.at(j).wavelength_grid[wlv->first].J;
                epsilon(j) = grid.at(j).wavelength_grid[wlv->first].epsilon;
            }
            rhs = Lambda_star*J_old;
            for (unsigned int j = 0; j < n_depth_pts; ++j) {
                rhs(j) *= (1.0 - epsilon(j));
            }
            rhs = J_fs - rhs;

            std::vector< Eigen::Triplet<double> > tripletList;
            tripletList.reserve(n_depth_pts);

            for (unsigned int k = 0; k < n_depth_pts; ++k) {
                tripletList.push_back(Eigen::Triplet<double> (k, k, 1.0 - (1.0 - epsilon(k))*(Lambda_star(k, k))));
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

            rmsd = calc_rmsd(J_old, J_new);
            J_old = J_new;
            for (unsigned int j = 0; j < n_depth_pts; ++j) {
                grid.at(j).wavelength_grid[wlv->first].J = J_old(j);
            }
            if (rmsd < max_tol)
                break;

            if (config["print_every_iter"].as<bool>()) {
                for (GridVoxel& gv: grid) {
                    for (auto &wlp: gv.wavelength_grid) {
                        moments_file << std::setw(16) << gv.z << std::setw(15) << gv.rho << std::setw(15) << wlp.second.lambda << std::setw(15) << wlp.second.J << std::setw(15) << wlp.second.H << std::setw(15) << wlp.second.K << std::setw(15) << planck_function(wlp.second.lambda, gv.temperature) << std::endl;
                    }
                }
                moments_file << std::endl;
            }
            iter++;
        } while (rmsd > max_tol);
        thread_buf << " Converged to relative error: " << rmsd << " after " << iter << " iterations." << std::endl;
#pragma omp critical
        log_file << thread_buf.rdbuf();
        } // #pragma omp single nowait
    }
    } // #pragma omp parallel


}
