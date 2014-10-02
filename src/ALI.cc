#include <fstream>
#include <iomanip>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <viennacl/vector.hpp>
#include <viennacl/compressed_matrix.hpp>
#include "viennacl/linalg/cg.hpp"

#include "calc_ALO.hh"
#include "globals.hh"
#include "grid.hh"
#include "misc/planck_function.hh"
#include "rmsd.hh"
#include "wavelength_grid.hh"

void do_ALI() {

    const unsigned int n_depth_pts = grid.size();

    double rmsd;
    const double max_tol = 1.0e-8;

    unsigned int i;
#pragma omp parallel private (i)
    {
    for (i = 0; i < wavelength_values.size(); ++i) {
#pragma omp single nowait
        {
        unsigned int iter = 0;
        std::vector< std::map< unsigned int, double> >Lambda_star = calc_ALO(i);

        viennacl::vector<double> J_old(n_depth_pts);
        viennacl::vector<double> J_new(n_depth_pts);
        viennacl::vector<double> J_fs(n_depth_pts);
        for (unsigned int j = 0; j < n_depth_pts; ++j) {
            J_old(j) = grid.at(j).wavelength_grid.at(i).J;
        }
        viennacl::vector<double> rhs(n_depth_pts);
        viennacl::vector<double> epsilon(n_depth_pts);

        viennacl::vector<double> vcl_rhs(n_depth_pts);
        viennacl::compressed_matrix<double> vcl_sparsemat(n_depth_pts, n_depth_pts);

        viennacl::copy(Lambda_star, vcl_sparsemat);

        do {
            for (std::vector<Ray>::iterator r = rays.begin(); r != rays.end(); ++r) {
                for (std::vector<GridVoxel>::iterator gv = grid.begin(); gv != grid.end(); ++gv) {
                    gv->calc_source_fn(i);
                }
                r->formal_soln(i);
            }
            for (std::vector<GridVoxel>::iterator gv = grid.begin(); gv != grid.end(); ++gv) {
                gv->calc_J(i);
                gv->calc_H(i);
                gv->calc_K(i);
            }
            for (unsigned int j = 0; j < n_depth_pts; ++j) {
                J_fs(j) = grid.at(j).wavelength_grid.at(i).J;
                epsilon(j) = grid.at(j).wavelength_grid.at(i).epsilon;
            }
            rhs = viennacl::linalg::prod(Lambda_star, J_old);
            for (unsigned int j = 0; j < n_depth_pts; ++j) {
                rhs(j) *= (1.0 - epsilon(j));
            }
            rhs = J_fs - rhs;

            // WARNING: only works for diagonal ALO!
            std::vector< std::map< unsigned int, double> >mtx(n_depth_pts);

            for (unsigned int k = 0; k < n_depth_pts; ++k) {
                mtx[k][k] = 1.0 - (1.0 - epsilon(k))*(Lambda_star[k][k]);
            }

            viennacl::copy(rhs, vcl_rhs);
            viennacl::copy(mtx, vcl_sparsemat);

            // Iterate the conjugate gradient method on the ALI linear system until it converges.
            J_new = viennacl::linalg::solve(vcl_sparsemat, vcl_rhs, viennacl::linalg::cg_tag());

            rmsd = calc_rmsd(J_old, J_new);
            J_old = J_new;
            for (unsigned int j = 0; j < n_depth_pts; ++j) {
                grid.at(j).wavelength_grid.at(i).J = J_old(j);
            }
            if (rmsd < max_tol)
                break;

            if (config["print_every_iter"].as<bool>()) {
                for (std::vector<GridVoxel>::const_iterator gv = grid.begin(); gv != grid.end(); ++gv) {
                    for (std::vector<GridWavelengthPoint>::const_iterator wlp = gv->wavelength_grid.begin(); wlp != gv->wavelength_grid.end(); ++wlp) {
                        moments_file << std::setw(16) << gv->z << std::setw(15) << gv->rho << std::setw(15) << wlp->lambda << std::setw(15) << wlp->J << std::setw(15) << wlp->H << std::setw(15) << wlp->K << std::setw(15) << planck_function(wlp->lambda, gv->temperature) << std::endl;
                    }
                }
                moments_file << std::endl;
            }
            iter++;
        } while (rmsd > max_tol);
    } // #pragma omp single nowait
    } // #pragma omp parallel
    }


}
