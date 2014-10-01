#include <iostream>
#include <limits>
#include <cmath>
#include <iomanip>
#include <fstream>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "mpi.h"

#include <yaml-cpp/yaml.h>

#include "ALI.hh"
#include "EOS/atoms.hh"
#include "EOS/LTE_EOS.hh"
#include "grid.hh"
#include "ray.hh"
#include "globals.hh"
#include "constants.hh"
#include "misc/planck_function.hh"
#include "initialize_rays.hh"
#include "build_internal_model.hh"
#include "read_mesa_model.hh"
#include "misc/atomic_symbols.hh"
#include "misc/calc_Delta_T.hh"
#include "EOS/insert_ions.hh"
#include "misc/calc_emergent_spectrum.hh"

std::vector<class GridVoxel> grid;
std::vector<Ray> rays;
YAML::Node config;
std::vector<double> tot_wavelength_values;
std::vector<double> wavelength_values;
std::ofstream log_file;
std::ofstream moments_file;
std::map<std::string, unsigned int> atomic_symbols;
std::vector<Atom> tmp_atoms;

int main(int argc, char *argv[]) {

    if (argc == 1) {
        std::cerr << "Usage: <executable name> <YAML control file>" << std::endl;
        exit(1);
    }

    MPI::Init(argc, argv);

    const unsigned int num_mpi_proc = MPI::COMM_WORLD.Get_size();
    const int mpi_rank = MPI::COMM_WORLD.Get_rank();

    const int mpi_master_rank = 0;

    std::cout << "# of MPI processes: " << num_mpi_proc << std::endl;
    std::cout << "my MPI rank: " << mpi_rank << std::endl;

    config = YAML::LoadFile(argv[1]);

    const std::string log_file_name = config["log_file"].as<std::string>();
    log_file.open(log_file_name.c_str());
    log_file << std::scientific;

    log_file << "PARAMETERS USED:" << std::endl;
    log_file << config << std::endl << std::endl;

    set_up_atomic_symbols();

    insert_ions(tmp_atoms);

    if (config["read_mesa_model"].as<bool>()) {
        log_file << "Reading MESA model file: " << config["TAMS_mesa_model_name"].as<std::string>() << " ... ";
        std::flush(log_file);
        read_mesa_model(config["TAMS_mesa_model_name"].as<std::string>());
        log_file << "done." << std::endl << std::endl;
    } else {
        log_file << "Building internal model ... ";
        std::flush(log_file);
        build_internal_model();
        log_file << "done." << std::endl << std::endl;;
    }

    // set up the wavelength grid
    log_file << "Building wavelength grid ... ";
    std::flush(log_file);
    const unsigned int n_wavelength_pts = config["n_wavelength_pts"].as<int>();
    const double wl_min = config["wl_min"].as<double>();
    const double wl_max = config["wl_max"].as<double>();
    if (wl_min > wl_max) {
        std::cerr << "ERROR: wl_min > wl_max!" << std::endl;
        exit(1);
    }

    for (unsigned int i = 0; i < n_wavelength_pts; ++i) {
        tot_wavelength_values.push_back(wl_min + double(i) * (wl_max - wl_min) / double(n_wavelength_pts-1));
    }
    wavelength_values.resize(n_wavelength_pts/num_mpi_proc);

    MPI_Scatter(&tot_wavelength_values.front(), n_wavelength_pts/num_mpi_proc, MPI_DOUBLE,
                &wavelength_values.front(),     n_wavelength_pts/num_mpi_proc, MPI_DOUBLE,
                mpi_master_rank, MPI::COMM_WORLD);

    for (std::vector<GridVoxel>::iterator gv = grid.begin(); gv != grid.end(); ++gv) {
        gv->wavelength_grid.resize(wavelength_values.size());
        for (unsigned int i = 0; i < wavelength_values.size(); ++i) {
            GridWavelengthPoint gwlp_tmp;
            gwlp_tmp.lambda = wavelength_values.at(i);
            gv->wavelength_grid.at(i) = gwlp_tmp;
        }
    }
    log_file << "done." << std::endl << std::endl;

    std::cout << std::scientific;

    log_file << "Setting model atom pointers ... ";
    std::flush(log_file);
    std::vector<GridVoxel>::iterator gv;
    #pragma omp parallel private (gv)
    {
        for (gv = grid.begin(); gv != grid.end(); ++gv) {
            #pragma omp single nowait
            {
                for (std::vector<Atom>::iterator atom = gv->atoms.begin(); atom != gv->atoms.end(); ++atom) {
                    atom->set_pointers();
                }
            }
        }
    }
    log_file << "done." << std::endl << std::endl;

    const std::string moments_file_name = config["moments_file"].as<std::string>();
    moments_file.open(moments_file_name.c_str());
    moments_file << std::scientific;
    moments_file << "#" << std::setw(15) << "z" << std::setw(15) << "rho" << std::setw(15) << "lambda" << std::setw(15) << "J_lam" << std::setw(15) << "H_lam" << std::setw(15) << "K_lam" << std::setw(15) << "B_lam" << std::endl;

    for (std::vector<GridVoxel>::iterator gv = grid.begin(); gv != grid.end(); ++gv) {
        for (std::vector<Atom>::iterator atom = gv->atoms.begin(); atom != gv->atoms.end(); ++atom) {
            for (std::vector<Ion>::iterator ion = atom->ions.begin(); ion != atom->ions.end(); ++ion) {
                for (std::vector<AtomicLine>::iterator line = ion->lines.begin(); line != ion->lines.end(); ++line) {
                    line->set_line_width(gv->temperature, *atom);
                }
            }
        }
    }

    log_file << "Initializing rays ... ";
    std::flush(log_file);
    initialize_rays();
    log_file << "done." << std::endl << std::endl;

    log_file << "Finding nearby lines each wavelength point ... ";
    std::flush(log_file);
    #pragma omp parallel private (gv)
    {
        for (gv = grid.begin(); gv != grid.end(); ++gv) {
            #pragma omp single nowait
            {
                for (unsigned int i = 0; i < wavelength_values.size(); ++i) {
                    gv->find_nearby_lines(i);
                }
            }
        }
    }
    log_file << "done." << std::endl;

    for (unsigned int i = 0; i < config["num_temp_correction_iterations"].as<unsigned int>(); ++i) {
        log_file << std::endl;
        log_file << "Doing temperature correction iteration " << i+1 << " ..." << std::endl;

        std::ostringstream convert_loop_index_to_string;
        convert_loop_index_to_string << i;

        for (gv = grid.begin(); gv != grid.end(); ++gv) {
            for (std::vector<Atom>::iterator atom = gv->atoms.begin(); atom != gv->atoms.end(); ++atom) {
                for (std::vector<Ion>::iterator ion = atom->ions.begin(); ion != atom->ions.end(); ++ion) {
                    ion->calc_partition_function(gv->temperature);
                }
            }
        }

        log_file << std::endl;
        log_file << "Setting matter to LTE ... ";
        std::flush(log_file);
        #pragma omp parallel private (gv)
        {
            for (gv = grid.begin(); gv != grid.end(); ++gv) {
                #pragma omp single nowait
                {
                    calc_n_e_LTE(*gv);
                    gv->calc_LTE_populations();
                }
            }
        }
        log_file << "done." << std::endl;

        log_file << std::endl;
        log_file << "Calculating opacities ... ";
        std::flush(log_file);
        for (gv = grid.begin(); gv != grid.end(); ++gv) {
            unsigned int i;
            #pragma omp parallel private (i)
            {
                for (i = 0; i < gv->wavelength_grid.size(); ++i) {
                    #pragma omp single nowait
                    {
                        gv->calculate_emissivity_and_opacity(i);
                    }
                }
            }
            gv->calc_Rosseland_mean_opacity();
        }
        log_file << "done." << std::endl;

        for (std::vector<Ray>::iterator r = rays.begin(); r != rays.end(); ++r) {
            for (unsigned int i = 0; i < wavelength_values.size(); ++i) {
                r->calc_tau(i);
                r->calc_SC_coeffs(i);
            }
        }

        do_ALI();

        if (config["print_every_iter"].as<bool>()) {
            for (std::vector<GridVoxel>::const_iterator gv = grid.begin(); gv != grid.end(); ++gv) {
                for (std::vector<GridWavelengthPoint>::const_iterator wlp = gv->wavelength_grid.begin(); wlp != gv->wavelength_grid.end(); ++wlp) {
                    moments_file << std::setw(16) << gv->z << std::setw(15) << gv->rho << std::setw(15) << wlp->lambda << std::setw(15) << wlp->J << std::setw(15) << wlp->H << std::setw(15) << wlp->K << std::setw(15) << planck_function(wlp->lambda, gv->temperature) << std::endl;
                }
            }
            moments_file << std::endl;
        }

        // A bunch of post-processing integrals we need to do temperature corrections.
        log_file << std::endl << "Calculating post-processing integrals ... ";
        std::flush(log_file);
        for (gv = grid.begin(); gv != grid.end(); ++gv) {
            gv->calc_J_wl_integral();
            gv->calc_H_wl_integral();
            gv->calc_K_wl_integral();
            gv->calc_chi_H();
            gv->calc_kappa_J();
            gv->calc_kappa_B();
            gv->calc_Eddington_factor_f();
            gv->calc_eta_minus_chi_J_wl_integral();
        }
        log_file << "done." << std::endl;

        // Set up temporary vectors to store wavelength-integrated data which
        // needs to be sent to MPI master rank for reduction and processing.
        // This lets me be lazy and not construct MPI data types to do the
        // communication since these data lie within complicated GridVoxel
        // classes.
        std::vector<double> J_wl_integral(grid.size());
        std::vector<double> H_wl_integral(grid.size());
        std::vector<double> K_wl_integral(grid.size());
        std::vector<double> chi_H(grid.size());
        std::vector<double> kappa_J(grid.size());
        std::vector<double> kappa_B(grid.size());
        std::vector<double> Eddington_factor_f(grid.size());
        std::vector<double> eta_minus_chi_J_wl_integral(grid.size());

        // These will store the MPI_Reduce()-ed values, which will just be
        // MPI_SUM-ed. Only the master rank needs to worry about these.
        std::vector<double> J_wl_integral_tot(grid.size());
        std::vector<double> H_wl_integral_tot(grid.size());
        std::vector<double> K_wl_integral_tot(grid.size());
        std::vector<double> chi_H_tot(grid.size());
        std::vector<double> kappa_J_tot(grid.size());
        std::vector<double> kappa_B_tot(grid.size());
        std::vector<double> Eddington_factor_f_tot(grid.size());
        std::vector<double> eta_minus_chi_J_wl_integral_tot(grid.size());

        for (unsigned int j = 0; j < grid.size(); ++j) {
            J_wl_integral.at(j) = grid.at(j).J_wl_integral;
            H_wl_integral.at(j) = grid.at(j).H_wl_integral;
            K_wl_integral.at(j) = grid.at(j).K_wl_integral;
            chi_H.at(j) = grid.at(j).chi_H;
            kappa_J.at(j) = grid.at(j).kappa_J;
            kappa_B.at(j) = grid.at(j).kappa_B;
            Eddington_factor_f.at(j) = grid.at(j).Eddington_factor_f;
            eta_minus_chi_J_wl_integral.at(j) = grid.at(j).eta_minus_chi_J_wl_integral;
        }
        // Because all of these data are just integrals, it's perfectly valid
        // to have each MPI task do the integral over only its own section of
        // the wavelength grid, and then add up the results.
        MPI_Reduce(&J_wl_integral.front(), &J_wl_integral_tot.front(), grid.size(), MPI_DOUBLE, MPI_SUM, mpi_master_rank, MPI::COMM_WORLD);
        MPI_Reduce(&H_wl_integral.front(), &H_wl_integral_tot.front(), grid.size(), MPI_DOUBLE, MPI_SUM, mpi_master_rank, MPI::COMM_WORLD);
        MPI_Reduce(&K_wl_integral.front(), &K_wl_integral_tot.front(), grid.size(), MPI_DOUBLE, MPI_SUM, mpi_master_rank, MPI::COMM_WORLD);
        MPI_Reduce(&chi_H.front(), &chi_H_tot.front(), grid.size(), MPI_DOUBLE, MPI_SUM, mpi_master_rank, MPI::COMM_WORLD);
        MPI_Reduce(&kappa_J.front(), &kappa_J_tot.front(), grid.size(), MPI_DOUBLE, MPI_SUM, mpi_master_rank, MPI::COMM_WORLD);
        MPI_Reduce(&kappa_B.front(), &kappa_B_tot.front(), grid.size(), MPI_DOUBLE, MPI_SUM, mpi_master_rank, MPI::COMM_WORLD);
        MPI_Reduce(&Eddington_factor_f.front(), &Eddington_factor_f_tot.front(), grid.size(), MPI_DOUBLE, MPI_SUM, mpi_master_rank, MPI::COMM_WORLD);
        MPI_Reduce(&eta_minus_chi_J_wl_integral.front(), &eta_minus_chi_J_wl_integral_tot.front(), grid.size(), MPI_DOUBLE, MPI_SUM, mpi_master_rank, MPI::COMM_WORLD);

        // Only the master rank will compute Delta_T, then will broadcast it to everyone else.
        std::vector<double> Delta_T(grid.size());
        if (mpi_rank == mpi_master_rank) {
            for (unsigned int j = 0; j < grid.size(); ++j) {
                grid.at(j).J_wl_integral = J_wl_integral_tot.at(j);
                grid.at(j).H_wl_integral = H_wl_integral_tot.at(j);
                grid.at(j).K_wl_integral = K_wl_integral_tot.at(j);
                grid.at(j).chi_H = chi_H_tot.at(j);
                grid.at(j).kappa_J = kappa_J_tot.at(j);
                grid.at(j).kappa_B = kappa_B_tot.at(j);
                grid.at(j).Eddington_factor_f = Eddington_factor_f_tot.at(j);
                grid.at(j).eta_minus_chi_J_wl_integral = eta_minus_chi_J_wl_integral_tot.at(j);
            }

            if (config["do_temperature_corrections"].as<bool>()) {
                std::cout << std::endl;
                for (unsigned int j = 0; j < grid.size(); ++j) {
                    Delta_T.at(j) = calc_Delta_T(grid.at(j));
                }
            }
        }
        // Now broadcast the computed temperature correction to everybody, and they will apply it themselves.
        MPI_Bcast(&Delta_T.front(), grid.size(), MPI_DOUBLE, mpi_master_rank, MPI::COMM_WORLD);

        for (unsigned int j = 0; j < grid.size(); ++j) {
            if (std::abs(Delta_T.at(j) / grid.at(j).temperature) > 0.1)
                // If the requested temperature change is large, damp it to at most 20% of the current temperature.
                Delta_T.at(j) *= (0.2 * grid.at(j).temperature / std::abs(Delta_T.at(j)));
            if (mpi_rank == mpi_master_rank) {
                std::cout << "z: "<< grid.at(j).z << " current T: " << grid.at(j).temperature << " new T: " << grid.at(j).temperature + Delta_T.at(j) << " Delta T: " << Delta_T.at(j) << std::endl;
            }
            grid.at(j).temperature += Delta_T.at(j);
        }


        // Print the emergent spectrum.
        std::ofstream spectrum_file;
        const std::string spectrum_file_name = config["spectrum_file"].as<std::string>() + "_" + convert_loop_index_to_string.str();
        spectrum_file.open(spectrum_file_name.c_str());
        spectrum_file << std::scientific;

        // Find the surface layer.
        std::vector<GridVoxel>::iterator surface_gv = grid.begin();
        for (surface_gv = grid.begin(); surface_gv != grid.end()-1; ++surface_gv) {
            std::vector<GridVoxel>::iterator next_gv = surface_gv;
            if (next_gv->z > surface_gv->z) surface_gv = next_gv;
        }

        std::vector<double> emergent_spectrum(wavelength_values.size());
        #pragma omp parallel private (i)
        {
            for (i = 0; i < wavelength_values.size(); ++i) {
                #pragma omp single nowait
                {
                    emergent_spectrum.at(i) = calc_emergent_spectrum(&(*surface_gv), i);
                }
            }
        }
        for (unsigned int i = 0; i < wavelength_values.size(); ++i) {
            spectrum_file << wavelength_values.at(i) * 1.0e+8 << " " << emergent_spectrum.at(i) << std::endl;
        }
        spectrum_file.close();

    }

    for (std::vector<Ray>::iterator ray = rays.begin(); ray != rays.end(); ++ray) {
        ray->calc_tau_Rosseland();
        ray->print_ray_data(wavelength_values.size()/2);
    }

    moments_file.close();
    log_file.close();

    MPI::Finalize();

    return(0);
}
