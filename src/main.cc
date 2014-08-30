#include <iostream>
#include <limits>
#include <cmath>
#include <iomanip>
#include <fstream>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <yaml-cpp/yaml.h>
#include <boost/functional/hash.hpp>

#include "ALI.hh"
#include "EOS/atoms.hh"
#include "EOS/LTE_EOS.hh"
#include "grid.hh"
#include "ray.hh"
#include "globals.hh"
#include "constants.hh"
#include "planck_function.hh"
#include "initialize_rays.hh"
#include "build_internal_model.hh"
#include "read_mesa_model.hh"
#include "misc/atomic_symbols.hh"
#include "misc/calc_Delta_T.hh"

std::vector<class GridVoxel> grid;
std::vector<Ray> rays;
YAML::Node config;
std::map<std::size_t, double> wavelength_values;
std::ofstream log_file;
std::ofstream moments_file;
std::map<std::string, unsigned int> atomic_symbols;

int main(int argc, char *argv[]) {

    if (argc == 1) {
        std::cerr << "Usage: <executable name> <YAML control file>" << std::endl;
        exit(1);
    }

    config = YAML::LoadFile(argv[1]);

    const std::string log_file_name = config["log_file"].as<std::string>();
    log_file.open(log_file_name.c_str());
    log_file << std::scientific;

    log_file << "PARAMETERS USED:" << std::endl;
    log_file << config << std::endl << std::endl;

    set_up_atomic_symbols();

    if (config["read_mesa_model"].as<bool>()) {
        log_file << "Reading MESA model file: " << config["TAMS_mesa_model_name"].as<std::string>() << " ... ";
        std::flush(log_file);
        read_mesa_model(config["TAMS_mesa_model_name"].as<std::string>());
        log_file << "done." << std::endl;
    } else {
        log_file << "Building internal model ... ";
        std::flush(log_file);
        build_internal_model();
        log_file << "done." << std::endl;
    }

    // set up the wavelength grid
    const unsigned int n_wavelength_pts = config["n_wavelength_pts"].as<int>();
    const double wl_min = config["wl_min"].as<double>();
    const double wl_max = config["wl_max"].as<double>();
    if (wl_min > wl_max) {
        std::cerr << "ERROR: wl_min > wl_max!" << std::endl;
        exit(1);
    }
    // Assign a unique integer hash to each wavelength point.
    for (unsigned int i = 0; i < n_wavelength_pts; ++i) {
        const double wl_value = wl_min + double(i) * (wl_max - wl_min) / double(n_wavelength_pts-1);
        boost::hash<double> double_hash;
        const std::size_t wl_value_hash = double_hash(wl_value);
        wavelength_values[wl_value_hash] = wl_value;
    }

    for (std::vector<GridVoxel>::iterator gv = grid.begin(); gv != grid.end(); ++gv) {
        for (std::map<std::size_t, double>::const_iterator wlv = wavelength_values.begin(); wlv != wavelength_values.end(); ++wlv) {
            GridWavelengthPoint gwlp_tmp;
            gwlp_tmp.lambda = wlv->second;
            gv->wavelength_grid[wlv->first] = gwlp_tmp;
        }
    }

    std::cout << std::scientific;

    for (std::vector<GridVoxel>::iterator gv = grid.begin(); gv != grid.end(); ++gv) {
        for (std::vector<Atom>::iterator atom = gv->atoms.begin(); atom != gv->atoms.end(); ++atom) {
            for (std::vector<Ion>::iterator ion = atom->ions.begin(); ion != atom->ions.end(); ++ion) {
                const bool continuum_ion_only = ((ion->ionization_stage == atom->max_ionization_stage+1) ? true : false);
                std::cout << "Reading atomic data for Z = " << atom->atomic_number << " I = " << ion->ionization_stage << " ";
                if (continuum_ion_only) {
                    std::cout << "Setting only continuum ion." << std::endl;
                } else {
                    std::cout << "Reading atomic data." << std::endl;
                }
                ion->read_atomic_data(continuum_ion_only);
            }
        }
    }

    for (std::vector<GridVoxel>::iterator gv = grid.begin(); gv != grid.end(); ++gv) {
        for (std::vector<Atom>::iterator atom = gv->atoms.begin(); atom != gv->atoms.end(); ++atom) {
            atom->set_continuum_pointers();
        }
    }

    const std::string moments_file_name = config["moments_file"].as<std::string>();
    moments_file.open(moments_file_name.c_str());
    moments_file << std::scientific;
    moments_file << "#" << std::setw(15) << "z" << std::setw(15) << "rho" << std::setw(15) << "lambda" << std::setw(15) << "J_lam" << std::setw(15) << "H_lam" << std::setw(15) << "K_lam" << std::setw(15) << "B_lam" << std::endl;

    for (std::vector<GridVoxel>::iterator gv = grid.begin(); gv != grid.end(); ++gv) {
        for (std::vector<Atom>::iterator atom = gv->atoms.begin(); atom != gv->atoms.end(); ++atom) {
            for (std::vector<Ion>::iterator ion = atom->ions.begin(); ion != atom->ions.end(); ++ion) {
                for (std::vector<AtomicLine>::iterator line = ion->lines.begin(); line != ion->lines.end(); ++line) {
                    line->set_line_width(gv->temperature);
                }
            }
        }
    }

    log_file << std::endl;
    log_file << "Initializing rays ... ";
    std::flush(log_file);
    initialize_rays();
    log_file << "done." << std::endl;


    for (unsigned int i = 0; i < config["num_temp_correction_iterations"].as<unsigned int>(); ++i) {
        log_file << std::endl;
        log_file << "Doing temperature correction iteration " << i+1 << " ..." << std::endl;

        std::ostringstream convert_loop_index_to_string;
        convert_loop_index_to_string << i;

        for (std::vector<GridVoxel>::iterator gv = grid.begin(); gv != grid.end(); ++gv) {
            for (std::vector<Atom>::iterator atom = gv->atoms.begin(); atom != gv->atoms.end(); ++atom) {
                for (std::vector<Ion>::iterator ion = atom->ions.begin(); ion != atom->ions.end(); ++ion) {
                    ion->calc_partition_function(gv->temperature);
                }
            }
        }

        log_file << std::endl;
        log_file << "Setting matter to LTE ... ";
        std::flush(log_file);
        std::vector<GridVoxel>::iterator gv;
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
        #pragma omp parallel private (gv)
        {
            for (gv = grid.begin(); gv != grid.end(); ++gv) {
                #pragma omp single nowait
                {
                    for (std::map<std::size_t, GridWavelengthPoint>::const_iterator wlp = gv->wavelength_grid.begin(); wlp != gv->wavelength_grid.end(); ++wlp) {
                        gv->calculate_emissivity_and_opacity(wlp->first);
                    }
                }
            }
        }
        log_file << "done." << std::endl;

        for (std::vector<Ray>::iterator r = rays.begin(); r != rays.end(); ++r) {
            for (std::map<std::size_t, double>::const_iterator wlv = wavelength_values.begin(); wlv != wavelength_values.end(); ++wlv) {
                r->calc_tau(wlv->first);
                r->calc_SC_coeffs(wlv->first);
            }
        }

        do_ALI();

        if (config["print_every_iter"].as<bool>()) {
            for (std::vector<GridVoxel>::const_iterator gv = grid.begin(); gv != grid.end(); ++gv) {
                for (std::map<std::size_t, GridWavelengthPoint>::const_iterator wlp = gv->wavelength_grid.begin(); wlp != gv->wavelength_grid.end(); ++wlp) {
                    moments_file << std::setw(16) << gv->z << std::setw(15) << gv->rho << std::setw(15) << wlp->second.lambda << std::setw(15) << wlp->second.J << std::setw(15) << wlp->second.H << std::setw(15) << wlp->second.K << std::setw(15) << planck_function(wlp->second.lambda, gv->temperature) << std::endl;
                }
            }
            moments_file << std::endl;
        }

        // A bunch of post-processing integrals we need to do temperature corrections.
        log_file << std::endl << "Calculating post-processing integrals ... ";
        std::flush(log_file);
        for (std::vector<GridVoxel>::iterator gv = grid.begin(); gv != grid.end(); ++gv) {
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

        std::cout << std::endl;
        for (std::vector<GridVoxel>::iterator gv = grid.begin(); gv != grid.end(); ++gv) {
            double Delta_T = calc_Delta_T(*gv);
            if (std::abs(Delta_T / gv->temperature) > 0.1)
                // If the requested temperature change is large, damp it to at most 20% of the current temperature.
                Delta_T *= (0.2 * gv->temperature / std::abs(Delta_T));
            std::cout << "z: "<< gv->z << " current T: " << gv->temperature << " new T: " << gv->temperature + Delta_T << " Delta T: " << Delta_T <<  std::endl;
            gv->temperature += Delta_T;
        }


        // Print the emergent spectrum.
        std::ofstream spectrum_file;
        const std::string spectrum_file_name = config["spectrum_file"].as<std::string>() + "_" + convert_loop_index_to_string.str();
        spectrum_file.open(spectrum_file_name.c_str());
        spectrum_file << std::scientific;

        std::vector<GridVoxel>::const_iterator surface_gv = grid.begin();
        for (surface_gv = grid.begin(); surface_gv != grid.end()-1; ++surface_gv) {
            std::vector<GridVoxel>::const_iterator next_gv = surface_gv;
            if (next_gv->z > surface_gv->z) surface_gv = next_gv;
        }
        for (std::map<std::size_t, GridWavelengthPoint>::const_iterator gwlp = surface_gv->wavelength_grid.begin(); gwlp != surface_gv->wavelength_grid.end(); ++gwlp) {
            spectrum_file << gwlp->second.lambda * 1.0e+8 << " " << 4.0 * pi * gwlp->second.H << std::endl;
        }
        spectrum_file.close();

    }

    moments_file.close();
    log_file.close();

    return(0);
}
