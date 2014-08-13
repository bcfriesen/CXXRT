#include <iostream>
#include <limits>
#include <cmath>
#include <iomanip>
#include <fstream>
#ifdef _OPENMP
#include <omp.h>
#endif

#include <yaml-cpp/yaml.h>

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

std::vector<class GridVoxel> grid;
std::vector<Ray> rays;
YAML::Node config;
std::vector<double> wavelength_values;
std::ofstream log_file;
std::ofstream moments_file;

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


    if (config["read_mesa_model"].as<bool>()) {
        log_file << "Reading MESA model file: " << config["TAMS_mesa_model_name"].as<std::string>() << " ... ";
        read_mesa_model(config["TAMS_mesa_model_name"].as<std::string>());
        log_file << "done." << std::endl;
    } else {
        log_file << "Building internal model ... ";
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
    for (unsigned int i = 0; i < n_wavelength_pts; ++i) {
        wavelength_values.push_back(wl_min + double(i) * (wl_max - wl_min) / double(n_wavelength_pts-1));
    }

    for (GridVoxel& gv: grid) {
        for (double &wlv: wavelength_values) {
            GridWavelengthPoint gwlp_tmp;
            gwlp_tmp.lambda = &wlv;
            gv.wavelength_grid.push_back(gwlp_tmp);
        }
    }

    std::cout << std::scientific;

    initialize_rays();

    for (auto gv = grid.begin(); gv != grid.end(); ++gv) {
        // TODO: this works only for hydrogen! fix when adding more elements!!
        gv->n_g = gv->rho / H_mass;
    }

    for (auto gv = grid.begin(); gv != grid.end(); ++gv) {
        Atom H(1);
        gv->atoms.push_back(H);
    }

    for (auto gv = grid.begin(); gv != grid.end(); ++gv) {
        for (auto atom = gv->atoms.begin(); atom != gv->atoms.end(); ++atom) {
            // TODO: make this variable when we add more than 1 element
            atom->number_fraction = 1.0;
        }
    }

    // Normalize number fractions. They must add up to 1.
    for (auto gv = grid.begin(); gv != grid.end(); ++gv) {
        double tmp = 0.0;
        for (auto atom = gv->atoms.begin(); atom != gv->atoms.end(); ++atom) {
            tmp += atom->number_fraction;
        }
        for (auto atom = gv->atoms.begin(); atom != gv->atoms.end(); ++atom) {
            atom->number_fraction /= tmp;
        }
    }

    for (auto gv = grid.begin(); gv != grid.end(); ++gv) {
        for (auto atom = gv->atoms.begin(); atom != gv->atoms.end(); ++atom) {
            for (auto ion = atom->ions.begin(); ion != atom->ions.end(); ++ion) {
                ion->read_atomic_data();
            }
        }
    }

    for (auto gv = grid.begin(); gv != grid.end(); ++gv) {
        for (auto atom = gv->atoms.begin(); atom != gv->atoms.end(); ++atom) {
            atom->set_continuum_pointers();
        }
    }

    const std::string moments_file_name = config["moments_file"].as<std::string>();
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



    //---------------------END OF INITIALIZATION---------------------//

    for (auto gv = grid.begin(); gv != grid.end(); ++gv) {
        for (auto atom = gv->atoms.begin(); atom != gv->atoms.end(); ++atom) {
            for (auto ion = atom->ions.begin(); ion != atom->ions.end(); ++ion) {
                for (auto line = ion->lines.begin(); line != ion->lines.end(); ++line) {
                    line->set_line_width(gv->temperature);
                }
            }
        }
    }

    for (auto gv = grid.begin(); gv != grid.end(); ++gv) {
        for (auto atom = gv->atoms.begin(); atom != gv->atoms.end(); ++atom) {
            for (auto ion = atom->ions.begin(); ion != atom->ions.end(); ++ion) {
                ion->calc_partition_function(gv->temperature);
            }
        }
    }

    log_file << std::endl;
    log_file << "GRID VALUES:" << std::endl;
    log_file << std::setw(15) << "rho" << std::setw(15) << "temperature" << std::setw(15) << "n_e" << std::endl;
#pragma omp parallel
    {
        for (auto gv = grid.begin(); gv != grid.end(); ++gv) {
#pragma omp single nowait
            calc_n_e_LTE(*gv);
            log_file << std::setw(15) << gv->rho << std::setw(15) << gv->temperature << std::setw(15) << gv->n_e << std::endl;
        }
    }

#pragma omp parallel
    {
        for (auto gv = grid.begin(); gv != grid.end(); ++gv) {
#pragma omp single nowait
            gv->calc_LTE_populations();
            for (auto wlp = gv->wavelength_grid.begin(); wlp != gv->wavelength_grid.end(); ++wlp) {
                gv->calculate_emissivity_and_opacity(*(wlp->lambda));
            }
        }
    }

#pragma omp parallel
    {
        for (auto r = rays.begin(); r != rays.end(); ++r) {
#pragma omp single nowait
            for (auto wlv = wavelength_values.begin(); wlv != wavelength_values.end(); ++wlv) {
                r->calc_tau(*wlv);
                r->calc_SC_coeffs(*wlv);
            }
        }
    }

    do_ALI();

    // Print the emergent spectrum.
    std::ofstream spectrum_file;
    const std::string spectrum_file_name = config["spectrum_file"].as<std::string>();
    spectrum_file.open(spectrum_file_name.c_str());
    spectrum_file << std::scientific;

    auto surface_gv = grid.begin();
    for (surface_gv = grid.begin(); surface_gv != grid.end()-1; ++surface_gv) {
        auto next_gv = surface_gv;
        if (next_gv->z > surface_gv->z) surface_gv = next_gv;
    }
    for (auto gwlp = surface_gv->wavelength_grid.begin(); gwlp != surface_gv->wavelength_grid.end(); ++gwlp) {
        spectrum_file << *(gwlp->lambda) * 1.0e+8 << " " << 4.0 * pi * gwlp->H << std::endl;
    }


    spectrum_file.close();
    moments_file.close();
    log_file.close();

    return(0);
}
