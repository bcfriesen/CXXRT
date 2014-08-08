#include <iostream>
#include <limits>
#include <cmath>
#include <iomanip>
#include <fstream>

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
    const std::string moments_file_name = config["moments_file"].as<std::string>();
    const std::string log_file_name = config["log_file"].as<std::string>();
    const unsigned int n_depth_pts = config["n_depth_pts"].as<int>();
    const double log10_rho_min = config["log10_rho_min"].as<double>();
    const double log10_rho_max = config["log10_rho_max"].as<double>();

    double log10_rho = log10_rho_min;
    const double log10_delta_rho = (log10_rho_max - log10_rho_min) / double(n_depth_pts-1);
    grid.resize(n_depth_pts);

    // set up the wavelength grid
    const unsigned int n_wavelength_pts = config["n_wavelength_pts"].as<int>();
    const double wl_min = config["wl_min"].as<double>();
    const double wl_max = config["wl_max"].as<double>();
    for (unsigned int i = 0; i < n_wavelength_pts; ++i) {
        wavelength_values.push_back(wl_min + double(i) * (wl_max - wl_min) / double(n_wavelength_pts-1));
    }

    log_file.open(log_file_name.c_str());
    log_file << std::scientific;

    log_file << "PARAMETERS USED:" << std::endl;
    log_file << config << std::endl << std::endl;

    for (GridVoxel& gv: grid) {
        gv.rho = std::pow(10.0, log10_rho);
        gv.temperature = config["blackbody_temperature"].as<double>();
        log10_rho += log10_delta_rho;
        for (double &wlv: wavelength_values) {
            GridWavelengthPoint gwlp_tmp;
            gwlp_tmp.lambda = &wlv;
            gv.wavelength_grid.push_back(gwlp_tmp);
        }
    }

    // TODO: add a switch that lets the model span the radius limits either linearly or logarithmically
    const double radius_min = config["radius_min"].as<double>();
    const double radius_max = config["radius_max"].as<double>();
    unsigned int i = 0;
    for (auto it = grid.rbegin(); it != grid.rend(); ++it) {
        it->z = radius_min + double(i) * (radius_max - radius_min) / double(n_depth_pts-1);
        i++;
    }

    std::cout << std::scientific;

    initialize_rays();

    for (auto &gv: grid) {
        // TODO: this works only for hydrogen! fix when adding more elements!!
        gv.n_g = gv.rho / H_mass;
    }

    for (auto &gv: grid) {
        Atom H(1);
        gv.atoms.push_back(H);
    }

    for (auto &gv: grid) {
        for (auto &atom: gv.atoms) {
            // TODO: make this variable when we add more than 1 element
            atom.number_fraction = 1.0;
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

    for (auto &gv: grid) {
        for (auto &atom: gv.atoms) {
            for (auto &ion: atom.ions) {
                ion.read_atomic_data();
            }
        }
    }

    for (auto &gv: grid) {
        for (auto &atom: gv.atoms) {
            atom.set_continuum_pointers();
        }
    }

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

    do_ALI();

    for (auto &gv: grid) {
        for (auto &atom: gv.atoms) {
            for (auto &ion: atom.ions) {
                for (auto &line: ion.lines) {
                    line.set_line_width(gv.temperature);
                }
            }
        }
    }

    for (auto &gv: grid) {
        for (auto &atom: gv.atoms) {
            for (auto &ion: atom.ions) {
                ion.calc_partition_function(gv.temperature);
            }
        }
    }

    log_file << std::endl;
    log_file << "GRID VALUES:" << std::endl;
    log_file << std::setw(15) << "rho" << std::setw(15) << "temperature" << std::setw(15) << "n_e" << std::endl;
    for (auto &gv: grid) {
        calc_n_e_LTE(gv);
        log_file << std::setw(15) << gv.rho << std::setw(15) << gv.temperature << std::setw(15) << gv.n_e << std::endl;
    }

    for (auto &gv: grid) {
        gv.calc_LTE_populations();
        for (auto wlp: gv.wavelength_grid) {
            gv.calculate_emissivity_and_opacity(*(wlp.lambda));
        }
    }

    for (auto &r: rays) {
        for (auto wlv: wavelength_values) {
            r.calc_tau(wlv);
            r.calc_SC_coeffs(wlv);
        }
    }

    do_ALI();

    // Print the emergent spectrum.
    std::ofstream spectrum_file;
    const std::string spectrum_file_name = config["spectrum_file"].as<std::string>();
    spectrum_file.open(spectrum_file_name.c_str());
    spectrum_file << std::scientific;

    for (auto gv: grid) {
        if (std::abs(gv.z - radius_max) < std::numeric_limits<double>::epsilon()) {
            for (auto gwlp: gv.wavelength_grid) {
                spectrum_file << *(gwlp.lambda) * 1.0e+8 << " " << 4.0 * pi * gwlp.H << std::endl;
            }
        }
    }


    spectrum_file.close();
    moments_file.close();
    log_file.close();

    return(0);
}
