#include <fstream>

#include "constants.hh"
#include "globals.hh"
#include "grid.hh"

void build_internal_model() {

    const unsigned int n_depth_pts = config["n_depth_pts"].as<int>();

    if (config["n_depth_pts"].as<int>() <= 0) {
        std::cerr << "ERROR: n_depth_pts is zero or negative!" << std::endl;
        exit(1);
    }
    const double log10_rho_min = config["log10_rho_min"].as<double>();
    const double log10_rho_max = config["log10_rho_max"].as<double>();
    if (log10_rho_min > log10_rho_max) {
        std::cerr << "ERROR: log10_rho_min > log10_rho_max!" << std::endl;
        exit(1);
    }

    if (config["blackbody_temperature"].as<double>() <= 0.0) {
        std::cerr << "ERROR: black body temperature is zero or negative!" << std::endl;
        exit(1);
    }

    double log10_rho = log10_rho_min;
    const double log10_delta_rho = (log10_rho_max - log10_rho_min) / double(n_depth_pts-1);

    grid.resize(n_depth_pts);

    for (std::vector<GridVoxel>::iterator gv = grid.begin(); gv != grid.end(); ++gv) {
        gv->rho = std::pow(10.0, log10_rho);
        gv->temperature = config["blackbody_temperature"].as<double>();
        log10_rho += log10_delta_rho;
    }

    // TODO: add a switch that lets the model span the radius limits either linearly or logarithmically
    const double radius_min = config["radius_min"].as<double>();
    const double radius_max = config["radius_max"].as<double>();
    if (radius_min > radius_max) {
        std::cerr << "ERROR: radius_min > radius_max!" << std::endl;
        exit(1);
    }
    if (radius_min <= 0.0 || radius_max <= 0.0) {
        std::cerr << "ERROR: radius coordinate cannot be zero or negative!" << std::endl;
        exit(1);
    }
    unsigned int i = 0;
    for (std::vector<GridVoxel>::reverse_iterator it = grid.rbegin(); it != grid.rend(); ++it) {
        it->z = radius_min + double(i) * (radius_max - radius_min) / double(n_depth_pts-1);
        i++;
    }

    for (std::vector<GridVoxel>::iterator gv = grid.begin(); gv != grid.end(); ++gv) {
        Atom H(1, 1);
        Atom He(2, 2);
        gv->atoms.push_back(H);
        gv->atoms.push_back(He);
    }

    // Read solar abundance data.
    const std::string solar_abundance_data_file_name = config["atomic_data_root_folder"].as<std::string>() + "/" + "solar_abundances.dat";
    std::ifstream solar_abundance_data_file;
    std::string one_line;
    solar_abundance_data_file.open(solar_abundance_data_file_name.c_str()); // in C++11 the argument of open() can be a string
    if (!solar_abundance_data_file.good()) {
        std::cerr << "ERROR: could not open solar abundance data file for reading!: " << solar_abundance_data_file_name << std::endl;
        exit(1);
    }
    // Skip first line, just has comments.
    std::getline(solar_abundance_data_file, one_line);
    while (std::getline(solar_abundance_data_file, one_line)) {
        std::string element;
        double number_fraction;
        std::istringstream iss(one_line);
        iss >> element >> number_fraction;
        for (std::vector<GridVoxel>::iterator gv = grid.begin(); gv != grid.end(); ++gv) {
            for (std::vector<Atom>::iterator atom = gv->atoms.begin(); atom != gv->atoms.end(); ++atom) {
                if (atom->atomic_number == atomic_symbols[element]) {
                    atom->number_fraction = number_fraction; // These number fractions are relative to hydrogen so we have to re-normalize them after we read them all in.
                }
            }
            // Re-normalize the number fractions.
            double tmp = 0.0;
            for (std::vector<Atom>::const_iterator atom = gv->atoms.begin(); atom != gv->atoms.end(); ++atom) {
                tmp += atom->number_fraction;
            }
            for (std::vector<Atom>::iterator atom = gv->atoms.begin(); atom != gv->atoms.end(); ++atom) {
                atom->number_fraction /= tmp;
            }
        }
    }
    solar_abundance_data_file.close();

    for (std::vector<GridVoxel>::iterator gv = grid.begin(); gv != grid.end(); ++gv) {
        // TODO: get rid of this hacky stuff and make it automagic
        const double h1_mass_frac = gv->atoms.at(0).number_fraction * H_molar_mass / N_A;
        const double he4_mass_frac = gv->atoms.at(1).number_fraction * He_molar_mass / N_A;
        const double H_number_density  = gv->rho * h1_mass_frac   * N_A / H_molar_mass;
        const double He_number_density = gv->rho * he4_mass_frac  * N_A / He_molar_mass;
        grid.at(i).n_g = H_number_density     *  2.0
                         + He_number_density  *  3.0;
    }
}
