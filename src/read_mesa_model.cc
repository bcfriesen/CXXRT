#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include "globals.hh"
#include "grid.hh"
#include "constants.hh"

void read_mesa_model(const std::string model_name) {

    std::ifstream model_file;
    model_file.open(model_name.c_str());
    if (!model_file.good()) {
        std::cerr << "ERROR: could not open MESA TAMS model file: " << model_name << std::endl;
        exit(1);
    }

    std::string one_line;

    double density;
    std::string density_string;
    double temperature;
    std::string temperature_string;
    double radius;
    std::string radius_string;
    unsigned int shell_index;
    std::string luminosity_string; // We currently don't use MESA's calculation of luminosity.
    std::string mass_coordinate_string; // We currently don't use MESA's calculation of mass coordinate.
    double h1_mass_frac;
    std::string h1_mass_frac_string;
    std::string he3_mass_frac_string; // We currently don't use MESA's calculation of He3 mass.
    double he4_mass_frac;
    std::string he4_mass_frac_string;
    double c12_mass_frac;
    std::string c12_mass_frac_string;
    double n14_mass_frac;
    std::string n14_mass_frac_string;
    double o16_mass_frac;
    std::string o16_mass_frac_string;
    double ne20_mass_frac;
    std::string ne20_mass_frac_string;
    double mg24_mass_frac;
    std::string mg24_mass_frac_string;

    std::string tmp;

    // skip first 4 lines, which are comments
    for (unsigned int i = 0; i < 4; ++i) {
        std::getline(model_file, one_line);
    }

    // TODO: maybe actually use the next 4 lines?
    // for now, skip the next 4 lines as well
    for (unsigned int i = 0; i < 4; ++i) {
        std::getline(model_file, one_line);
    }

    // read number of shells
    unsigned int n_shells;
    std::getline(model_file, one_line);
    std::istringstream iss(one_line);
    iss >> tmp >> n_shells;

    grid.resize(n_shells);

    for (std::vector<GridVoxel>::iterator gv = grid.begin(); gv != grid.end(); ++gv) {
        Atom H(1);
        Atom He(2);
        gv->atoms.push_back(H);
        gv->atoms.push_back(He);
    }

    // skip next 4 lines
    for (unsigned int i = 0; i < 4; ++i) {
        std::getline(model_file, one_line);
    }

    // read shell data
    for (unsigned int i = 0; i < n_shells; ++i) {
        std::getline(model_file, one_line);
        std::istringstream iss(one_line);

        /* Read to strings first because we have to replace stupid "D" with "E"
         * in Fortran's stupid scientific notation. C++ can't read "D". */
        iss >> shell_index
            >> density_string
            >> temperature_string
            >> radius_string
            >> luminosity_string
            >> mass_coordinate_string
            >> h1_mass_frac_string
            >> he3_mass_frac_string
            >> he4_mass_frac_string
            >> c12_mass_frac_string
            >> n14_mass_frac_string
            >> o16_mass_frac_string
            >> ne20_mass_frac_string
            >> mg24_mass_frac_string;

        unsigned int where_is_d(density_string.find_first_of("Dd"));
        if (where_is_d != std::string::npos)
            density_string[where_is_d] = 'E';
        density = atof(density_string.c_str()); // use std::stod when more compilers are C++11-compliant

        where_is_d = temperature_string.find_first_of("Dd");
        if (where_is_d != std::string::npos)
            temperature_string[where_is_d] = 'E';
        temperature = atof(temperature_string.c_str()); // use std::stod when more compilers are C++11-compliant

        where_is_d = radius_string.find_first_of("Dd");
        if (where_is_d != std::string::npos)
            radius_string[where_is_d] = 'E';
        radius = atof(radius_string.c_str()); // use std::stod when more compilers are C++11-compliant

        where_is_d = h1_mass_frac_string.find_first_of("Dd");
        if (where_is_d != std::string::npos)
            h1_mass_frac_string[where_is_d] = 'E';
        h1_mass_frac = atof(h1_mass_frac_string.c_str()); // use std::stod when more compilers are C++11-compliant

        where_is_d = he4_mass_frac_string.find_first_of("Dd");
        if (where_is_d != std::string::npos)
            he4_mass_frac_string[where_is_d] = 'E';
        he4_mass_frac = atof(he4_mass_frac_string.c_str()); // use std::stod when more compilers are C++11-compliant

        where_is_d = c12_mass_frac_string.find_first_of("Dd");
        if (where_is_d != std::string::npos)
            c12_mass_frac_string[where_is_d] = 'E';
        c12_mass_frac = atof(c12_mass_frac_string.c_str()); // use std::stod when more compilers are C++11-compliant

        where_is_d = n14_mass_frac_string.find_first_of("Dd");
        if (where_is_d != std::string::npos)
            n14_mass_frac_string[where_is_d] = 'E';
        n14_mass_frac = atof(n14_mass_frac_string.c_str()); // use std::stod when more compilers are C++11-compliant

        where_is_d = o16_mass_frac_string.find_first_of("Dd");
        if (where_is_d != std::string::npos)
            o16_mass_frac_string[where_is_d] = 'E';
        o16_mass_frac = atof(o16_mass_frac_string.c_str()); // use std::stod when more compilers are C++11-compliant

        where_is_d = ne20_mass_frac_string.find_first_of("Dd");
        if (where_is_d != std::string::npos)
            ne20_mass_frac_string[where_is_d] = 'E';
        ne20_mass_frac = atof(ne20_mass_frac_string.c_str()); // use std::stod when more compilers are C++11-compliant

        where_is_d = mg24_mass_frac_string.find_first_of("Dd");
        if (where_is_d != std::string::npos)
            mg24_mass_frac_string[where_is_d] = 'E';
        mg24_mass_frac = atof(mg24_mass_frac_string.c_str()); // use std::stod when more compilers are C++11-compliant

        density = std::exp(density);
        temperature = std::exp(temperature);
        radius = std::exp(radius);

        grid.at(i).rho = density;
        grid.at(i).temperature = temperature;
        grid.at(i).z = radius;

        // Make sure the mass fractions are normalized.
        const double total_mass_fraction = h1_mass_frac
                                         + he4_mass_frac
                                         + c12_mass_frac
                                         + n14_mass_frac
                                         + o16_mass_frac
                                         + ne20_mass_frac
                                         + mg24_mass_frac;
        h1_mass_frac /= total_mass_fraction;
        he4_mass_frac /= total_mass_fraction;
        c12_mass_frac /= total_mass_fraction;
        n14_mass_frac /= total_mass_fraction;
        o16_mass_frac /= total_mass_fraction;
        ne20_mass_frac /= total_mass_fraction;
        mg24_mass_frac /= total_mass_fraction;

        const double H_number_density  = grid.at(i).rho * h1_mass_frac   * N_A / H_molar_mass;
        const double He_number_density = grid.at(i).rho * he4_mass_frac  * N_A / He_molar_mass;
        const double C_number_density  = grid.at(i).rho * c12_mass_frac  * N_A / C_molar_mass;
        const double N_number_density  = grid.at(i).rho * n14_mass_frac  * N_A / N_molar_mass;
        const double O_number_density  = grid.at(i).rho * o16_mass_frac  * N_A / O_molar_mass;
        const double Ne_number_density = grid.at(i).rho * ne20_mass_frac * N_A / Ne_molar_mass;
        const double Mg_number_density = grid.at(i).rho * mg24_mass_frac * N_A / Mg_molar_mass;

        // Number density of nuclei (WITHOUT electrons).
        const double N_N = H_number_density
                         + He_number_density
                         + C_number_density
                         + N_number_density
                         + O_number_density
                         + Ne_number_density
                         + Mg_number_density;

        // Total particle density: nuclei plus electrons.
        grid.at(i).n_g = H_number_density     *  2.0
                         + He_number_density  *  3.0
                         + C_number_density   *  7.0
                         + N_number_density   *  8.0
                         + O_number_density   *  9.0
                         + Ne_number_density  * 11.0
                         + Mg_number_density  * 13.0;

        // TODO: make this less clunky; maybe put atoms in a map as well instead of a vector?
        for (std::vector<Atom>::iterator atom = grid.at(i).atoms.begin(); atom != grid.at(i).atoms.end(); ++atom) {
            if (atom->atomic_number == atomic_symbols["H"]) {
                atom->number_fraction = H_number_density / N_N;
            } else if (atom->atomic_number == atomic_symbols["He"]) {
                atom->number_fraction = He_number_density / N_N;
            } else if (atom->atomic_number == atomic_symbols["C"]) {
                atom->number_fraction = C_number_density / N_N;
            } else if (atom->atomic_number == atomic_symbols["N"]) {
                atom->number_fraction = N_number_density / N_N;
            } else if (atom->atomic_number == atomic_symbols["O"]) {
                atom->number_fraction = O_number_density / N_N;
            } else if (atom->atomic_number == atomic_symbols["Ne"]) {
                atom->number_fraction = Ne_number_density / N_N;
            } else if (atom->atomic_number == atomic_symbols["Mg"]) {
                atom->number_fraction = Mg_number_density / N_N;
            }
        }
    }

    model_file.close();
}
