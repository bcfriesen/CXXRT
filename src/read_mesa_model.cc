#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include "globals.hh"
#include "grid.hh"
#include "constants.hh"
#include "EOS/insert_ions.hh"

void read_mesa_model(const std::string model_name) {

    std::ifstream model_file;
    model_file.open(model_name.c_str());
    if (!model_file.good()) {
        std::cerr << "ERROR: could not open MESA TAMS model file: " << model_name << std::endl;
        exit(1);
    }

    std::string one_line;

    double density;
    double temperature;
    double radius;
    double neutron_mass_frac = 0.0;
    double h1_mass_frac = 0.0;
    double h2_mass_frac = 0.0;
    double proton_mass_frac = 0.0;
    double he3_mass_frac = 0.0;
    double he4_mass_frac = 0.0;
    double li7_mass_frac = 0.0;
    double be7_mass_frac = 0.0;
    double b8_mass_frac = 0.0;
    double c12_mass_frac = 0.0;
    double c13_mass_frac = 0.0;
    double n13_mass_frac = 0.0;
    double n14_mass_frac = 0.0;
    double n15_mass_frac = 0.0;
    double o14_mass_frac = 0.0;
    double o15_mass_frac = 0.0;
    double o16_mass_frac = 0.0;
    double o17_mass_frac = 0.0;
    double o18_mass_frac = 0.0;
    double f17_mass_frac = 0.0;
    double f18_mass_frac = 0.0;
    double f19_mass_frac = 0.0;
    double ne18_mass_frac = 0.0;
    double ne19_mass_frac = 0.0;
    double ne20_mass_frac = 0.0;
    double mg22_mass_frac = 0.0;
    double mg24_mass_frac = 0.0;
    double si28_mass_frac = 0.0;
    double s32_mass_frac = 0.0;
    double ar36_mass_frac = 0.0;
    double ca40_mass_frac = 0.0;
    double ti44_mass_frac = 0.0;
    double cr48_mass_frac = 0.0;
    double cr56_mass_frac = 0.0;
    double fe52_mass_frac = 0.0;
    double fe54_mass_frac = 0.0;
    double fe56_mass_frac = 0.0;
    double ni56_mass_frac = 0.0;

    // read number of shells
    unsigned int n_shells;
    std::getline(model_file, one_line);
    std::istringstream iss(one_line);
    iss >> n_shells;

    grid.resize(n_shells);

    for (std::vector<GridVoxel>::iterator gv = grid.begin(); gv != grid.end(); ++gv) {
        for (unsigned int i = 0; i < tmp_atoms.size(); ++i) {
            gv->atoms.push_back(tmp_atoms.at(i));
        }
    }

    // skip one line
    std::getline(model_file, one_line);

    // read shell data
    for (unsigned int i = 0; i < n_shells; ++i) {
        std::getline(model_file, one_line);
        std::istringstream iss(one_line);

        if (config["mesa_network"].as<std::string>() == "example_solar") {
            /* Read to strings first because we have to replace stupid "D" with "E"
             * in Fortran's stupid scientific notation. C++ can't read "D". */
            iss >> temperature
                >> density
                >> radius
                >> h1_mass_frac
                >> h2_mass_frac
                >> he3_mass_frac
                >> he4_mass_frac
                >> li7_mass_frac
                >> be7_mass_frac
                >> b8_mass_frac
                >> c12_mass_frac
                >> c13_mass_frac
                >> n13_mass_frac
                >> n14_mass_frac
                >> n15_mass_frac
                >> o14_mass_frac
                >> o15_mass_frac
                >> o16_mass_frac
                >> o17_mass_frac
                >> o18_mass_frac
                >> f17_mass_frac
                >> f18_mass_frac
                >> f19_mass_frac
                >> ne18_mass_frac
                >> ne19_mass_frac
                >> ne20_mass_frac
                >> mg22_mass_frac
                >> mg24_mass_frac;

            // Add isotope mass fractions together. Isotopic shifts in lines are small.
            h1_mass_frac += h2_mass_frac;
            he4_mass_frac += he3_mass_frac;
            c12_mass_frac += c13_mass_frac;
            n14_mass_frac += n13_mass_frac + n15_mass_frac;
            o16_mass_frac += o14_mass_frac + o15_mass_frac + o17_mass_frac + o18_mass_frac;
            f18_mass_frac += f17_mass_frac + f19_mass_frac;
            ne20_mass_frac += ne18_mass_frac + ne19_mass_frac;
            mg24_mass_frac += mg22_mass_frac;

        } else {
            std::cerr << "ERROR: unknown MESA nuclear network: " << config["mesa_network"].as<std::string>() << std::endl;
            exit(1);
        }

        grid.at(i).rho = density;
        grid.at(i).temperature = temperature;
        grid.at(i).z = radius * rad_sun;

        // Make sure the mass fractions are normalized.
        const double total_mass_fraction = neutron_mass_frac
                                         + h1_mass_frac
                                         + he4_mass_frac
                                         + li7_mass_frac
                                         + be7_mass_frac
                                         + b8_mass_frac
                                         + c12_mass_frac
                                         + n14_mass_frac
                                         + o16_mass_frac
                                         + f18_mass_frac
                                         + ne20_mass_frac
                                         + mg24_mass_frac
                                         + si28_mass_frac
                                         + s32_mass_frac
                                         + ar36_mass_frac
                                         + ca40_mass_frac
                                         + ti44_mass_frac
                                         + cr48_mass_frac
                                         + fe56_mass_frac
                                         + ni56_mass_frac;
        neutron_mass_frac /= total_mass_fraction;
        h1_mass_frac /= total_mass_fraction;
        he4_mass_frac /= total_mass_fraction;
        li7_mass_frac /= total_mass_fraction;
        be7_mass_frac /= total_mass_fraction;
        b8_mass_frac /= total_mass_fraction;
        c12_mass_frac /= total_mass_fraction;
        n14_mass_frac /= total_mass_fraction;
        o16_mass_frac /= total_mass_fraction;
        f18_mass_frac /= total_mass_fraction;
        ne20_mass_frac /= total_mass_fraction;
        mg24_mass_frac /= total_mass_fraction;
        si28_mass_frac /= total_mass_fraction;
        s32_mass_frac /= total_mass_fraction;
        ar36_mass_frac /= total_mass_fraction;
        ca40_mass_frac /= total_mass_fraction;
        ti44_mass_frac /= total_mass_fraction;
        cr48_mass_frac /= total_mass_fraction;
        fe56_mass_frac /= total_mass_fraction;
        ni56_mass_frac /= total_mass_fraction;

        const double H_number_density  = grid.at(i).rho * h1_mass_frac   * N_A / H_molar_mass;
        const double He_number_density = grid.at(i).rho * he4_mass_frac  * N_A / He_molar_mass;
        const double Li_number_density = grid.at(i).rho * li7_mass_frac  * N_A / Li_molar_mass;
        const double Be_number_density = grid.at(i).rho * be7_mass_frac  * N_A / Be_molar_mass;
        const double B_number_density  = grid.at(i).rho * b8_mass_frac   * N_A / B_molar_mass;
        const double C_number_density  = grid.at(i).rho * c12_mass_frac  * N_A / C_molar_mass;
        const double N_number_density  = grid.at(i).rho * n14_mass_frac  * N_A / N_molar_mass;
        const double O_number_density  = grid.at(i).rho * o16_mass_frac  * N_A / O_molar_mass;
        const double F_number_density  = grid.at(i).rho * f18_mass_frac  * N_A / F_molar_mass;
        const double Ne_number_density = grid.at(i).rho * ne20_mass_frac * N_A / Ne_molar_mass;
        const double Mg_number_density = grid.at(i).rho * mg24_mass_frac * N_A / Mg_molar_mass;
        const double Si_number_density = grid.at(i).rho * si28_mass_frac * N_A / Si_molar_mass;
        const double S_number_density  = grid.at(i).rho * s32_mass_frac  * N_A / S_molar_mass;
        const double Ar_number_density = grid.at(i).rho * ar36_mass_frac * N_A / Ar_molar_mass;
        const double Ca_number_density = grid.at(i).rho * ca40_mass_frac * N_A / Ca_molar_mass;
        const double Ti_number_density = grid.at(i).rho * ti44_mass_frac * N_A / Ti_molar_mass;
        const double Cr_number_density = grid.at(i).rho * cr48_mass_frac * N_A / Cr_molar_mass;
        const double Fe_number_density = grid.at(i).rho * fe56_mass_frac * N_A / Fe_molar_mass;
        const double Ni_number_density = grid.at(i).rho * ni56_mass_frac * N_A / Ni_molar_mass;

        // Number density of nuclei (WITHOUT electrons).
        const double N_N = H_number_density
                         + He_number_density
                         + Li_number_density
                         + Be_number_density
                         + B_number_density
                         + C_number_density
                         + N_number_density
                         + O_number_density
                         + F_number_density
                         + Ne_number_density
                         + Mg_number_density
                         + Si_number_density
                         + S_number_density
                         + Ar_number_density
                         + Ca_number_density
                         + Ti_number_density
                         + Cr_number_density
                         + Fe_number_density
                         + Ni_number_density;

        // Total particle density: nuclei plus electrons.
        grid.at(i).n_g = H_number_density     *  2.0
                         + He_number_density  *  3.0
                         + Li_number_density  *  4.0
                         + Be_number_density  *  5.0
                         + B_number_density   *  6.0
                         + C_number_density   *  7.0
                         + N_number_density   *  8.0
                         + O_number_density   *  9.0
                         + F_number_density   * 10.0
                         + Ne_number_density  * 11.0
                         + Mg_number_density  * 13.0
                         + Si_number_density  * 15.0
                         + S_number_density   * 17.0
                         + Ar_number_density  * 19.0
                         + Ca_number_density  * 21.0
                         + Ti_number_density  * 23.0
                         + Cr_number_density  * 25.0
                         + Fe_number_density  * 27.0
                         + Ni_number_density  * 29.0;

        // TODO: make this less clunky; maybe put atoms in a map as well instead of a vector?
        for (std::vector<Atom>::iterator atom = grid.at(i).atoms.begin(); atom != grid.at(i).atoms.end(); ++atom) {
            if (atom->atomic_number == atomic_symbols["H"]) {
                atom->number_fraction = H_number_density / N_N;
            } else if (atom->atomic_number == atomic_symbols["He"]) {
                atom->number_fraction = He_number_density / N_N;
            } else if (atom->atomic_number == atomic_symbols["Li"]) {
                atom->number_fraction = Li_number_density / N_N;
            } else if (atom->atomic_number == atomic_symbols["Be"]) {
                atom->number_fraction = Be_number_density / N_N;
            } else if (atom->atomic_number == atomic_symbols["B"]) {
                atom->number_fraction = B_number_density / N_N;
            } else if (atom->atomic_number == atomic_symbols["C"]) {
                atom->number_fraction = C_number_density / N_N;
            } else if (atom->atomic_number == atomic_symbols["N"]) {
                atom->number_fraction = N_number_density / N_N;
            } else if (atom->atomic_number == atomic_symbols["O"]) {
                atom->number_fraction = O_number_density / N_N;
            } else if (atom->atomic_number == atomic_symbols["F"]) {
                atom->number_fraction = F_number_density / N_N;
            } else if (atom->atomic_number == atomic_symbols["Ne"]) {
                atom->number_fraction = Ne_number_density / N_N;
            } else if (atom->atomic_number == atomic_symbols["Mg"]) {
                atom->number_fraction = Mg_number_density / N_N;
            } else if (atom->atomic_number == atomic_symbols["Si"]) {
                atom->number_fraction = Si_number_density / N_N;
            } else if (atom->atomic_number == atomic_symbols["S"]) {
                atom->number_fraction = S_number_density / N_N;
            } else if (atom->atomic_number == atomic_symbols["Ar"]) {
                atom->number_fraction = Ar_number_density / N_N;
            } else if (atom->atomic_number == atomic_symbols["Ca"]) {
                atom->number_fraction = Ca_number_density / N_N;
            } else if (atom->atomic_number == atomic_symbols["Ti"]) {
                atom->number_fraction = Ti_number_density / N_N;
            } else if (atom->atomic_number == atomic_symbols["Cr"]) {
                atom->number_fraction = Cr_number_density / N_N;
            } else if (atom->atomic_number == atomic_symbols["Fe"]) {
                atom->number_fraction = Fe_number_density / N_N;
            } else if (atom->atomic_number == atomic_symbols["Ni"]) {
                atom->number_fraction = Ni_number_density / N_N;
            }
        }
    }

    model_file.close();
}
