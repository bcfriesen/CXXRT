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
    std::string density_string;
    double temperature;
    std::string temperature_string;
    double radius;
    std::string radius_string;
    unsigned int shell_index;
    std::string luminosity_string; // We currently don't use MESA's calculation of luminosity.
    std::string mass_coordinate_string; // We currently don't use MESA's calculation of mass coordinate.
    double neutron_mass_frac = 0.0;
    std::string neutron_mass_frac_string;
    double h1_mass_frac = 0.0;
    std::string h1_mass_frac_string;
    double h2_mass_frac = 0.0;
    std::string h2_mass_frac_string;
    double proton_mass_frac = 0.0;
    std::string proton_mass_frac_string;
    double he3_mass_frac = 0.0;
    std::string he3_mass_frac_string; // We currently don't use MESA's calculation of He3 mass.
    double he4_mass_frac = 0.0;
    std::string he4_mass_frac_string;
    double li7_mass_frac = 0.0;
    std::string li7_mass_frac_string;
    double be7_mass_frac = 0.0;
    std::string be7_mass_frac_string;
    double b8_mass_frac = 0.0;
    std::string b8_mass_frac_string;
    double c12_mass_frac = 0.0;
    std::string c12_mass_frac_string;
    double c13_mass_frac = 0.0;
    std::string c13_mass_frac_string;
    double n13_mass_frac = 0.0;
    std::string n13_mass_frac_string;
    double n14_mass_frac = 0.0;
    std::string n14_mass_frac_string;
    double n15_mass_frac = 0.0;
    std::string n15_mass_frac_string;
    double o14_mass_frac = 0.0;
    std::string o14_mass_frac_string;
    double o15_mass_frac = 0.0;
    std::string o15_mass_frac_string;
    double o16_mass_frac = 0.0;
    std::string o16_mass_frac_string;
    double o17_mass_frac = 0.0;
    std::string o17_mass_frac_string;
    double o18_mass_frac = 0.0;
    std::string o18_mass_frac_string;
    double f17_mass_frac = 0.0;
    std::string f17_mass_frac_string;
    double f18_mass_frac = 0.0;
    std::string f18_mass_frac_string;
    double f19_mass_frac = 0.0;
    std::string f19_mass_frac_string;
    double ne18_mass_frac = 0.0;
    std::string ne18_mass_frac_string;
    double ne19_mass_frac = 0.0;
    std::string ne19_mass_frac_string;
    double ne20_mass_frac = 0.0;
    std::string ne20_mass_frac_string;
    double mg22_mass_frac = 0.0;
    std::string mg22_mass_frac_string;
    double mg24_mass_frac = 0.0;
    std::string mg24_mass_frac_string;
    double si28_mass_frac = 0.0;
    std::string si28_mass_frac_string;
    double s32_mass_frac = 0.0;
    std::string s32_mass_frac_string;
    double ar36_mass_frac = 0.0;
    std::string ar36_mass_frac_string;
    double ca40_mass_frac = 0.0;
    std::string ca40_mass_frac_string;
    double ti44_mass_frac = 0.0;
    std::string ti44_mass_frac_string;
    double cr48_mass_frac = 0.0;
    std::string cr48_mass_frac_string;
    double cr56_mass_frac = 0.0;
    std::string cr56_mass_frac_string;
    double fe52_mass_frac = 0.0;
    std::string fe52_mass_frac_string;
    double fe54_mass_frac = 0.0;
    std::string fe54_mass_frac_string;
    double fe56_mass_frac = 0.0;
    std::string fe56_mass_frac_string;
    double ni56_mass_frac = 0.0;
    std::string ni56_mass_frac_string;

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
        for (unsigned int i = 0; i < tmp_atoms.size(); ++i) {
            gv->atoms.push_back(tmp_atoms.at(i));
        }
    }

    // skip next 4 lines
    for (unsigned int i = 0; i < 4; ++i) {
        std::getline(model_file, one_line);
    }

    // read shell data
    for (unsigned int i = 0; i < n_shells; ++i) {
        std::getline(model_file, one_line);
        std::istringstream iss(one_line);

        // TODO: This is getting really clunky. In the future I should probably
        // parse the MESA scripts into my own format so this doesn't have to be
        // so long.
        if (config["mesa_network"].as<std::string>() == "basic") {
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

            where_is_d = he3_mass_frac_string.find_first_of("Dd");
            if (where_is_d != std::string::npos)
                he3_mass_frac_string[where_is_d] = 'E';
            he3_mass_frac = atof(he3_mass_frac_string.c_str()); // use std::stod when more compilers are C++11-compliant

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

            // Add 3He to 4He mass fraction. Isotopic shifts in lines are small.
            he4_mass_frac += he3_mass_frac;

        } else if (config["mesa_network"].as<std::string>() == "example_solar") {
            /* Read to strings first because we have to replace stupid "D" with "E"
             * in Fortran's stupid scientific notation. C++ can't read "D". */
            iss >> shell_index
                >> density_string
                >> temperature_string
                >> radius_string
                >> luminosity_string
                >> mass_coordinate_string
                >> h1_mass_frac_string
                >> h2_mass_frac_string
                >> he3_mass_frac_string
                >> he4_mass_frac_string
                >> li7_mass_frac_string
                >> be7_mass_frac_string
                >> b8_mass_frac_string
                >> c12_mass_frac_string
                >> c13_mass_frac_string
                >> n13_mass_frac_string
                >> n14_mass_frac_string
                >> n15_mass_frac_string
                >> o14_mass_frac_string
                >> o15_mass_frac_string
                >> o16_mass_frac_string
                >> o17_mass_frac_string
                >> o18_mass_frac_string
                >> f17_mass_frac_string
                >> f18_mass_frac_string
                >> f19_mass_frac_string
                >> ne18_mass_frac_string
                >> ne19_mass_frac_string
                >> ne20_mass_frac_string
                >> mg22_mass_frac_string
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

            where_is_d = h2_mass_frac_string.find_first_of("Dd");
            if (where_is_d != std::string::npos)
                h2_mass_frac_string[where_is_d] = 'E';
            h2_mass_frac = atof(h2_mass_frac_string.c_str()); // use std::stod when more compilers are C++11-compliant

            where_is_d = he3_mass_frac_string.find_first_of("Dd");
            if (where_is_d != std::string::npos)
                he3_mass_frac_string[where_is_d] = 'E';
            he3_mass_frac = atof(he3_mass_frac_string.c_str()); // use std::stod when more compilers are C++11-compliant

            where_is_d = he4_mass_frac_string.find_first_of("Dd");
            if (where_is_d != std::string::npos)
                he4_mass_frac_string[where_is_d] = 'E';
            he4_mass_frac = atof(he4_mass_frac_string.c_str()); // use std::stod when more compilers are C++11-compliant

            where_is_d = li7_mass_frac_string.find_first_of("Dd");
            if (where_is_d != std::string::npos)
                li7_mass_frac_string[where_is_d] = 'E';
            li7_mass_frac = atof(li7_mass_frac_string.c_str()); // use std::stod when more compilers are C++11-compliant

            where_is_d = be7_mass_frac_string.find_first_of("Dd");
            if (where_is_d != std::string::npos)
                be7_mass_frac_string[where_is_d] = 'E';
            be7_mass_frac = atof(be7_mass_frac_string.c_str()); // use std::stod when more compilers are C++11-compliant

            where_is_d = b8_mass_frac_string.find_first_of("Dd");
            if (where_is_d != std::string::npos)
                b8_mass_frac_string[where_is_d] = 'E';
            b8_mass_frac = atof(b8_mass_frac_string.c_str()); // use std::stod when more compilers are C++11-compliant

            where_is_d = c12_mass_frac_string.find_first_of("Dd");
            if (where_is_d != std::string::npos)
                c12_mass_frac_string[where_is_d] = 'E';
            c12_mass_frac = atof(c12_mass_frac_string.c_str()); // use std::stod when more compilers are C++11-compliant

            where_is_d = c13_mass_frac_string.find_first_of("Dd");
            if (where_is_d != std::string::npos)
                c13_mass_frac_string[where_is_d] = 'E';
            c13_mass_frac = atof(c13_mass_frac_string.c_str()); // use std::stod when more compilers are C++11-compliant

            where_is_d = n13_mass_frac_string.find_first_of("Dd");
            if (where_is_d != std::string::npos)
                n13_mass_frac_string[where_is_d] = 'E';
            n13_mass_frac = atof(n13_mass_frac_string.c_str()); // use std::stod when more compilers are C++11-compliant

            where_is_d = n14_mass_frac_string.find_first_of("Dd");
            if (where_is_d != std::string::npos)
                n14_mass_frac_string[where_is_d] = 'E';
            n14_mass_frac = atof(n14_mass_frac_string.c_str()); // use std::stod when more compilers are C++11-compliant

            where_is_d = n15_mass_frac_string.find_first_of("Dd");
            if (where_is_d != std::string::npos)
                n15_mass_frac_string[where_is_d] = 'E';
            n15_mass_frac = atof(n15_mass_frac_string.c_str()); // use std::stod when more compilers are C++11-compliant

            where_is_d = o14_mass_frac_string.find_first_of("Dd");
            if (where_is_d != std::string::npos)
                o14_mass_frac_string[where_is_d] = 'E';
            o14_mass_frac = atof(o14_mass_frac_string.c_str()); // use std::stod when more compilers are C++11-compliant

            where_is_d = o15_mass_frac_string.find_first_of("Dd");
            if (where_is_d != std::string::npos)
                o15_mass_frac_string[where_is_d] = 'E';
            o15_mass_frac = atof(o15_mass_frac_string.c_str()); // use std::stod when more compilers are C++11-compliant

            where_is_d = o16_mass_frac_string.find_first_of("Dd");
            if (where_is_d != std::string::npos)
                o16_mass_frac_string[where_is_d] = 'E';
            o16_mass_frac = atof(o16_mass_frac_string.c_str()); // use std::stod when more compilers are C++11-compliant

            where_is_d = o17_mass_frac_string.find_first_of("Dd");
            if (where_is_d != std::string::npos)
                o17_mass_frac_string[where_is_d] = 'E';
            o17_mass_frac = atof(o17_mass_frac_string.c_str()); // use std::stod when more compilers are C++11-compliant

            where_is_d = o18_mass_frac_string.find_first_of("Dd");
            if (where_is_d != std::string::npos)
                o18_mass_frac_string[where_is_d] = 'E';
            o18_mass_frac = atof(o18_mass_frac_string.c_str()); // use std::stod when more compilers are C++11-compliant

            where_is_d = f17_mass_frac_string.find_first_of("Dd");
            if (where_is_d != std::string::npos)
                f17_mass_frac_string[where_is_d] = 'E';
            f17_mass_frac = atof(f17_mass_frac_string.c_str()); // use std::stod when more compilers are C++11-compliant

            where_is_d = f18_mass_frac_string.find_first_of("Dd");
            if (where_is_d != std::string::npos)
                f18_mass_frac_string[where_is_d] = 'E';
            f18_mass_frac = atof(f18_mass_frac_string.c_str()); // use std::stod when more compilers are C++11-compliant

            where_is_d = f19_mass_frac_string.find_first_of("Dd");
            if (where_is_d != std::string::npos)
                f19_mass_frac_string[where_is_d] = 'E';
            f19_mass_frac = atof(f19_mass_frac_string.c_str()); // use std::stod when more compilers are C++11-compliant

            where_is_d = ne18_mass_frac_string.find_first_of("Dd");
            if (where_is_d != std::string::npos)
                ne18_mass_frac_string[where_is_d] = 'E';
            ne18_mass_frac = atof(ne18_mass_frac_string.c_str()); // use std::stod when more compilers are C++11-compliant

            where_is_d = ne19_mass_frac_string.find_first_of("Dd");
            if (where_is_d != std::string::npos)
                ne19_mass_frac_string[where_is_d] = 'E';
            ne19_mass_frac = atof(ne19_mass_frac_string.c_str()); // use std::stod when more compilers are C++11-compliant

            where_is_d = ne20_mass_frac_string.find_first_of("Dd");
            if (where_is_d != std::string::npos)
                ne20_mass_frac_string[where_is_d] = 'E';
            ne20_mass_frac = atof(ne20_mass_frac_string.c_str()); // use std::stod when more compilers are C++11-compliant

            where_is_d = mg22_mass_frac_string.find_first_of("Dd");
            if (where_is_d != std::string::npos)
                mg22_mass_frac_string[where_is_d] = 'E';
            mg22_mass_frac = atof(mg22_mass_frac_string.c_str()); // use std::stod when more compilers are C++11-compliant

            where_is_d = mg24_mass_frac_string.find_first_of("Dd");
            if (where_is_d != std::string::npos)
                mg24_mass_frac_string[where_is_d] = 'E';
            mg24_mass_frac = atof(mg24_mass_frac_string.c_str()); // use std::stod when more compilers are C++11-compliant

            // Add isotope mass fractions together. Isotopic shifts in lines are small.
            h1_mass_frac += h2_mass_frac;
            he4_mass_frac += he3_mass_frac;
            c12_mass_frac += c13_mass_frac;
            n14_mass_frac += n13_mass_frac + n15_mass_frac;
            o16_mass_frac += o14_mass_frac + o15_mass_frac + o17_mass_frac + o18_mass_frac;
            f18_mass_frac += f17_mass_frac + f19_mass_frac;
            ne20_mass_frac += ne18_mass_frac + ne19_mass_frac;
            mg24_mass_frac += mg22_mass_frac;

        } else if (config["mesa_network"].as<std::string>() == "approx21") {
            iss >> shell_index
                >> density_string
                >> temperature_string
                >> radius_string
                >> luminosity_string
                >> mass_coordinate_string
                >> neutron_mass_frac_string
                >> h1_mass_frac_string
                >> proton_mass_frac_string
                >> he3_mass_frac_string
                >> he4_mass_frac_string
                >> c12_mass_frac_string
                >> n14_mass_frac_string
                >> o16_mass_frac_string
                >> ne20_mass_frac_string
                >> mg24_mass_frac_string
                >> si28_mass_frac_string
                >> s32_mass_frac_string
                >> ar36_mass_frac_string
                >> ca40_mass_frac_string
                >> ti44_mass_frac_string
                >> cr48_mass_frac_string
                >> cr56_mass_frac_string
                >> fe52_mass_frac_string
                >> fe54_mass_frac_string
                >> fe56_mass_frac_string
                >> ni56_mass_frac_string;

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

            where_is_d = neutron_mass_frac_string.find_first_of("Dd");
            if (where_is_d != std::string::npos)
                neutron_mass_frac_string[where_is_d] = 'E';
            neutron_mass_frac = atof(neutron_mass_frac_string.c_str()); // use std::stod when more compilers are C++11-compliant

            where_is_d = h1_mass_frac_string.find_first_of("Dd");
            if (where_is_d != std::string::npos)
                h1_mass_frac_string[where_is_d] = 'E';
            h1_mass_frac = atof(h1_mass_frac_string.c_str()); // use std::stod when more compilers are C++11-compliant

            where_is_d = proton_mass_frac_string.find_first_of("Dd");
            if (where_is_d != std::string::npos)
                proton_mass_frac_string[where_is_d] = 'E';
            proton_mass_frac = atof(proton_mass_frac_string.c_str()); // use std::stod when more compilers are C++11-compliant

            where_is_d = he3_mass_frac_string.find_first_of("Dd");
            if (where_is_d != std::string::npos)
                he3_mass_frac_string[where_is_d] = 'E';
            he3_mass_frac = atof(he3_mass_frac_string.c_str()); // use std::stod when more compilers are C++11-compliant

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

            where_is_d = si28_mass_frac_string.find_first_of("Dd");
            if (where_is_d != std::string::npos)
                si28_mass_frac_string[where_is_d] = 'E';
            si28_mass_frac = atof(si28_mass_frac_string.c_str()); // use std::stod when more compilers are C++11-compliant

            where_is_d = s32_mass_frac_string.find_first_of("Dd");
            if (where_is_d != std::string::npos)
                s32_mass_frac_string[where_is_d] = 'E';
            s32_mass_frac = atof(s32_mass_frac_string.c_str()); // use std::stod when more compilers are C++11-compliant

            where_is_d = ar36_mass_frac_string.find_first_of("Dd");
            if (where_is_d != std::string::npos)
                ar36_mass_frac_string[where_is_d] = 'E';
            ar36_mass_frac = atof(ar36_mass_frac_string.c_str()); // use std::stod when more compilers are C++11-compliant

            where_is_d = ca40_mass_frac_string.find_first_of("Dd");
            if (where_is_d != std::string::npos)
                ca40_mass_frac_string[where_is_d] = 'E';
            ca40_mass_frac = atof(ca40_mass_frac_string.c_str()); // use std::stod when more compilers are C++11-compliant

            where_is_d = ti44_mass_frac_string.find_first_of("Dd");
            if (where_is_d != std::string::npos)
                ti44_mass_frac_string[where_is_d] = 'E';
            ti44_mass_frac = atof(ti44_mass_frac_string.c_str()); // use std::stod when more compilers are C++11-compliant

            where_is_d = cr48_mass_frac_string.find_first_of("Dd");
            if (where_is_d != std::string::npos)
                cr48_mass_frac_string[where_is_d] = 'E';
            cr48_mass_frac = atof(cr48_mass_frac_string.c_str()); // use std::stod when more compilers are C++11-compliant

            where_is_d = cr56_mass_frac_string.find_first_of("Dd");
            if (where_is_d != std::string::npos)
                cr56_mass_frac_string[where_is_d] = 'E';
            cr56_mass_frac = atof(cr56_mass_frac_string.c_str()); // use std::stod when more compilers are C++11-compliant

            where_is_d = fe52_mass_frac_string.find_first_of("Dd");
            if (where_is_d != std::string::npos)
                fe52_mass_frac_string[where_is_d] = 'E';
            fe52_mass_frac = atof(fe52_mass_frac_string.c_str()); // use std::stod when more compilers are C++11-compliant

            where_is_d = fe54_mass_frac_string.find_first_of("Dd");
            if (where_is_d != std::string::npos)
                fe54_mass_frac_string[where_is_d] = 'E';
            fe54_mass_frac = atof(fe54_mass_frac_string.c_str()); // use std::stod when more compilers are C++11-compliant

            where_is_d = fe56_mass_frac_string.find_first_of("Dd");
            if (where_is_d != std::string::npos)
                fe56_mass_frac_string[where_is_d] = 'E';
            fe56_mass_frac = atof(fe56_mass_frac_string.c_str()); // use std::stod when more compilers are C++11-compliant

            where_is_d = ni56_mass_frac_string.find_first_of("Dd");
            if (where_is_d != std::string::npos)
                ni56_mass_frac_string[where_is_d] = 'E';
            ni56_mass_frac = atof(ni56_mass_frac_string.c_str()); // use std::stod when more compilers are C++11-compliant

            // Add 3He to 4He mass fraction. Isotopic shifts in lines are small.
            he4_mass_frac += he3_mass_frac;

            // Add protons to 1H mass. A proton is just ionized hydrogen and we
            // will calculate the ionization balance much better than MESA did.
            h1_mass_frac += proton_mass_frac;

            // Add 56Cr to 48Cr mass fraction.
            cr48_mass_frac += cr56_mass_frac;

            // Add all iron isotopes to 56Fe.
            fe56_mass_frac += fe52_mass_frac + fe54_mass_frac;
        } else {
            std::cerr << "ERROR: unknown MESA nuclear network: " << config["mesa_network"].as<std::string>() << std::endl;
            exit(1);
        }

        density = std::exp(density);
        temperature = std::exp(temperature);
        radius = std::exp(radius);

        grid.at(i).rho = density;
        grid.at(i).temperature = temperature;
        grid.at(i).z = radius;

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
