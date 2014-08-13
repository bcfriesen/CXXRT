#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include "globals.hh"
#include "grid.hh"

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
        iss >> shell_index >> density_string >> temperature_string >> radius_string;
        auto e_d(density_string.find_first_of("Dd"));
        if (e_d != std::string::npos)
            density_string[e_d] = 'E';
        density = std::stod(density_string);
        auto e_t(temperature_string.find_first_of("Dd"));
        if (e_t != std::string::npos)
            temperature_string[e_t] = 'E';
        temperature = std::stod(temperature_string);
        auto e_r(radius_string.find_first_of("Dd"));
        if (e_r != std::string::npos)
            radius_string[e_r] = 'E';
        radius = std::stod(radius_string);

        density = std::exp(density);
        temperature = std::exp(temperature);
        radius = std::exp(radius);

        grid.at(i).rho = density;
        grid.at(i).temperature = temperature;
        grid.at(i).z = radius;
    }

    model_file.close();
}
