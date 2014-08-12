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

    for (GridVoxel& gv: grid) {
        gv.rho = std::pow(10.0, log10_rho);
        gv.temperature = config["blackbody_temperature"].as<double>();
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
    for (auto it = grid.rbegin(); it != grid.rend(); ++it) {
        it->z = radius_min + double(i) * (radius_max - radius_min) / double(n_depth_pts-1);
        i++;
    }

}
