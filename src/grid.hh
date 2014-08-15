#ifndef GRID_HH
#define GRID_HH

#include <map>

#include "ray.hh"
#include "wavelength_grid.hh"
#include "EOS/atoms.hh"

class GridVoxel {
  public:
    double z;
    double rho;
    double n_g; // number density of all particles
    double n_e;
    double temperature;
    std::vector<struct RayIntersectionData> ray_intersection_data;
    std::vector<Atom> atoms;
    std::map<std::size_t, GridWavelengthPoint> wavelength_grid;
    void calc_J(const std::size_t wl_value_hash);
    void calc_H(const std::size_t wl_value_hash);
    void calc_K(const std::size_t wl_value_hash);
    void calc_LTE_populations();
    void calculate_emissivity_and_opacity(const std::size_t wl_value_hash);
};

#endif
