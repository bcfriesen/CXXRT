#ifndef GRID_HH
#define GRID_HH

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
    std::vector<GridWavelengthPoint> wavelength_grid;
    void calc_J(const double lambda);
    void calc_H(const double lambda);
    void calc_K(const double lambda);
};

#endif
