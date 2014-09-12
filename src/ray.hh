#ifndef RAY_HH
#define RAY_HH

#include <list>
#include <iostream>
#include <vector>

#include "grid.hh"
#include "wavelength_grid.hh"

class RayData {
  public:
    RayData();
    class GridVoxel *gridvoxel;
    double mu;
    std::vector<RayWavelengthPoint> wavelength_grid;
};

class Ray {
  public:
    std::vector<class RayData> raydata;
    void bind_to_grid(const double mu);
    void calc_tau(const unsigned int wl_index);
    void formal_soln(const unsigned int wl_index);
    void calc_SC_coeffs(const unsigned int wl_index);
    void print_ray_data(const unsigned int wl_index);
};

struct RayIntersectionData {
    class Ray* ray;
    std::size_t intersection_point;
};

#endif
