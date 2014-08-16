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
    std::map<std::size_t, RayWavelengthPoint> wavelength_grid;
    void calc_source_fn(const std::size_t wl_value_hash);
};

class Ray {
  public:
    std::vector<class RayData> raydata;
    void bind_to_grid(const double mu);
    void calc_tau(const std::size_t wl_value_hash);
    void formal_soln(const std::size_t wl_value_hash);
    void calc_SC_coeffs(const std::size_t wl_value_hash);
    void print_ray_data(const std::size_t wl_value_hash);
};

struct RayIntersectionData {
    class Ray* ray;
    std::size_t intersection_point;
};

#endif
