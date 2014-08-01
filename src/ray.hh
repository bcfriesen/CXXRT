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
    const class GridVoxel* gridvoxel;
    double mu;
    std::vector<RayWavelengthPoint> wavelength_grid;
    void calc_source_fn(const double lambda);
};

class Ray {
  public:
    std::vector<class RayData> raydata;
    void bind_to_grid(const double mu);
    void calc_tau(const double lambda);
    void formal_soln(const double lambda);
    void calc_SC_coeffs(const double lambda);
};

struct RayIntersectionData {
    class Ray* ray;
    std::size_t intersection_point;
};

std::ostream& operator<<(std::ostream& os, const Ray& r);

#endif
