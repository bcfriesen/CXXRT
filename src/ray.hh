#ifndef RAY_HH
#define RAY_HH

#include <list>
#include <iostream>
#include <vector>

#include "grid.hh"

class RayData {
  public:
    const struct GridVoxel* gridvoxel;
    double I_lam;
    double tau;
    double chi;
    double source_fn;
};

class Ray {
  public:
    std::list<struct RayData> raydata;
    double mu;
    double lambda;
    void bind_to_grid(const std::vector<struct GridVoxel> grid);
    void set_to_LTE(const double temperature);
};

#endif
