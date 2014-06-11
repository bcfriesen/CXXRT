#ifndef GRID_HH
#define GRID_HH

#include "ray.hh"

class GridVoxel {
  public:
    double z;
    double rho;
    double temperature;
    std::vector<struct RayIntersectionData> ray_intersection_data;
    double J_lam;
    double H_lam; // TODO: make this a vector in less symmetric geometries
    void calc_J();
    void calc_H();
};

#endif
