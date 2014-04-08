#ifndef GRID_HH
#define GRID_HH

#include "ray.hh"

struct GridVoxel {
    double z;
    double rho;
    double temperature;
    std::vector<struct RayIntersectionData> ray_intersection_data;
    double J_old;
    double J_new;
    double J_fs;
};

#endif
