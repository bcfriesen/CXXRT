#ifndef GRID_HH
#define GRID_HH

#include "ray.hh"

struct GridVoxel {
    double z;
    double rho;
    double temperature;
    std::vector<std::list<class RayData>::iterator> intersecting_raydata;
    double J_lam;
};

#endif
