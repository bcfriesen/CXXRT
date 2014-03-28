#ifndef GRID_HH
#define GRID_HH

#include "ray.hh"

struct GridVoxel {
    double z;
    double rho;
    double temperature;
    std::vector<class Ray*> intersecting_rays;
};

#endif
