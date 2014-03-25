#ifndef RAY_HH
#define RAY_HH

#include <list>

#include "grid.hh"

class RayData {
  public:
    struct GridVoxel* gridvoxel;
    double I_lam;
    double tau;
    double chi;
    double source_fn;
};

class Ray {
  public:
    std::list<struct RayData> raydata;
    double mu;
};

#endif
