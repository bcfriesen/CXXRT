#ifndef RAY_HH
#define RAY_HH

#include <list>
#include <iostream>
#include <vector>

#include "grid.hh"

class RayData {
  public:
    RayData();
    const struct GridVoxel* gridvoxel;
    double I_lam;
    double tau;
    double chi;
    double source_fn;
    double mu;
    double lambda;
};

class Ray {
  public:
    std::list<class RayData> raydata;
    void bind_to_grid(const double mu);
    void set_to_LTE();
    void calc_chi();
    void calc_tau();
    void formal_soln();
};

std::ostream& operator<<(std::ostream& os, const Ray& r);

#endif
