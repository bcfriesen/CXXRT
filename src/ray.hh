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
    double alpha;
    double beta;
    double gamma;
    double Delta_tau;
};

class Ray {
  public:
    std::vector<class RayData> raydata;
    void bind_to_grid(const double mu);
    void set_to_LTE();
    void calc_chi();
    void calc_tau();
    void formal_soln();
    void calc_SC_coeffs();
};

struct RayIntersectionData {
    class Ray* ray;
    std::size_t intersection_point;
};

std::ostream& operator<<(std::ostream& os, const Ray& r);

#endif
