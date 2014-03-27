#include <algorithm>
#include <iomanip>

#include "planck_function.hh"
#include "ray.hh"
#include "globals.hh"

void Ray::bind_to_grid() {
  raydata.resize(grid.size());
  int i = 0;
  if (mu < 0.0) {
    for (auto it = raydata.begin(); it != raydata.end(); ++it) {
      it->gridvoxel = &grid[i];
      ++i;
    }
  } else {
    for (auto it = raydata.rbegin(); it != raydata.rend(); ++it) {
      it->gridvoxel = &grid[i];
      ++i;
    }
  }
}

void Ray::set_to_LTE() {
  // We get the temperature from the grid, so first check that the ray is bound
  // to the grid.
  if (std::none_of(raydata.begin(), raydata.end(), [](struct RayData d) {return (d.gridvoxel);})) {
          std::cout << "ERROR: cannot set ray data to LTE because ray is not bound to grid!" << std::endl;
          exit(0);
      }
  for (RayData& d: raydata) {
    d.source_fn = planck_function(lambda, d.gridvoxel->temperature);
  }
}

std::ostream& operator<<(std::ostream& os, const Ray& r) {
  os << "mu: " << r.mu << std::endl;
  os << "lambda: " << r.lambda << std::endl;
  os << "RAY DATA:" << std::endl;
  for (auto it = r.raydata.begin(); it != r.raydata.end(); ++it) {
    os << it->I_lam << std::setw(10) << it->tau << std::setw(10) << it->chi << std::setw(10) << it->source_fn << std::endl;
  }
  return os;
}

RayData::RayData() {
    gridvoxel = nullptr;
    I_lam = 0.0;
    tau = 0.0;
    chi = 0.0;
    source_fn = 0.0;
}

void Ray::calc_chi() {
    for (RayData& rd: raydata) {
        // TODO: calculate opacity the right way, not by just using density as a proxy.
        rd.chi = rd.gridvoxel->rho;
    }
}
