#include <algorithm>

#include "planck_function.hh"
#include "ray.hh"

#include <iomanip>

void Ray::bind_to_grid(const std::vector<struct GridVoxel> grid) {
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
  if (std::any_of(raydata.begin(), raydata.end(), [](struct RayData d) {return d.gridvoxel == nullptr;})) {
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
