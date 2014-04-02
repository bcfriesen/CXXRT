#include <algorithm>
#include <iterator>

#include "ray.hh"

bool ray_angle_sort_function(class RayData* rd1, class RayData* rd2) {
  return (rd1->mu < rd2->mu);
}

void calc_J(GridVoxel gv) {
  std::sort(gv.intersecting_raydata.begin(), gv.intersecting_raydata.end(), ray_angle_sort_function);

  double result = 0.0;
  for (auto it = gv.intersecting_raydata.begin(); it != gv.intersecting_raydata.end(); ++it) {
    const auto it_next = std::next(it, 1);
    result += 0.5 * ((*it)->I_lam + (*it_next)->I_lam) * ((*it_next)->mu - (*it)->mu);
  }

  gv.J_lam = result;
}
