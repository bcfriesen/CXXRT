#include <algorithm>
#include <iterator>

#include "ray.hh"

bool ray_angle_sort_function(const struct RayIntersectionData rd1, const struct RayIntersectionData rd2) {
  const auto it1 = rd1.ray->raydata.begin() + rd1.intersection_point;
  const auto it2 = rd2.ray->raydata.begin() + rd2.intersection_point;
  return (it1->mu < it2->mu);
}

void calc_J(GridVoxel &gv) {
  std::sort(gv.ray_intersection_data.begin(), gv.ray_intersection_data.end(), ray_angle_sort_function);

  double result = 0.0;
  for (auto it = gv.ray_intersection_data.begin(); it != gv.ray_intersection_data.end()-1; ++it) {
    const auto real_it = it->ray->raydata.begin() + it->intersection_point;
    const auto it_next = std::next(real_it, 1);
    result += 0.5 * (real_it->I_lam + it_next->I_lam) * (it_next->mu - real_it->mu);
  }

  gv.J_lam = 0.5 * result;
}
