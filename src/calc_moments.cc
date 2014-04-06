#include <algorithm>
#include <iterator>

#include "ray.hh"

bool ray_angle_sort_function(struct RayIntersectionData rd1, struct RayIntersectionData rd2) {
  return (rd1.data->mu < rd2.data->mu);
}

void calc_J(GridVoxel &gv) {
  std::sort(gv.ray_intersection_data.begin(), gv.ray_intersection_data.end(), ray_angle_sort_function);

  double result = 0.0;
  for (auto it = gv.ray_intersection_data.begin(); it != gv.ray_intersection_data.end()-1; ++it) {
    const auto it_next = std::next(it, 1);
    result += 0.5 * (it->data->I_lam + it_next->data->I_lam) * (it_next->data->mu - it->data->mu);
  }

  gv.J_lam = 0.5 * result;
}
