#include <algorithm>
#include <cmath>

#include "grid.hh"

bool ray_angle_sort_function(const struct RayIntersectionData rd1, const struct RayIntersectionData rd2) {
  const auto it1 = rd1.ray->raydata.begin() + rd1.intersection_point;
  const auto it2 = rd2.ray->raydata.begin() + rd2.intersection_point;
  return (it1->mu < it2->mu);
}


void GridVoxel::calc_J() {
  std::sort(ray_intersection_data.begin(), ray_intersection_data.end(), ray_angle_sort_function);
  double result = 0.0;
  for (auto it = ray_intersection_data.begin(); it != ray_intersection_data.end()-1; ++it) {
    const auto real_it = it->ray->raydata.begin() + it->intersection_point;
    auto it_next = it; std::advance(it_next, 1); // use std::next when more C++ compilers are C++11-compliant
    const auto real_it_next = it_next->ray->raydata.begin() + it_next->intersection_point;
    result += 0.5 * (real_it->I_lam + real_it_next->I_lam) * (real_it_next->mu - real_it->mu);
  }
  J_lam = 0.5 * result;
}

void GridVoxel::calc_H() {
  std::sort(ray_intersection_data.begin(), ray_intersection_data.end(), ray_angle_sort_function);
  double result = 0.0;
  for (auto it = ray_intersection_data.begin(); it != ray_intersection_data.end()-1; ++it) {
    const auto real_it = it->ray->raydata.begin() + it->intersection_point;
    auto it_next = it; std::advance(it_next, 1); // use std::next when more C++ compilers are C++11-compliant
    const auto real_it_next = it_next->ray->raydata.begin() + it_next->intersection_point;
    result += 0.5 * (real_it->mu * real_it->I_lam + real_it_next->mu * real_it_next->I_lam) * (real_it_next->mu - real_it->mu);
  }
  H_lam = 0.5 * result;
}

void GridVoxel::calc_K() {
  std::sort(ray_intersection_data.begin(), ray_intersection_data.end(), ray_angle_sort_function);
  double result = 0.0;
  for (auto it = ray_intersection_data.begin(); it != ray_intersection_data.end()-1; ++it) {
    const auto real_it = it->ray->raydata.begin() + it->intersection_point;
    auto it_next = it; std::advance(it_next, 1); // use std::next when more C++ compilers are C++11-compliant
    const auto real_it_next = it_next->ray->raydata.begin() + it_next->intersection_point;
    result += 0.5 * (std::pow(real_it->mu, 2) * real_it->I_lam + std::pow(real_it_next->mu, 2) * real_it_next->I_lam) * (real_it_next->mu - real_it->mu);
  }
  K_lam = 0.5 * result;
}
