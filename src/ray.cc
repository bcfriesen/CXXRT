#include <algorithm>
#include <iomanip>
#include <iterator>
#include <cmath>

#include "planck_function.hh"
#include "ray.hh"
#include "globals.hh"

void Ray::bind_to_grid(const double mu) {
  raydata.resize(grid.size());
  for (RayData& rd: raydata) {
    rd.mu = mu;
  }
  int i = 0;
  for (auto it = raydata.begin(); it != raydata.end(); ++it) {
      struct RayIntersectionData ray_intersection_data;
      if (raydata.front().mu < 0.0) {
        it->gridvoxel = &grid[i];
        ray_intersection_data.intersection_point = it - raydata.begin();
      } else {
        it->gridvoxel = &grid[(grid.size()-1)-i];
        ray_intersection_data.intersection_point = raydata.end() - it;
      }
      ray_intersection_data.ray = this;
      grid.at(i).ray_intersection_data.push_back(ray_intersection_data);
      ++i;
  }
}

void Ray::set_to_LTE() {
  // We get the temperature from the grid, so first check that the ray is bound
  // to the grid.
  if (std::none_of(raydata.begin(), raydata.end(), [](class RayData d) {return (d.gridvoxel);})) {
          std::cout << "ERROR: cannot set ray data to LTE because ray is not bound to grid!" << std::endl;
          exit(0);
  }
  for (RayData& d: raydata) {
    d.source_fn = planck_function(d.lambda, d.gridvoxel->temperature);
  }
}

std::ostream& operator<<(std::ostream& os, const Ray& r) {
  os << "RAY DATA:" << std::endl;
  os << std::setw(15) << "z";
  os << std::setw(15) << "mu";
  os << std::setw(15) << "lambda";
  os << std::setw(15) << "I_lam";
  os << std::setw(15) << "tau";
  os << std::setw(15) << "chi";
  os << std::setw(15) << "source_fn";
  os << std::endl;
  for (const RayData& rd: r.raydata) {
    os << std::setw(15) << std::scientific << rd.gridvoxel->z;
    os << std::setw(15) << std::scientific << rd.mu;
    os << std::setw(15) << std::scientific << rd.lambda;
    os << std::setw(15) << std::scientific << rd.I_lam;
    os << std::setw(15) << std::scientific << rd.tau;
    os << std::setw(15) << std::scientific << rd.chi;
    os << std::setw(15) << std::scientific << rd.source_fn;
    os << std::endl;
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
  // First check that the ray is bound to the grid.
  if (std::none_of(raydata.begin(), raydata.end(), [](class RayData d) {return (d.gridvoxel);})) {
          std::cout << "ERROR: cannot calculate chi because ray is not bound to grid!" << std::endl;
          exit(0);
  }
  for (RayData& rd: raydata) {
      // TODO: calculate opacity the right way, not by just using density as a proxy.
      rd.chi = rd.gridvoxel->rho;
  }
}

void Ray::calc_tau() {
  // First check that the ray is bound to the grid.
  if (std::none_of(raydata.begin(), raydata.end(), [](class RayData d) {return (d.gridvoxel);})) {
          std::cout << "ERROR: cannot calculate tau because ray is not bound to grid!" << std::endl;
          exit(0);
  }
  for (auto it = raydata.begin(); it != raydata.end(); ++it) {
    if (it == raydata.begin()) {
      it->tau = 0.0;
    } else {
      const auto it_prev = std::prev(it, 1);
      it->tau = it_prev->tau + (0.5 * (it_prev->chi + it->chi) * std::abs(it->gridvoxel->z - it_prev->gridvoxel->z) / std::abs(it->mu));
    }
  }
}

void Ray::calc_SC_coeffs() {
    for (auto it = raydata.begin(); it != raydata.end(); ++it) {
        if (it == raydata.begin()) {
            it->alpha = 0.0;
            it->beta = 0.0;
            it->gamma = 0.0;
        } else {
            const auto it_prev = std::prev(it, 1);
            it->Delta_tau = it->tau - it_prev->tau;
            it->alpha = 1.0 - std::exp(-it->Delta_tau) - ((it->Delta_tau - 1.0 + std::exp(-it->Delta_tau)) / it->Delta_tau);
            it->beta = (it->Delta_tau - 1.0 + std::exp(-it->Delta_tau)) / it->Delta_tau;
            // TODO: fill in gamma for parabolic interpolation
            it->gamma = 0.0;
        }
    }
}

void Ray::formal_soln() {
  for (auto it = raydata.begin(); it != raydata.end(); ++it) {
    if (it == raydata.begin()) {
      if (it->mu > 0.0) {
        it->I_lam = planck_function(it->lambda, it->gridvoxel->temperature);
      } else {
        it->I_lam = 0.0;
      }
    } else {
      const auto it_prev = std::prev(it, 1);
      const double Delta_I = (it->alpha * it_prev->source_fn) + (it->beta * it->source_fn);
      it->I_lam = it_prev->I_lam * std::exp(-it->Delta_tau) + Delta_I;
    }
  }
}
