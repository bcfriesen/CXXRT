#include "planck_function.hh"
#include "ray.hh"

void Ray::set_to_LTE(const double temperature) {
  for (RayData& d: raydata) {
    d.source_fn = planck_function(lambda, temperature);
  }
}
