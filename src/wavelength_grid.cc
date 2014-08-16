#include "wavelength_grid.hh"
#include "planck_function.hh"


void RayWavelengthPoint::set_to_LTE(const double temperature) {
    source_fn = planck_function(lambda, temperature);
}
