#include "wavelength_grid.hh"
#include "planck_function.hh"


void RayWavelengthPoint::set_to_LTE(const double temperature) {
    source_fn = planck_function(lambda, temperature);
}

void RayWavelengthPoint::calc_chi(const double rho, const double lambda) {
    // TODO: calculate opacity the right way, not by just using density as a proxy.
    chi = rho;
}

