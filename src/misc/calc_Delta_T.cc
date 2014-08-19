#include <cmath>

#include "../constants.hh"
#include "../grid.hh"
#include "../globals.hh"

double calc_Delta_T(GridVoxel gv) {
    double result = 0.0;

    // We need to integrate from the inside out, but the grid is (likely) ordered outside in, so use rbegin() and rend() to go in reverse order.
    for (auto gv_it = grid.rbegin(); gv_it != grid.rend()-1; ++gv_it) {
        // There's probably a more elegant way to say "do the integral until we get to the requested voxel".
        if (gv_it->z <= gv.z * (1.0 + std::numeric_limits<double>::epsilon())) {
            // Make sure we call calc_chi_H(), calc_H_target(), and calc_H_wl_integral() before this!
            const double integrand = -gv_it->chi_H * (gv_it->H_target - gv_it->H_wl_integral);
            auto gv_it_next = gv_it;
            std::advance(gv_it_next, +1);
            const double integrand_next = -gv_it_next->chi_H * (gv_it_next->H_target - gv_it_next->H_wl_integral);
            result += 0.5 * (gv_it_next->z - gv_it->z) * (integrand_next - integrand);
        }
    }

    const double B = sigma_stefan * std::pow(gv.temperature, 4) / pi;

    // Make sure we call calc_kappa_J(), calc_kappa_B(), and calc_eddington_factor_f() before this!
    result *= (gv.temperature / (4.0 * B)) * (gv.kappa_J / gv.kappa_B) / gv.Eddington_factor_f;

    return result;
}
