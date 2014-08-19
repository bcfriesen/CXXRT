#include "../globals.hh"
#include "../grid.hh"

void calc_H_target(GridVoxel &gv) {
    gv.H_target = 0.0;

    // We need to integrate from the inside out, but the grid is (likely) ordered outside in, so use rbegin() and rend() to go in reverse order.
    for (auto gv_it = grid.rbegin(); gv_it != grid.rend()-1; ++gv_it) {
        // There's probably a more elegant way to say "do the integral until we get to the requested voxel".
        if (gv_it->z <= gv.z * (1.0 + std::numeric_limits<double>::epsilon())) {
            // Make sure we calculate calc_eta_minus_chi_j_wl_integral() before this!
            const double wl_integral = gv_it->eta_minus_chi_J_wl_integral;
            auto gv_it_next = gv_it;
            std::advance(gv_it_next, +1);
            // Make sure we calculate calc_eta_minus_chi_j_wl_integral() before this!
            const double wl_integral_next = gv_it_next->eta_minus_chi_J_wl_integral;
            gv.H_target += 0.5 * (gv_it_next->z - gv_it->z) * (wl_integral_next - wl_integral);
        }
    }
}
