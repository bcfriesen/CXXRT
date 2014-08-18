#include "../globals.hh"
#include "../grid.hh"

double calc_H_target(const GridVoxel gv) {
    double spatial_integral = 0.0;

    // We need to integrate from the inside out, but the grid is (likely) ordered outside in, so use rbegin() and rend() to go in reverse order.
    for (auto gv_it = grid.rbegin(); gv_it != grid.rend()-1; ++gv_it) {
        // There's probably a more elegant way to say "do the integral until we get to the requested voxel".
        if (gv_it->z <= gv.z * (1.0 + std::numeric_limits<double>::epsilon())) {
            const double wl_integral = gv_it->calc_eta_minus_chi_J_wl_integral();
            auto gv_it_next = gv_it;
            std::advance(gv_it_next, +1);
            const double wl_integral_next = gv_it_next->calc_eta_minus_chi_J_wl_integral();
            spatial_integral += 0.5 * (gv_it_next->z - gv_it->z) * (wl_integral_next - wl_integral);
        }
    }

    return spatial_integral;
}
