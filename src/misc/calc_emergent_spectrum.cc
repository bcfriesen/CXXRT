#include "../grid.hh"
#include "../constants.hh"

double calc_emergent_spectrum(GridVoxel *gv, const unsigned int wl_index) {
    double result = 0.0;

    for (std::vector<RayIntersectionData>::const_iterator it = gv->ray_intersection_data.begin(); it != gv->ray_intersection_data.end()-1; ++it) {
        const std::vector<RayData>::iterator real_it = it->ray->raydata.begin() + it->intersection_point;

        std::vector<RayIntersectionData>::const_iterator it_next = it;
        std::advance(it_next, 1); // use std::next when more C++ compilers are C++11-compliant
        const std::vector<RayData>::iterator real_it_next = it_next->ray->raydata.begin() + it_next->intersection_point;

        // We only see outward-pointing rays.
        if (real_it->mu > 0.0 && real_it_next->mu > 0.0) {
            result += 0.5 * (real_it->mu * real_it->wavelength_grid.at(wl_index).I + real_it_next->mu * real_it_next->wavelength_grid.at(wl_index).I) * (real_it_next->mu - real_it->mu);
        }
    }
    return 4.0 * pi * (0.5 * result);
}
