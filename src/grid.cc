#include <algorithm>
#include <cmath>
#include <limits>

#include "grid.hh"

bool ray_angle_sort_function(const struct RayIntersectionData rd1, const struct RayIntersectionData rd2) {
    const auto it1 = rd1.ray->raydata.begin() + rd1.intersection_point;
    const auto it2 = rd2.ray->raydata.begin() + rd2.intersection_point;
    return (it1->mu < it2->mu);
}


void GridVoxel::calc_J(const double lambda) {
    std::sort(ray_intersection_data.begin(), ray_intersection_data.end(), ray_angle_sort_function);
    double result = 0.0;

    // Find the requested wavelength point on the grid voxel.
    // TODO: make this faster than a crude linear search.
    std::vector<GridWavelengthPoint>::iterator grid_wlp;
    for (grid_wlp = wavelength_grid.begin(); grid_wlp != wavelength_grid.end(); ++grid_wlp) {
        if (std::abs(*(grid_wlp->lambda) - lambda) < std::numeric_limits<double>::epsilon())
            break;
    }

    for (auto it = ray_intersection_data.begin(); it != ray_intersection_data.end()-1; ++it) {
        const auto real_it = it->ray->raydata.begin() + it->intersection_point;

        // Find the requested wavelength point on the rays.
        std::vector<RayWavelengthPoint>::const_iterator wlp_it;
        for (wlp_it = real_it->wavelength_grid.begin(); wlp_it != real_it->wavelength_grid.end(); ++wlp_it) {
            if (std::abs(*(wlp_it->lambda) - lambda) < std::numeric_limits<double>::epsilon())
                break;
        }

        auto it_next = it;
        std::advance(it_next, 1); // use std::next when more C++ compilers are C++11-compliant
        const auto real_it_next = it_next->ray->raydata.begin() + it_next->intersection_point;

        std::vector<RayWavelengthPoint>::const_iterator wlp_it_next;
        for (wlp_it_next = real_it_next->wavelength_grid.begin(); wlp_it_next != real_it_next->wavelength_grid.end(); ++wlp_it_next) {
            if (std::abs(*(wlp_it_next->lambda) - lambda) < std::numeric_limits<double>::epsilon())
                break;
        }

        result += 0.5 * (wlp_it->I + wlp_it_next->I) * (real_it_next->mu - real_it->mu);
    }
    grid_wlp->J = 0.5 * result;
}


void GridVoxel::calc_H(const double lambda) {
    std::sort(ray_intersection_data.begin(), ray_intersection_data.end(), ray_angle_sort_function);
    double result = 0.0;

    // Find the requested wavelength point on the grid voxel.
    // TODO: make this faster than a crude linear search.
    std::vector<GridWavelengthPoint>::iterator grid_wlp;
    for (grid_wlp = wavelength_grid.begin(); grid_wlp != wavelength_grid.end(); ++grid_wlp) {
        if (std::abs(*(grid_wlp->lambda) - lambda) < std::numeric_limits<double>::epsilon())
            break;
    }

    for (auto it = ray_intersection_data.begin(); it != ray_intersection_data.end()-1; ++it) {
        const auto real_it = it->ray->raydata.begin() + it->intersection_point;

        // Find the requested wavelength point on the rays.
        std::vector<RayWavelengthPoint>::const_iterator wlp_it;
        for (wlp_it = real_it->wavelength_grid.begin(); wlp_it != real_it->wavelength_grid.end(); ++wlp_it) {
            if (std::abs(*(wlp_it->lambda) - lambda) < std::numeric_limits<double>::epsilon())
                break;
        }

        auto it_next = it;
        std::advance(it_next, 1); // use std::next when more C++ compilers are C++11-compliant
        const auto real_it_next = it_next->ray->raydata.begin() + it_next->intersection_point;

        std::vector<RayWavelengthPoint>::const_iterator wlp_it_next;
        for (wlp_it_next = real_it_next->wavelength_grid.begin(); wlp_it_next != real_it_next->wavelength_grid.end(); ++wlp_it_next) {
            if (std::abs(*(wlp_it_next->lambda) - lambda) < std::numeric_limits<double>::epsilon())
                break;
        }

        result += 0.5 * (real_it->mu * wlp_it->I + real_it_next->mu * wlp_it_next->I) * (real_it_next->mu - real_it->mu);
    }
    grid_wlp->H = 0.5 * result;
}


void GridVoxel::calc_K(const double lambda) {
    std::sort(ray_intersection_data.begin(), ray_intersection_data.end(), ray_angle_sort_function);
    double result = 0.0;

    // Find the requested wavelength point on the grid voxel.
    // TODO: make this faster than a crude linear search.
    std::vector<GridWavelengthPoint>::iterator grid_wlp;
    for (grid_wlp = wavelength_grid.begin(); grid_wlp != wavelength_grid.end(); ++grid_wlp) {
        if (std::abs(*(grid_wlp->lambda) - lambda) < std::numeric_limits<double>::epsilon())
            break;
    }

    for (auto it = ray_intersection_data.begin(); it != ray_intersection_data.end()-1; ++it) {
        const auto real_it = it->ray->raydata.begin() + it->intersection_point;

        // Find the requested wavelength point on the rays.
        std::vector<RayWavelengthPoint>::const_iterator wlp_it;
        for (wlp_it = real_it->wavelength_grid.begin(); wlp_it != real_it->wavelength_grid.end(); ++wlp_it) {
            if (std::abs(*(wlp_it->lambda) - lambda) < std::numeric_limits<double>::epsilon())
                break;
        }

        auto it_next = it;
        std::advance(it_next, 1); // use std::next when more C++ compilers are C++11-compliant
        const auto real_it_next = it_next->ray->raydata.begin() + it_next->intersection_point;

        std::vector<RayWavelengthPoint>::const_iterator wlp_it_next;
        for (wlp_it_next = real_it_next->wavelength_grid.begin(); wlp_it_next != real_it_next->wavelength_grid.end(); ++wlp_it_next) {
            if (std::abs(*(wlp_it_next->lambda) - lambda) < std::numeric_limits<double>::epsilon())
                break;
        }

        result += 0.5 * (std::pow(real_it->mu, 2) * wlp_it->I + std::pow(real_it_next->mu, 2) * wlp_it_next->I) * (real_it_next->mu - real_it->mu);
    }
    grid_wlp->K = 0.5 * result;
}
