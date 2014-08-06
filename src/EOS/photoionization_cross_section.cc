#include <cmath>

#include "photoionization_cross_section.hh"

// This formula is Eq. 1 of Verner et al, The Astrophysical Journal, 465:487-498, 1996 July 1
double photo_xs(const Ion ion,
                const double sigma_0,
                const double x,
                const double y_w,
                const double y,
                const double P,
                const double y_a) {

    return (sigma_0 * F_y(x, y_w, y, P, y_a) * 1.0e-18);
}

double F_y(const double x,
           const double y_w,
           const double y,
           const double P,
           const double y_a) {
    return (std::pow(x - 1.0, 2) + std::pow(y_w, 2)) * std::pow(y, 0.5 * P - 5.5) * std::pow(1.0 + std::sqrt(y / y_a), -P);
}
