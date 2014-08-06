#ifndef PHOTOIONIZATION_CROSS_SECTION_HH
#define PHOTOIONIZATION_CROSS_SECTION_HH

#include "atoms.hh"

double photo_xs(const Ion ion,
                const double sigma_0,
                const double x,
                const double y_w,
                const double y,
                const double P,
                const double y_a);

double F_y(const double x,
           const double y_w,
           const double y,
           const double P,
           const double y_a);

#endif
