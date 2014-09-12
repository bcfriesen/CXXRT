#ifndef WAVELENGTH_GRID_HH
#define WAVELENGTH_GRID_HH

#include "EOS/atoms.hh"

class RayWavelengthPoint {
  public:
    double lambda; // wavelength

    double I; // specific intensity

    double tau; // optical depth

    // interpolation coefficients for short characteristic method
    double alpha;
    double beta;
    double gamma;

    double Delta_tau; // change in optical depth between this point and the previous one along a given ray
};


class GridWavelengthPoint {
  public:
    double lambda; // wavelength

    // first three moments of the radiation field
    double J;
    double H; // TODO: make this a vector in less symmetric geometries
    double K;
    double epsilon; // thermalization parameter in two-level-atom formalism

    double chi; // total opacity (absorption + scattering)
    double eta; // total emissivity
    double kappa; // total absorption opacity
    double source_fn; // source function ( = eta / chi )

    // Lines which likely contribute to opacity at this wavelength point.
    // Formally we should some over *all* lines at every wavelength points, but
    // that becomes an expensive operation when you have lots of lines and lots
    // of wavelength points. So we can set a cutoff criterion for lines
    // included in the opacity calculation at a given wavelength point, such as
    // "any line whose rest wavelength is within +/-5% of the current
    // wavelength point." Ideally we should include some actual physics in this
    // calculation, i.e., something related to the intrinsic line width.
    std::vector<class AtomicLine*> nearby_lines;
};


#endif
