#ifndef WAVELENGTH_GRID_HH
#define WAVELENGTH_GRID_HH

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
};


#endif
