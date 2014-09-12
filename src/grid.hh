#ifndef GRID_HH
#define GRID_HH

#include "ray.hh"
#include "wavelength_grid.hh"
#include "EOS/atoms.hh"

class GridVoxel {
  public:
    double z;
    double rho;
    double n_g; // number density of all particles
    double n_e;
    double temperature;
    double J_wl_integral;
    double H_wl_integral;
    double K_wl_integral;
    double chi_H;
    double kappa_J;
    double kappa_B;
    double Rosseland_mean_opacity;
    double Eddington_factor_f;
    double eta_minus_chi_J_wl_integral;
    double H_target;
    double sigma; // scattering opacity (independent of wavelength if we use only Thomson scattering)
    std::vector<struct RayIntersectionData> ray_intersection_data;
    std::vector<Atom> atoms;
    std::vector<GridWavelengthPoint> wavelength_grid;
    void calc_J(const unsigned int wl_index);
    void calc_H(const unsigned int wl_index);
    void calc_K(const unsigned int wl_index);
    void calc_source_fn(const unsigned int wl_index);
    void calc_LTE_populations();
    void calculate_emissivity_and_opacity(const unsigned int wl_index);
    void calc_J_wl_integral();
    void calc_H_wl_integral();
    void calc_K_wl_integral();
    void calc_chi_H();
    void calc_kappa_J();
    void calc_kappa_B();
    void calc_Eddington_factor_f();
    void calc_eta_minus_chi_J_wl_integral();
    void calc_Rosseland_mean_opacity();
};

#endif
