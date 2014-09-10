#ifndef GRID_HH
#define GRID_HH

#include <map>

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
    std::vector<struct RayIntersectionData> ray_intersection_data;
    std::vector<Atom> atoms;
    std::map<std::size_t, GridWavelengthPoint> wavelength_grid;
    void calc_J(const std::size_t wl_value_hash);
    void calc_H(const std::size_t wl_value_hash);
    void calc_K(const std::size_t wl_value_hash);
    void calc_LTE_populations();
    void calculate_emissivity_and_opacity(const std::size_t wl_value_hash);
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
