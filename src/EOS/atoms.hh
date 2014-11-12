#ifndef ATOMS_HH
#define ATOMS_HH

#include <vector>
#include <string>

#include "../wavelength_grid.hh"

class Ion {
  public:
    Ion (const unsigned int atomic_number_in, const unsigned int ionization_stage_in);
    unsigned int atomic_number;
    unsigned int ionization_stage;
    double atomic_weight;
    void calc_partition_function(const double temperature);
    double partition_function;
    double ionization_potential;
    std::string atomic_symbol;
    std::vector<class AtomicLevel> levels;
    std::vector<class AtomicLine> lines;
    AtomicLevel* ground_state;
    AtomicLevel* continuum_state;
    Ion* next_ion;
    void read_atomic_data(const bool continuum_ion_only);
    double alpha(const double lambda, const AtomicLevel level) const; // bound-free radiative cross-section
    double eta(const double lambda, const AtomicLevel level, const double n_e, const double temperature) const; // bound-free emissivity (recombination)
    double kappa(const double lambda, const AtomicLevel level, const double n_e, const double temperature) const; // bound-free opacity (photoionization)
    double LTE_number_density(const AtomicLevel level, const double n_e, const double temperature) const;
};

class Atom {
  public:
    Atom(const unsigned int atomic_number_in, const unsigned int max_ionization_stage_in);
    unsigned int atomic_number;
    double atomic_weight;
    std::string atomic_symbol;
    std::vector<Ion> ions;
    double number_fraction;
    void set_pointers();
    unsigned int max_ionization_stage; // Not the maximum possible, but rather the maximum that we have atomic data for.
};

class AtomicLine {
  public:
    double wavelength;
    double oscillator_strength;
    class AtomicLevel* lower_level;
    class AtomicLevel* upper_level;
    double Einstein_A() const;
    double Einstein_B_ij() const;
    double Einstein_B_ji() const;
    double Delta_lambda; // Intrinsic line width.
    double (*line_profile) (const double lambda, const double lambda_0, const double Delta_lambda); // Function pointer to normalized line profile, e.g., Gaussian, Lorentzian, etc.
    double alpha_ij(const double lambda) const;
    double alpha_ji(const double lambda) const;
    void set_line_width(const double temperature, const Atom &atom);
    double radiative_rate_absorption(const std::vector<class GridWavelengthPoint> wavelength_grid) const;
    double radiative_rate_emission(const std::vector<GridWavelengthPoint> wavelength_grid, const double temperature) const;
    double collisional_rate_absorption(const double n_e, const double temperature) const;
    double eta(const double lambda) const; // bound-bound emissivity
    double kappa(const double lambda) const; // bound-bound opacity
};

class AtomicLevel {
  public:
    double energy;
    unsigned int g;
    double J;
    std::vector<class AtomicLine*> lines;
    double number_density;
};

std::ostream& operator<<(std::ostream& os, const Ion& ion);

#endif
