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
    void read_atomic_data();
};

class Atom {
  public:
    Atom(const unsigned int atomic_number_in);
    unsigned int atomic_number;
    double atomic_weight;
    std::string atomic_symbol;
    std::vector<Ion> ions;
    double number_fraction;
    void set_continuum_pointers();
};

class AtomicLine {
  public:
    double wavelength;
    double oscillator_strength;
    class AtomicLevel* lower_level;
    class AtomicLevel* upper_level;
    double Einstein_A();
    double Einstein_B();
    double Delta_lambda; // Intrinsic line width.
    double (*line_profile) (const double lambda, const double lambda_0, const double Delta_lambda); // Function pointer to normalized line profile, e.g., Gaussian, Lorentzian, etc.
    double alpha(const double lambda); // Radiative cross-section. We assume complete redistribution so there is no need to distinguish between alpha_{i->j} and alpha_{j->i}.
    void set_line_width(const double temperature);
    double radiative_rate_absorption(const std::vector<GridWavelengthPoint> wavelength_grid);
    double radiative_rate_emission(const std::vector<GridWavelengthPoint> wavelength_grid, const double temperature);
    double collisional_rate_absorption(const double n_e, const double temperature);
};

class AtomicLevel {
  public:
    double energy;
    unsigned int g;
    double J;
    std::vector<class AtomicLine*> lines;
};

std::ostream& operator<<(std::ostream& os, const Ion& ion);

#endif
