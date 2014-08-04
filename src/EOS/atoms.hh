#ifndef ATOMS_HH
#define ATOMS_HH

#include <vector>
#include <string>

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
};

class Atom {
  public:
    Atom(const unsigned int atomic_number_in);
    unsigned int atomic_number;
    double atomic_weight;
    std::string atomic_symbol;
    std::vector<Ion> ions;
    double number_fraction;
};

class AtomicLine {
  public:
    double wavelength;
    double oscillator_strength;
    class AtomicLevel* lower_level;
    class AtomicLevel* upper_level;
    double Einstein_B();
    double (*line_profile) (const double lambda, const double lambda_0, const double Delta_lambda); // Function pointer to normalized line profile, e.g., Gaussian, Lorentzian, etc.
    double alpha(const double lambda); // Radiative cross-section. We assume complete redistribution so there is no need to distinguish between alpha_{i->j} and alpha_{j->i}.
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
