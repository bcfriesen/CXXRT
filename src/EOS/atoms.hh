#ifndef ATOMS_HH
#define ATOMS_HH

#include <vector>
#include <string>

class Atom {
  public:
    unsigned int atomic_number;
    double atomic_weight;
    std::string atomic_symbol;
};

class Ion : public Atom {
  public:
    Ion (const std::string atomic_symbol_in, const unsigned int atomic_number_in, const unsigned int ionization_stage_in, const double atomic_weight_in);
    unsigned int ionization_stage;
    double partition_function;
    std::vector<class AtomicLevel> levels;
    std::vector<class AtomicLine> lines;
};

class AtomicLine {
  public:
    double wavelength;
    double oscillator_strength;
    class AtomicLevel* lower_level;
    class AtomicLevel* upper_level;
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
