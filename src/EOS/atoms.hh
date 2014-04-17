#ifndef ATOMS_HH
#define ATOMS_HH

#include <vector>
#include <string>

class Atom {
  public:
    unsigned int atomic_number;
    double atomic_weight;
};

class Ion : public Atom {
  public:
    Ion (const unsigned int atomic_number, const unsigned int ionization_stage, const double atomic_weight);
    unsigned int ionization_stage;
    double partition_function;
    std::vector<class AtomicLevel> levels;
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
    std::string term;
    std::string configuration;
    unsigned int g;
    double J;
    std::vector<class AtomicLine*> lines;
};

std::ostream& operator<<(std::ostream& os, const Ion& ion);

#endif
