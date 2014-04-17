#ifndef ATOMS_HH
#define ATOMS_HH

#include <vector>

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
    unsigned int degeneracy;
    std::vector<class AtomicLine*> lines;
};

std::ostream& operator<<(std::ostream& os, const Ion& ion);

#endif
