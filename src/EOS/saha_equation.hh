#ifndef SAHA_EQUATION_HH
#define SAHA_EQUATION_HH

#include "atoms.hh"

double saha_equation(Atom* atom, const unsigned int lower_ionization_stage, const double n_e, const double temperature);

#endif
