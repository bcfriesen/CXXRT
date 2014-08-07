#ifndef SAHA_EQUATION_HH
#define SAHA_EQUATION_HH

#include "atoms.hh"

double saha_equation(const Atom atom, const Ion lower_ion, const double n_e, const double temperature);

#endif
