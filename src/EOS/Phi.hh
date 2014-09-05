#ifndef PHI_HH
#define PHI_HH

#include "atoms.hh"

double Phi(const AtomicLevel &level, const Ion &ion, const Atom &atom, const double &temperature);
double Phi_tilde(const AtomicLevel &level, const Ion &ion, const Atom &atom, const double &temperature);
double Phi_tilde(const Ion &ion, const Atom &atom, const double &temperature);

#endif
