#ifndef LTE_EOS_HH
#define LTE_EOS_HH

#include "../grid.hh"

// the fraction of the population of atom i in ionization stage j
double f_ij (const Atom atom, const double n_e, const double temperature);

double RHS(GridVoxel &gv);

#endif
