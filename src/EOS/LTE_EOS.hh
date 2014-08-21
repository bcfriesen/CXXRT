#ifndef LTE_EOS_HH
#define LTE_EOS_HH

#include "../grid.hh"

double f_ij (const Atom atom, const Ion ion, const double n_e, const double temperature);

double RHS(const GridVoxel gv);

void calc_n_e_LTE(GridVoxel &gv);

#endif
