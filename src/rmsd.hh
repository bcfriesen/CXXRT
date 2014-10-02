#ifndef RMSD_HH
#define RMSD_HH

#include <viennacl/vector.hpp>

// This isn't the raw RMSD, it's the RMSD of the relative change among the
// elements of two vectors.
double calc_rmsd(const viennacl::vector<double> vec1, const viennacl::vector<double> vec2);

#endif
