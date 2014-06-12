#include <Eigen/Dense>

// This isn't the raw RMSD, it's the RMSD of the relative change among the
// elements of two vectors.
double calc_rmsd(const Eigen::VectorXd vec1, const Eigen::VectorXd vec2);
