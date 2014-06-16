#include "saha_equation.hh"

double f_ij (Atom atom, const double n_e, const double temperature) {
    double numerator = 1.0;
    for (unsigned int i = 0; i < atom.ions.size()-1; ++i) {
        numerator *= saha_equation(atom, i, n_e, temperature);
    }

    double denominator = 0.0;
    for (unsigned int i = 0; i < atom.ions.size()-1; ++i) {
        double tmp = 1.0;
        for (unsigned int j = 0; j <= i; ++j) {
            tmp *= saha_equation(atom, j, n_e, temperature);
        }
        denominator += tmp;
    }
    denominator += 1.0;

    return (numerator/denominator);
}
