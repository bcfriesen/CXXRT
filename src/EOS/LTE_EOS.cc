#include "saha_equation.hh"
#include "atoms.hh"
#include "../grid.hh"

double f_ij (Atom atom, const double n_e, const double temperature) {
    double numerator = 1.0;
    for (auto ion_it: atom.ions) {
        // If the requested j+1'th ionization stage is more than fully ionized,
        // then skip it. This case should happen only once: when trying to
        // calculate the Saha equation where the lower ionization stage j is
        // the fully-ionized atom.
        if (ion_it.ionization_stage+1 > atom.atomic_number) {
            continue;
        } else {
            numerator *= saha_equation(atom, ion_it.ionization_stage, n_e, temperature);
        }
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


double RHS(GridVoxel &gv) {
    double result = 0.0;
    for (auto atom: gv.atoms) {
        double tmp = 0.0;
        for (auto ion: atom.ions) {
            tmp += ion.ionization_stage * f_ij(atom, gv.n_e, gv.temperature);
        }
        tmp *= atom.number_fraction;
        result += tmp;
    }
    result = 1.0 / result;
    result = result + 1.0;
    result = result * gv.n_e;
    result = result - gv.n_g;

    return result;
}
