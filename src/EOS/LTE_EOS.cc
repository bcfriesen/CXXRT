#include <algorithm>
#include <limits>
#include <cmath>

#include "saha_equation.hh"
#include "atoms.hh"
#include "../grid.hh"
#include "Phi.hh"

// The fraction of species k in ionization stage j relative to the total number
// of atoms of that species. (See Mihalas p 114). This is Mihalas Eq (5-17).
double f_ij (const Atom atom, const Ion ion, const double n_e, const double temperature) {
    double numerator = 1.0;
    for (std::vector<Ion>::const_iterator ion_it = atom.ions.begin(); ion_it != atom.ions.end(); ++ion_it) {
        // If the requested j+1'th ionization stage is more than fully ionized,
        // then skip it. This case should happen only once: when trying to
        // calculate the Saha equation where the lower ionization stage j is
        // the fully-ionized atom.
        if (ion_it->ionization_stage+1 > atom.atomic_number) {
            continue;
        } else if (ion_it->ionization_stage >= ion.ionization_stage) {
            numerator *= (n_e * Phi_tilde(*ion_it, atom, temperature));
        }
    }

    double denominator = 0.0;
    for (std::vector<Ion>::const_iterator ion_it = atom.ions.begin(); ion_it != atom.ions.end(); ++ion_it) {
        double tmp = 1.0;
        for (std::vector<Ion>::const_iterator ion_it2 = atom.ions.begin(); ion_it2 != atom.ions.end(); ++ion_it2) {
            if (ion_it2->ionization_stage >= ion_it->ionization_stage && ion_it2->ionization_stage < ion.atomic_number) {
                tmp *= (n_e * Phi_tilde(*ion_it2, atom, temperature));
            }
        }
        denominator += tmp;
    }
    return (numerator/denominator);
}


double RHS(const GridVoxel gv) {
    double result = 0.0;
    for (std::vector<Atom>::const_iterator atom = gv.atoms.begin(); atom != gv.atoms.end(); ++atom) {
        double tmp = 0.0;
        for (std::vector<Ion>::const_iterator ion = atom->ions.begin(); ion != atom->ions.end(); ++ion) {
            tmp += ion->ionization_stage * f_ij(*atom, *ion, gv.n_e, gv.temperature);
        }
        tmp *= atom->number_fraction;
        result += tmp;
    }
    result *= (gv.n_g - gv.n_e);
    result = gv.n_e - result;

    return result;
}


void calc_n_e_LTE(GridVoxel &gv) {
    // Set lower and upper limits for n_e.
    double max_n_e = (1.0 - std::numeric_limits<double>::epsilon()) * gv.n_g;
    double min_n_e = std::numeric_limits<double>::epsilon() * gv.n_g;
    double root;

    const double max_iter = 100;
    const double tol = 1.0e-9;

    for (unsigned int i = 0; i < max_iter; ++i) {
        gv.n_e = 0.5 * (max_n_e + min_n_e);
        root = RHS(gv);
        if (root > 0.0) {
            max_n_e = gv.n_e;
        } else {
            min_n_e = gv.n_e;
        }
        if (std::abs(root)/gv.n_g < tol) return;
    }
    std::cerr << "EOS: could not converge n_e in " << max_iter << " iterations! Quitting ..." << std::endl;
    exit(1);
}
