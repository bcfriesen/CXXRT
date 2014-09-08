#include <fstream>

#include "../globals.hh"
#include "../grid.hh"

void insert_ions(std::vector<Atom> &atoms) {
    std::vector<unsigned int> ions;
    for (YAML::const_iterator it = config["atoms_with_max_ion"].begin(); it != config["atoms_with_max_ion"].end(); ++it) {
        // Make sure we don't try to add an atom which is more than fully ionized.
        if (it->second.as<unsigned int>() >= atomic_symbols[it->first.as<std::string>()]) {
            std::cerr << "ERROR: you requested atom " << it->first.as<std::string>() << " with ionization stage " << it->second.as<unsigned int>() << std::endl;
            exit(1);
        }
        Atom one_atom(atomic_symbols[it->first.as<std::string>()], it->second.as<unsigned int>());
        atoms.push_back(one_atom);
    }

    log_file << "Reading atomic data ... ";
    std::flush(log_file);
    for (std::vector<Atom>::iterator atom = atoms.begin(); atom != atoms.end(); ++atom) {
        for (std::vector<Ion>::iterator ion = atom->ions.begin(); ion != atom->ions.end(); ++ion) {
            const bool continuum_ion_only = ((ion->ionization_stage == atom->max_ionization_stage+1) ? true : false);
            ion->read_atomic_data(continuum_ion_only);
        }
    }
    log_file << "done." << std::endl << std::endl;
}
