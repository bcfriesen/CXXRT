#include "../globals.hh"
#include "../grid.hh"

void insert_ions(GridVoxel* gv) {
    std::vector<unsigned int> ions;
    for (YAML::const_iterator it = config["atoms_with_max_ion"].begin(); it != config["atoms_with_max_ion"].end(); ++it) {
        // Make sure we don't try to add an atom which is more than fully ionized.
        if (it->second.as<unsigned int>() >= atomic_symbols[it->first.as<std::string>()]) {
            std::cerr << "ERROR: you requested atom " << it->first.as<std::string>() << " with ionization stage " << it->second.as<unsigned int>() << std::endl;
            exit(1);
        }
        Atom one_atom(atomic_symbols[it->first.as<std::string>()], it->second.as<unsigned int>());
        gv->atoms.push_back(one_atom);
    }
}
