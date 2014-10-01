#ifndef GLOBALS_HH
#define GLOBALS_HH

#include <vector>
#include <map>

#include "EOS/atoms.hh"

#include <yaml-cpp/yaml.h>

extern std::vector<class GridVoxel> grid;
extern std::vector<class Ray> rays;

extern YAML::Node config;

extern std::vector<double> tot_wavelength_values;
extern std::vector<double> wavelength_values;

extern std::ofstream log_file;

extern std::ofstream moments_file;

extern std::map<std::string, unsigned int> atomic_symbols;

extern std::vector<Atom> tmp_atoms;

#endif
