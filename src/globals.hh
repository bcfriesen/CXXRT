#ifndef GLOBALS_HH
#define GLOBALS_HH

#include <vector>

#include <yaml-cpp/yaml.h>

extern std::vector<class GridVoxel> grid;
extern std::vector<class Ray> rays;

extern YAML::Node config;

extern std::vector<double> wavelength_values;

#endif
