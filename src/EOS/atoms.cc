#include <cstdio>
#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <limits>

#include "atoms.hh"
#include "../globals.hh"
#include "../constants.hh"

Atom::Atom(const unsigned int atomic_number_in)
  : atomic_number(atomic_number_in)
{
  const std::string atomic_data_file_name = config["atomic_data_root_folder"].as<std::string>() + "/" + "atoms.dat";
  std::ifstream atomic_data_file;
  std::string one_line;
  atomic_data_file.open(atomic_data_file_name.c_str()); // in C++11 the argument of open() can be a string
  unsigned int atomic_number_from_file;
  while (std::getline(atomic_data_file, one_line)) {
    std::istringstream iss(one_line);
    iss >> atomic_number_from_file;
    if (atomic_number_from_file == atomic_number)
      iss >> atomic_weight >> atomic_symbol;
  }
  atomic_data_file.close();

  for (unsigned int i = 0; i <= atomic_number_in; ++i) {
    Ion ion(atomic_number, i);
    ions.push_back(ion);
  }
}

Ion::Ion(const unsigned int atomic_number_in, const unsigned int ionization_stage_in)
  : atomic_number(atomic_number_in),
    ionization_stage(ionization_stage_in)
{
  if (ionization_stage > atomic_number) {
    std::cerr << "ERROR! in atom " << atomic_number << " you tried to add an ion with ionization stage " << ionization_stage << std::endl;
    exit(1);
  } else if (ionization_stage < atomic_number) {
    std::ostringstream convert; // for compilers without std::to_string (a C++11 feature)
    convert << atomic_number*100 + ionization_stage;
    const std::string atomic_data_file_name = config["atomic_data_root_folder"].as<std::string>() + "/" + convert.str() + ".dat";
    std::ifstream atomic_data_file;
    atomic_data_file.open(atomic_data_file_name.c_str()); // in C++ the argument of open() can be a string
    std::string one_line;
    double wavelength, log_gf, first_energy_level, J_first, second_energy_level, J_second;
    bool duplicate;
    const double tolerance = 1.0e-30;

    // The first line is the ionization potential.
    std::getline(atomic_data_file, one_line);
    std::istringstream iss(one_line);
    iss >> ionization_potential;
    ionization_potential *= h_planck * c_light;

    while (std::getline(atomic_data_file, one_line)) {
      std::istringstream iss(one_line);
      iss >> wavelength >> log_gf >> first_energy_level >> J_first >> second_energy_level >> J_second;

      class AtomicLevel upper_level, lower_level;
      lower_level.energy = first_energy_level * h_planck * c_light;  // convert cm^-1 to erg
      lower_level.J = J_first;
      lower_level.g = int(2.0*lower_level.J + 1.0);
      upper_level.energy = second_energy_level * h_planck * c_light; // convert cm^-1 to erg
      upper_level.J = J_second;
      upper_level.g = int(2.0*upper_level.J + 1.0);

      // Add the lower level of the line to the model atom only if it hasn't been
      // added already. We detect duplicates by checking energy differences; if
      // two levels are similar to within some tolerance, they are defined as
      // duplicates.

      for (auto level: levels) {
        double difference;
        if (level.energy < std::numeric_limits<double>::epsilon()) {
          difference = std::abs(lower_level.energy - level.energy);
        } else {
          difference = std::abs(lower_level.energy - level.energy) / level.energy;
        }
        if (difference > tolerance) {
          duplicate = false;
        } else {
          duplicate = true;
          break;
        }
      }
      if (!duplicate) {
        levels.push_back(lower_level);
      }
      // Now add the upper level if necessary.
      for (auto level: levels) {
        const double difference = std::abs(upper_level.energy - level.energy) / level.energy;
        if (difference > tolerance) {
          duplicate = false;
        } else {
          duplicate = true;
          break;
        }
      }
      if (!duplicate) {
        levels.push_back(upper_level);
      }
    }

    atomic_data_file.close();

    // Now read the file again to get the lines.

    atomic_data_file.open(atomic_data_file_name.c_str()); // in C++ the argument of open() can be a string

    // The first line is the ionization potential so skip it this time.
    std::getline(atomic_data_file, one_line);

    while (std::getline(atomic_data_file, one_line)) {
      std::istringstream iss(one_line);
      iss >> wavelength >> log_gf >> first_energy_level >> J_first >> second_energy_level >> J_second;

      class AtomicLine atomic_line;
      atomic_line.wavelength = wavelength;
      const int g = int(2.0*J_first + 1.0);
      atomic_line.oscillator_strength = std::pow(10.0, log_gf) / double(g);

      // Now search through the model atom energy levels to find the lower and upper levels for this line.
      for (auto level = levels.begin(); level != levels.end(); ++level) {
        if (level->energy < std::numeric_limits<double>::epsilon()) {  // if we're comparing to the ground state (which has energy 0), don't divide by it
          if (std::abs(first_energy_level*h_planck*c_light - level->energy) < tolerance)
            atomic_line.lower_level = &(*level); // std::addressof(*level) would work but that's C++11
        } else if ((std::abs(first_energy_level*h_planck*c_light - level->energy) / level->energy) < tolerance) {
          atomic_line.lower_level = &(*level); // std::addressof(*level) would work but that's C++11
        }
      }
      for (auto level = levels.begin(); level != levels.end(); ++level) {
        if (level->energy < std::numeric_limits<double>::epsilon()) {  // if we're comparing to the ground state (which has energy 0), don't divide by it
          if (std::abs(second_energy_level*h_planck*c_light - level->energy) < tolerance)
            atomic_line.upper_level = &(*level); // std::addressof(*level) would work but that's C++11
        } else if ((std::abs(second_energy_level*h_planck*c_light - level->energy) / level->energy) < tolerance) {
          atomic_line.upper_level = &(*level); // std::addressof(*level) would work but that's C++11
        }
      }

      duplicate = false;
      // We assume the line lists don't have duplicate lines, but it can't hurt to check anyway.
      for (auto line: lines) {
        if ((std::abs(atomic_line.wavelength - line.wavelength) / line.wavelength) > tolerance) {
          duplicate = false;
        } else {
          duplicate = true;
          break;
        }
      }
      if (!duplicate) {
        lines.push_back(atomic_line);
      }
    }
    atomic_data_file.close();
  }
}


void Ion::calc_partition_function(const double temperature) {
    const double beta = 1.0 / (k_boltzmann * temperature);
    double result = 0.0;
    for (auto level: levels) {
        result += level.g * std::exp(-beta * level.energy);
    }
    partition_function = result;
}

std::ostream& operator<<(std::ostream& os, const Ion& ion) {
  os << std::setw(15) << "ION DATA: ";
  os << std::setw(20) << "atomic symbol" << std::setw(20) << "atomic number" << std::setw(20) << "ionization stage" << std::setw(20) << "atomic weight" << std::endl;
  os << std::setw(35) << ion.atomic_symbol << std::setw(20) << ion.atomic_number << std::setw(20) << ion.ionization_stage << std::setw(20) << ion.atomic_weight << std::endl;
  return os;
}
