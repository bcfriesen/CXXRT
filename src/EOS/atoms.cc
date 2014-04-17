#include <iostream>
#include <iomanip>
#include <string>

#include "atoms.hh"

Ion::Ion(const std::string atomic_symbol_in, const unsigned int atomic_number_in, const unsigned int ionization_stage_in, const double atomic_weight_in) {
  atomic_symbol = atomic_symbol_in;
  atomic_number = atomic_number_in;
  ionization_stage = ionization_stage_in;
  atomic_weight = atomic_weight_in;
}

std::ostream& operator<<(std::ostream& os, const Ion& ion) {
  os << std::setw(15) << "ION DATA: ";
  os << std::setw(20) << "atomic symbol" << std::setw(20) << "atomic number" << std::setw(20) << "ionization stage" << std::setw(20) << "atomic weight" << std::endl;
  os << std::setw(35) << ion.atomic_symbol << std::setw(20) << ion.atomic_number << std::setw(20) << ion.ionization_stage << std::setw(20) << ion.atomic_weight << std::endl;
  return os;
}
