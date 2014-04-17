#include <cstdio>
#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>

#include "atoms.hh"
#include "../globals.hh"

Ion::Ion(const std::string atomic_symbol_in, const unsigned int atomic_number_in, const unsigned int ionization_stage_in, const double atomic_weight_in) {
  atomic_symbol = atomic_symbol_in;
  atomic_number = atomic_number_in;
  ionization_stage = ionization_stage_in;
  atomic_weight = atomic_weight_in;

  const std::string atomic_data_file_name = config["atomic_data_root_folder"].as<std::string>() + "/" + atomic_symbol + "_" + std::to_string(ionization_stage+1) + ".dat";
  std::ifstream atomic_data_file;
  atomic_data_file.open(atomic_data_file_name);
  std::string one_line;
  double wavelength, log_gf, el_code, first_energy_level, J_first, second_energy_level, J_second, gam_rad, gam_stark, gam_vdw, hyperfine, iso_abund_frac;
  char *blank = nullptr, *label = nullptr, *ref = nullptr, *code = nullptr;
  int first_nlte_lvl_idx, second_nlte_lvl_idx, isotope_num, isotope_num2, first_hyper_shift, second_hyper_shift, first_hyperfine, second_hyperfine, icode, first_lande, second_lande, iso_shift;
  std::getline(atomic_data_file, one_line);
  while (std::getline(atomic_data_file, one_line)) {
    int ret = std::sscanf(one_line.c_str(), "%11lf%7lf%6lf%12lf%5lf%1s%10s%12lf%5lf%1s%10s%6lf%6lf%6lf%4s%2d%2d%3d%6lf%3d%6lf%5d%5d%1s%1d%1s%1s%1d%1s%1d%3s%5d%5d%6d",
                                             &wavelength,
                                             &log_gf,
                                             &el_code,
                                             &first_energy_level,
                                             &J_first,
                                             blank,
                                             label,
                                             &second_energy_level,
                                             &J_second,
                                             blank,
                                             label,
                                             &gam_rad,
                                             &gam_stark,
                                             &gam_vdw,
                                             ref,
                                             &first_nlte_lvl_idx,
                                             &second_nlte_lvl_idx,
                                             &isotope_num,
                                             &hyperfine,
                                             &isotope_num2,
                                             &iso_abund_frac,
                                             &first_hyper_shift,
                                             &second_hyper_shift,
                                             blank,
                                             &first_hyperfine,
                                             blank,
                                             blank,
                                             &second_hyperfine,
                                             blank,
                                             &icode,
                                             code,
                                             &first_lande,
                                             &second_lande,
                                             &iso_shift);

    std::cout << first_energy_level << " " << second_energy_level << std::endl;

    class AtomicLine atomic_line;
    class AtomicLevel upper_level, lower_level;
    lower_level.energy = first_energy_level;
    lower_level.J = J_first;
    lower_level.g = int(2.0*lower_level.J + 1);
    upper_level.energy = second_energy_level;
    upper_level.J = J_second;
    upper_level.g = int(2.0*upper_level.J + 1);

    atomic_line.wavelength = wavelength;
  }
  atomic_data_file.close();
}

std::ostream& operator<<(std::ostream& os, const Ion& ion) {
  os << std::setw(15) << "ION DATA: ";
  os << std::setw(20) << "atomic symbol" << std::setw(20) << "atomic number" << std::setw(20) << "ionization stage" << std::setw(20) << "atomic weight" << std::endl;
  os << std::setw(35) << ion.atomic_symbol << std::setw(20) << ion.atomic_number << std::setw(20) << ion.ionization_stage << std::setw(20) << ion.atomic_weight << std::endl;
  return os;
}
