set (misc_SOURCES
    atomic_symbols.cc
    calc_H_target.cc
    calc_emergent_spectrum.cc
    calc_Delta_T.cc
    planck_function.cc
    planck_function_temperature_derivative.cc
    )

add_library (misc ${misc_SOURCES})
