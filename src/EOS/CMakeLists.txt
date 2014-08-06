set (EOS_SOURCES
     LTE_EOS.cc
     atoms.cc
     photoionization_cross_section.cc
     saha_equation.cc
    )

add_library (EOS ${EOS_SOURCES})
