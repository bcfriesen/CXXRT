set (EOS_SOURCES
     LTE_EOS.cc
     atoms.cc
     insert_ions.cc
     saha_equation.cc
     Phi.cc
    )

add_library (EOS ${EOS_SOURCES})
