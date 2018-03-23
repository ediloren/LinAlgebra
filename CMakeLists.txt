set(LinAlgebra_Sources
    Linalg.h
    LUDecomp.cpp
    LUDecomp.h
    Mtx.cpp
    Mtx.h
    Range.h
    Vect.cpp
    Vect.h)

#set(LinAlgebra_Required_Libs
#   debug MSVCRTD.LIB
#   debug MSVCPRTD.LIB
#   optimized MSVCRT.LIB
#   optimized MSVCPRT.LIB)


add_library(LinAlgebra ${LinAlgebra_Sources})
# libraries that needs to be linked to the present library
#target_link_libraries(LinAlgebra ${LinAlgebra_Required_Libs})



