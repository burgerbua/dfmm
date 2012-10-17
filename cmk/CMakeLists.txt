# build type: (RELEASE || DEBUG)
#set (CMAKE_BUILD_TYPE RELEASE)
#set (CMAKE_BUILD_TYPE DEBUG)


set (OWN_R_FLAGS "-ffast-math") # evtl "-pg" for gprof
set (OWN_D_FLAGS "-Wall -Wno-sign-compare -fno-inline -Winvalid-pch -DDEBUG")

######################################################################
# if you do not really know what you do do not edit beyond this line #
######################################################################

#message ("\n\t*** Generating ${CMAKE_BUILD_TYPE} version of Makefiles. ***\n")

#---------------------------------------------------------------------
set (CMAKE_CXX_COMPILER g++)
project (dFMM)
enable_language (CXX Fortran)
message(STATUS "C++ compiler: ${CMAKE_CXX_COMPILER_ID}" )
#---------------------------------------------------------------------


# specify
cmake_minimum_required (VERSION 2.8)
if (COMMAND cmake_policy)
  cmake_policy (SET CMP0003 NEW)
endif (COMMAND cmake_policy)


# OS
set (CMAKE_SYSTEM LINUX)


# gprof
option (DFMM_WITH_GPROF "Set to ON to set compiler flags for GProf" OFF)
message (STATUS "Compiling for GProf: ${DFMM_BUILD_DEBUG}")
if (DFMM_WITH_GPROF)
  set (OWN_R_FLAGS "${OWN_R_FLAGS} -pg")
  set (OWN_D_FLAGS "${OWN_D_FLAGS} -pg")
endif()


# debugging symbols
option (DFMM_DEBUGGING_SYMBOLS "Set to ON to generate debugging symbols" OFF)
message (STATUS "Generating debugging symbols: ${DFMM_DEBUGGING_SYMBOLS}")
if (DFMM_DEBUGGING_SYMBOLS)
  set (OWN_R_FLAGS "${OWN_R_FLAGS} -g")
  set (OWN_D_FLAGS "${OWN_D_FLAGS} -g")
endif()


# debug or release
option (DFMM_BUILD_DEBUG "Set to ON to build in debug mode" OFF)
message (STATUS "Debug build: ${DFMM_BUILD_DEBUG}")
if (DFMM_BUILD_DEBUG)
  set (CMAKE_BUILD_TYPE DEBUG)
  set (CMAKE_CXX_FLAGS_DEBUG   "-O0 ${OpenMP_CXX_FLAGS} ${OWN_D_FLAGS}")
else()
  set (CMAKE_BUILD_TYPE RELEASE)
  set (CMAKE_C_FLAGS_RELEASE   "-O2 ${OpenMP_C_FLAGS}")
  set (CMAKE_CXX_FLAGS_RELEASE "-O2 ${OpenMP_CXX_FLAGS} ${OWN_R_FLAGS}")
endif()


# testing
option (DFMM_TESTING "Set to ON to enable testing" OFF)
if (DFMM_TESTING)
  enable_testing()
endif()


# count floating point operations
option (DFMM_COUNT_FLOPS "Set to ON to enable counting of flops" OFF)
if (DFMM_COUNT_FLOPS)
  add_definitions(-DDFMM_COUNT_FLOPS)
endif()

# Boost
include (FindBoost)

# OpenMP
include (FindOpenMP)
option (DFMM_USE_OMP "Set to ON to enable OpenMP parallelization" OFF)
if (DFMM_USE_OMP)
  add_definitions(-DDFMM_USE_OMP)
endif()


## accuracy
#set (DFMM_ACCURACY "5" CACHE STRING "Set the accuracy ACC (l,eps)=(ACC,1e-ACC)")
#if (DFMM_ACCURACY EQUAL 3)
# add_definitions (-DDFMM_EPSILON=1e-3 -DDFMM_ORDER=3)
#elseif (DFMM_ACCURACY EQUAL 4)
# add_definitions (-DDFMM_EPSILON=1e-4 -DDFMM_ORDER=4)
#elseif (DFMM_ACCURACY EQUAL 5)
# add_definitions (-DDFMM_EPSILON=1e-5 -DDFMM_ORDER=5)
#elseif (DFMM_ACCURACY EQUAL 6)
# add_definitions (-DDFMM_EPSILON=1e-6 -DDFMM_ORDER=6)
#elseif (DFMM_ACCURACY EQUAL 7)
# add_definitions (-DDFMM_EPSILON=1e-7 -DDFMM_ORDER=7)
#elseif (DFMM_ACCURACY EQUAL 8)
# add_definitions (-DDFMM_EPSILON=1e-8 -DDFMM_ORDER=8)
#elseif (DFMM_ACCURACY EQUAL 9)
# add_definitions (-DDFMM_EPSILON=1e-9  -DDFMM_ORDER=9)
#elseif (DFMM_ACCURACY EQUAL 10)
# add_definitions (-DDFMM_EPSILON=1e-10 -DDFMM_ORDER=10)
#elseif (DFMM_ACCURACY EQUAL 11)
# add_definitions (-DDFMM_EPSILON=1e-11 -DDFMM_ORDER=11)
#elseif (DFMM_ACCURACY EQUAL 12)
# add_definitions (-DDFMM_EPSILON=1e-12 -DDFMM_ORDER=12)
#elseif (DFMM_ACCURACY EQUAL 13)
# add_definitions (-DDFMM_EPSILON=1e-13 -DDFMM_ORDER=13)
#else ()
# message (FATAL_ERROR "Set accuracy ACC is not valid!")
#endif ()



# interpolation order
set (DFMM_INTERPOLATIONORDER "5" CACHE STRING "Set the interpolation order")
if (DFMM_INTERPOLATIONORDER EQUAL 3)
  add_definitions (-DDFMM_ORDER=3)
elseif (DFMM_INTERPOLATIONORDER EQUAL 4)
  add_definitions (-DDFMM_ORDER=4)
elseif (DFMM_INTERPOLATIONORDER EQUAL 5)
  add_definitions (-DDFMM_ORDER=5)
elseif (DFMM_INTERPOLATIONORDER EQUAL 6)
  add_definitions (-DDFMM_ORDER=6)
elseif (DFMM_INTERPOLATIONORDER EQUAL 7)
  add_definitions (-DDFMM_ORDER=7)
elseif (DFMM_INTERPOLATIONORDER EQUAL 8)
  add_definitions (-DDFMM_ORDER=8)
elseif (DFMM_INTERPOLATIONORDER EQUAL 9)
  add_definitions (-DDFMM_ORDER=9)
elseif (DFMM_INTERPOLATIONORDER EQUAL 10)
  add_definitions (-DDFMM_ORDER=10)
elseif (DFMM_INTERPOLATIONORDER EQUAL 11)
  add_definitions (-DDFMM_ORDER=11)
elseif (DFMM_INTERPOLATIONORDER EQUAL 12)
  add_definitions (-DDFMM_ORDER=12)
elseif (DFMM_INTERPOLATIONORDER EQUAL 13)
  add_definitions (-DDFMM_ORDER=13)
else ()
  message (FATAL_ERROR "Set interpolation order is not valid!")
endif ()

# epsilon
set (DFMM_LOWRANKEPSILON "5" CACHE STRING "Set the approximation accuracy")
if (DFMM_LOWRANKEPSILON EQUAL 1)
  add_definitions (-DDFMM_EPSILON=1e-1)
elseif (DFMM_LOWRANKEPSILON EQUAL 2)
  add_definitions (-DDFMM_EPSILON=1e-2)
elseif (DFMM_LOWRANKEPSILON EQUAL 3)
  add_definitions (-DDFMM_EPSILON=1e-3)
elseif (DFMM_LOWRANKEPSILON EQUAL 4)
  add_definitions (-DDFMM_EPSILON=1e-4)
elseif (DFMM_LOWRANKEPSILON EQUAL 5)
  add_definitions (-DDFMM_EPSILON=1e-5)
elseif (DFMM_LOWRANKEPSILON EQUAL 6)
  add_definitions (-DDFMM_EPSILON=1e-6)
elseif (DFMM_LOWRANKEPSILON EQUAL 7)
  add_definitions (-DDFMM_EPSILON=1e-7)
elseif (DFMM_LOWRANKEPSILON EQUAL 8)
  add_definitions (-DDFMM_EPSILON=1e-8)
elseif (DFMM_LOWRANKEPSILON EQUAL 9)
  add_definitions (-DDFMM_EPSILON=1e-9)
elseif (DFMM_LOWRANKEPSILON EQUAL 10)
  add_definitions (-DDFMM_EPSILON=1e-10)
elseif (DFMM_LOWRANKEPSILON EQUAL 11)
  add_definitions (-DDFMM_EPSILON=1e-11)
elseif (DFMM_LOWRANKEPSILON EQUAL 12)
  add_definitions (-DDFMM_EPSILON=1e-12)
elseif (DFMM_LOWRANKEPSILON EQUAL 13)
  add_definitions (-DDFMM_EPSILON=1e-13)
else ()
  message (FATAL_ERROR "Set epsilon is not valid!")
endif ()




# use mkl
option (DFMM_USE_MKL "Set to ON to use MKL" OFF )
if (DFMM_USE_MKL)
  add_definitions (-DDFMM_USE_MKL)
  set (BLAS_LIBRARIES "-L$ENV{MKLROOT}/lib -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm")
else () 
  include (FindBLAS)
  include (FindLAPACK)
endif()

# set directories
set (BASE_PATH        ${CMAKE_SOURCE_DIR}/..)
set (EXECUTABLE_PATH  ${CMAKE_SOURCE_DIR})
file (TO_CMAKE_PATH ~/bin BIN_PATH)


# add subdirectories
add_subdirectory (${BASE_PATH}/examples ./examples)
add_subdirectory (${BASE_PATH}/doc      ./doc)


