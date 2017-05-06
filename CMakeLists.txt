cmake_minimum_required(VERSION 3.6)

option(BUILD_PY "Build Python bindings" OFF)

set(CMAKE_CXX_STANDARD 11)
if(${CMAKE_COMPILER_IS_GNUCXX})
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
endif(${CMAKE_COMPILER_IS_GNUCXX})

find_package(Armadillo REQUIRED)

find_package(OpenMP)
if(${OPENMP_FOUND})
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
else()
  find_package(Threads)  # necessary when not using OpenMP
endif(${OPENMP_FOUND})

# find_package(HDF5 REQUIRED)

include_directories(
  include
  ${ARMADILLO_INCLUDE_DIRS}
  # ${HDF5_INCLUDE_DIRS}
)

add_subdirectory(libionmd)

if(BUILD_PY)
  add_subdirectory(pybind11)
  add_subdirectory(PyIonMD)
else()
  add_subdirectory(demo)
endif(BUILD_PY)
