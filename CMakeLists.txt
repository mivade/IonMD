cmake_minimum_required(VERSION 3.6)

option(BUILD_PY "Build Python bindings" OFF)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")  # TODO: fix for MSVC

# find_package(Threads)  # necessary when not using OpenMP
find_package(Armadillo REQUIRED)
find_package(OpenMP)

if(${OPENMP_FOUND})
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif(${OPENMP_FOUND})

include_directories(
  include
  ${ARMADILLO_INCLUDE_DIRS}
)

add_subdirectory(libionmd)

if(BUILD_PY)
  add_subdirectory(pybind11)
  add_subdirectory(PyIonMD)
else()
  add_subdirectory(demo)
endif(BUILD_PY)
