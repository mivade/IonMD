cmake_minimum_required(VERSION 3.6)

set(PYBIND ON)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")  # TODO: fix for MSVC
set(Boost_USE_STATIC_LIBS ON)

find_package(Armadillo REQUIRED)
find_package(Boost COMPONENTS log REQUIRED)
find_package(OpenMP REQUIRED)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

include_directories(
  ${ARMADILLO_INCLUDE_DIRS}
  ${BOOST_INCLUDE_DIRS}
)

add_subdirectory(ionmd)
add_subdirectory(demo)

if(PYBIND)
  add_subdirectory(pybind11)
  add_subdirectory(PyIonMD)
endif(PYBIND)