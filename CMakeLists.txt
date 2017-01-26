cmake_minimum_required(VERSION 3.6)

set(CMAKE_CXX_STANDARD 11)

find_package(Armadillo REQUIRED)
find_package(Boost REQUIRED)
find_package(OpenMP REQUIRED)

include_directories(
  ${ARMADILLO_INCLUDE_DIRS}
  ${BOOST_INCLUDE_DIRS}
)

add_subdirectory(ionmd)
add_subdirectory(demo)
