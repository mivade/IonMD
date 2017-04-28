cmake_minimum_required(VERSION 3.6)

# set(CMAKE_INSTALL_PREFIX ${CMAKE_SOURCE_DIR}/dist CACHE PATH "Install prefix")
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")  # TODO: fix for MSVC

# project(IonMD)

find_package(Armadillo REQUIRED)
find_package(OpenMP REQUIRED)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

include_directories(
  ${ARMADILLO_INCLUDE_DIRS}
)

add_subdirectory(ionmd)
add_subdirectory(demo)
#add_subdirectory(pybind11)
#add_subdirectory(PyIonMD)
