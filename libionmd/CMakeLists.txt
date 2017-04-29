project(libionmd)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")  # needed for pybind11

add_library(${PROJECT_NAME} ion.cpp simulation.cpp)
