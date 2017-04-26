project(ionmd)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")  # needed for pybind11

add_library(${PROJECT_NAME} ion.cpp simulation.cpp)
target_link_libraries(${PROJECT_NAME})
