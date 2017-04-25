project(pyionmd)

include_directories(${CMAKE_SOURCE_DIR})

pybind11_add_module(${PROJECT_NAME} ionmd.cpp)
target_link_libraries(${PROJECT_NAME} PRIVATE
  ionmd
  ${ARMADILLO_LIBRARIES}
)
