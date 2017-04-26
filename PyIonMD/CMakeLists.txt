project(pyionmd)

include_directories(${CMAKE_SOURCE_DIR})

configure_file(demo.py demo.py COPYONLY)

pybind11_add_module(${PROJECT_NAME} pyionmd.cpp)
target_link_libraries(${PROJECT_NAME} PRIVATE
  ionmd
  ${ARMADILLO_LIBRARIES}
)
