project(ionmd_demo)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

include_directories(${CMAKE_SOURCE_DIR})

add_executable(${PROJECT_NAME} demo.cpp)
target_link_libraries(${PROJECT_NAME}
  ionmd
  ${ARMADILLO_LIBRARIES}
)
