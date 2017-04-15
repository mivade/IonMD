project(ionmd_demo)

include_directories(${CMAKE_SOURCE_DIR})

add_executable(${PROJECT_NAME} demo.cpp)
target_link_libraries(${PROJECT_NAME}
  ionmd
  ${ARMADILLO_LIBRARIES}
)
