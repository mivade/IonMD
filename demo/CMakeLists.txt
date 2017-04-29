project(ionmd_demo)

add_executable(${PROJECT_NAME} demo.cpp)
target_link_libraries(${PROJECT_NAME}
  libionmd
  ${ARMADILLO_LIBRARIES}
)
