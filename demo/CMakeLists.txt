project(ionmd_demo)

add_executable(${PROJECT_NAME} demo.cpp)
target_link_libraries(${PROJECT_NAME}
  ${CMAKE_THREAD_LIBS_INIT}
  ${ARMADILLO_LIBRARIES}
  ${Boost_LIBRARIES}
  libionmd
)
