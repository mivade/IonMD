project(ionmd)

pybind11_add_module(${PROJECT_NAME} main.cpp)
target_link_libraries(
    ${PROJECT_NAME} PRIVATE
    ${CMAKE_THREAD_LIBS_INIT}
    ${ARMADILLO_LIBRARIES}
    libionmd
)
install(TARGETS ${PROJECT_NAME} DESTINATION ${CMAKE_INSTALL_PREFIX})
