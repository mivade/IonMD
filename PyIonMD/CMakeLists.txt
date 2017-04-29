project(ionmd)

pybind11_add_module(${PROJECT_NAME} main.cpp)
target_link_libraries(
    ${PROJECT_NAME} PRIVATE
    ${ARMADILLO_LIBRARIES} libionmd
)
