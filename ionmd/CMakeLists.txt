project(ionmd)

find_package(OpenMP REQUIRED)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

add_library(${PROJECT_NAME} STATIC ion.cpp simulation.cpp)
target_link_libraries(${PROJECT_NAME} ${Boost_LIBRARIES})
