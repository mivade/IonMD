project(ionmd)

add_library(${PROJECT_NAME} STATIC ion.cpp simulation.cpp)
target_link_libraries(${PROJECT_NAME} ${Boost_LIBRARIES})
