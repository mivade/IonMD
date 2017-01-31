project(ionmd)

add_library(${PROJECT_NAME} STATIC ion.cpp ionmd.cpp)
target_link_libraries(${PROJECT_NAME} ${Boost_LIBRARIES})
