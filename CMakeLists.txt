cmake_minimum_required(VERSION 3.6)

project(IonMD)

find_package(Armadillo REQUIRED)
find_package(Boost REQUIRED)

set(SRCDIR ${CMAKE_SOURCE_DIR}/src)

set(SRCS
  src/ion.cpp
  src/ionmd.cpp
)

include_directories(
  ${ARMADILLO_INCLUDE_DIRS}
  ${BOOST_INCLUDE_DIRS}
)

add_library(${PROJECT_NAME} SHARED ${SRCS})

target_link_libraries(${PROJECT_NAME} ${ARMADILLO_LIBRARIES})