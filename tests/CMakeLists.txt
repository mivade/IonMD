add_executable(tests test_data.cpp)
target_link_libraries(tests
    ${HDF5_LIBRARIES}
    ${ARMADILLO_LIBRARIES}
    libionmd
    ${Boost_LIBRARIES}
)
