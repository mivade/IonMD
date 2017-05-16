add_executable(tests test_data.cpp)
target_link_libraries(tests
    ${ARMADILLO_LIBRARIES}
    libionmd
    ${Boost_LIBRARIES}
)
