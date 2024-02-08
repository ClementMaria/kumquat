project(quantum_topology_example)

add_executable ( simple_quantum_data_example simple_quantum_data_example.cpp )
target_link_libraries( simple_quantum_data_example LINK_PUBLIC ${Boost_LIBRARIES} ${GMPXX_LIBRARIES} ${GMP_LIBRARIES})

add_executable( simple_quantum_group_example simple_quantum_group_example.cpp )
target_link_libraries( simple_quantum_group_example LINK_PUBLIC ${Boost_LIBRARIES} ${GMPXX_LIBRARIES} ${GMP_LIBRARIES})


add_executable( simple_jones_polynomial_example simple_jones_polynomial_example.cpp )
target_link_libraries( simple_jones_polynomial_example LINK_PUBLIC ${Boost_LIBRARIES} ${GMPXX_LIBRARIES} ${GMP_LIBRARIES})
