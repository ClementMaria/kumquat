project(algebra_example)

add_executable ( finite_ring_Z6_example finite_ring_Z6_example.cpp )

add_executable ( infinite_ring_Z_example infinite_ring_Z_example.cpp )
target_link_libraries( infinite_ring_Z_example LINK_PUBLIC ${Boost_LIBRARIES} ${GMPXX_LIBRARIES} ${GMP_LIBRARIES})

add_executable ( ring_Q_mod_Z_example ring_Q_mod_Z_example.cpp )
target_link_libraries( ring_Q_mod_Z_example LINK_PUBLIC ${Boost_LIBRARIES} ${GMPXX_LIBRARIES} ${GMP_LIBRARIES})

add_executable ( number_theory_example number_theory_example.cpp )
target_link_libraries( number_theory_example LINK_PUBLIC ${Boost_LIBRARIES} ${GMPXX_LIBRARIES} ${GMP_LIBRARIES})

add_executable ( laurent_polynomial_example laurent_polynomial_example.cpp )
target_link_libraries( laurent_polynomial_example LINK_PUBLIC ${Boost_LIBRARIES} ${GMPXX_LIBRARIES} ${GMP_LIBRARIES})

# add_executable ( simple_C_example simple_C_example.cpp )
# target_link_libraries( simple_C_example LINK_PUBLIC ${Boost_LIBRARIES} ${GMPXX_LIBRARIES} ${GMP_LIBRARIES})
