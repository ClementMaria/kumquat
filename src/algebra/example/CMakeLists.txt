project(algebra_example)

add_executable ( finite_ring_Z6_example finite_ring_Z6_example.cpp )

add_executable ( infinite_ring_Z_example infinite_ring_Z_example.cpp )
target_link_libraries( infinite_ring_Z_example LINK_PUBLIC ${Boost_LIBRARIES} ${GMPXX_LIBRARIES} ${GMP_LIBRARIES})

add_executable ( ring_Q_mod_Z_example ring_Q_mod_Z_example.cpp )
target_link_libraries( ring_Q_mod_Z_example LINK_PUBLIC ${Boost_LIBRARIES} ${GMPXX_LIBRARIES} ${GMP_LIBRARIES})