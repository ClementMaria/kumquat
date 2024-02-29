project(algebra_example)

add_executable ( laurent_polynomial_example laurent_polynomial_example.cpp )
target_link_libraries( laurent_polynomial_example LINK_PUBLIC ${Boost_LIBRARIES} ${GMPXX_LIBRARIES} ${GMP_LIBRARIES})

# add_executable ( simple_rational_function_flint simple_rational_function_flint.cpp )
# target_link_libraries( simple_rational_function_flint LINK_PUBLIC ${Boost_LIBRARIES} ${GMPXX_LIBRARIES} ${GMP_LIBRARIES} ${FLINT_LIBRARIES})

# add_executable ( simple_complex_number simple_complex_number.cpp )
# target_link_libraries( simple_complex_number LINK_PUBLIC ${Boost_LIBRARIES} ${GMPXX_LIBRARIES} ${GMP_LIBRARIES})
