project(matrix_example)

add_executable ( simple_dense_matrix_example simple_dense_matrix_example.cpp )

add_executable ( dense_matrix_pid_example dense_matrix_pid_example.cpp )
target_link_libraries( dense_matrix_pid_example LINK_PUBLIC ${Boost_LIBRARIES} ${GMPXX_LIBRARIES} ${GMP_LIBRARIES})