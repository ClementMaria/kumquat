project(matrix_example)

add_executable ( simple_dense_matrix_example simple_dense_matrix_example.cpp )

add_executable ( dense_matrix_pid_example dense_matrix_pid_example.cpp )
target_link_libraries( dense_matrix_pid_example LINK_PUBLIC ${Boost_LIBRARIES} ${GMPXX_LIBRARIES} ${GMP_LIBRARIES})

add_executable ( dense_matrix_Qp_mod_Z_gram_normalization_example dense_matrix_Qp_mod_Z_gram_normalization_example.cpp )
target_link_libraries( dense_matrix_Qp_mod_Z_gram_normalization_example LINK_PUBLIC ${Boost_LIBRARIES} ${GMPXX_LIBRARIES} ${GMP_LIBRARIES})

# add_executable( simple_dense_matrix_eigen simple_dense_matrix_eigen.cpp )
# target_link_libraries( simple_dense_matrix_eigen ${Boost_LIBRARIES} ${GMPXX_LIBRARIES} ${GMP_LIBRARIES} ${EIGEN_LIBRARIES})

add_executable ( dense_matrix_multiplication_example dense_matrix_multiplication_example.cpp )
target_link_libraries( dense_matrix_multiplication_example LINK_PUBLIC ${Boost_LIBRARIES} ${GMPXX_LIBRARIES} ${GMP_LIBRARIES}) 
#${TBB_LIBRARIES})

add_executable ( test_example test_example.cpp )
target_link_libraries( test_example ${TBB_LIBRARIES} TBB::tbb)