project(random_process_example)

add_executable ( simple_markov_decision_example simple_markov_decision_example.cpp )
target_link_libraries( simple_markov_decision_example LINK_PUBLIC ${Boost_LIBRARIES} ${GMPXX_LIBRARIES} ${GMP_LIBRARIES})
