project(random_process_example)

add_executable ( simple_markov_decision_example simple_markov_decision_example.cpp )
target_link_libraries( simple_markov_decision_example LINK_PUBLIC ${Boost_LIBRARIES} ${GMPXX_LIBRARIES} ${GMP_LIBRARIES} ${TBB_LIBRARIES} TBB::tbb)

add_executable ( knot_markov_decision knot_markov_decision.cpp )
target_link_libraries( knot_markov_decision LINK_PUBLIC ${Boost_LIBRARIES} ${GMPXX_LIBRARIES} ${GMP_LIBRARIES} ${TBB_LIBRARIES} TBB::tbb)

add_executable ( amplify_example amplify_example.cpp )
target_link_libraries( amplify_example LINK_PUBLIC ${Boost_LIBRARIES} ${GMPXX_LIBRARIES} ${GMP_LIBRARIES} ${TBB_LIBRARIES} TBB::tbb)
