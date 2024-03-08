project(random_process_example)

add_executable ( simple_markov_decision_example simple_markov_decision_example.cpp )
target_link_libraries( simple_markov_decision_example LINK_PUBLIC ${Boost_LIBRARIES} ${GMPXX_LIBRARIES} ${GMP_LIBRARIES} ${TBB_LIBRARIES} TBB::tbb)

#medea
# add_executable ( knot_markov_decision knot_markov_decision.cpp )
# target_link_libraries( knot_markov_decision LINK_PUBLIC ${Boost_LIBRARIES} ${GMPXX_LIBRARIES} ${GMP_LIBRARIES} ${TBB_LIBRARIES} TBB::tbb -L/usr/local/lib -lregina-engine  /Library/Developer/CommandLineTools/SDKs/MacOSX14.2.sdk/usr/lib/libxml2.tbd /usr/local/lib/libgmp.dylib /usr/local/lib/libgmpxx.dylib)


#MAC-
add_executable ( knot_markov_decision knot_markov_decision.cpp )
target_link_libraries( knot_markov_decision LINK_PUBLIC ${Boost_LIBRARIES} ${GMPXX_LIBRARIES} ${GMP_LIBRARIES} ${TBB_LIBRARIES} TBB::tbb -L/usr/local/lib -lregina-engine)# /Library/Developer/CommandLineTools/SDKs/MacOSX14.2.sdk/usr/lib/libxml2.tbd /opt/homebrew/lib/libgmp.dylib /opt/homebrew/lib/libgmpxx.dylib)

# add_executable ( test_markov_knot test_markov_knot.cpp )
# target_link_libraries( test_markov_knot LINK_PUBLIC ${Boost_LIBRARIES} ${GMPXX_LIBRARIES} ${GMP_LIBRARIES} ${TBB_LIBRARIES} TBB::tbb -L/usr/local/lib -lregina-engine /Library/Developer/CommandLineTools/SDKs/MacOSX14.2.sdk/usr/lib/libxml2.tbd /opt/homebrew/lib/libgmp.dylib /opt/homebrew/lib/libgmpxx.dylib)

add_executable ( amplify_example amplify_example.cpp )
target_link_libraries( amplify_example LINK_PUBLIC ${Boost_LIBRARIES} ${GMPXX_LIBRARIES} ${GMP_LIBRARIES} ${TBB_LIBRARIES} TBB::tbb)
