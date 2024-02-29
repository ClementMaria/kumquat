/*    This file is part of the KUMQUAT Library -
 *    https://kumquat.inria.fr/ 
 *    - which is a licence protected library. See file LICENSE 
 *    or go to https://kumquat.inria.fr/licensing/ for full 
 *    license details.
 *
 *    Author(s):       Clément Maria
 *
 *    Copyright (C) 2024 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include <kumquat/Markov_decision_quantum_knot.h>
#include <iostream>
#include <fstream>

using namespace kumquat;

typedef Markov_decision_quantum_knot::Proba_float P_float;

int main(int argc, char * argv[]) {

  int max_twists = 2;
  int num_strands = 4;//must be even
  int num_layers = 4;//better if odd
  int num_iterations = 10000;
/****************************************************************/
//
  Markov_decision_quantum_knot mdp(num_strands, num_layers, max_twists);
//


  std::vector<size_t> state = {1,19,4,23};
  auto b_simp = mdp.state_to_braid(state);
  auto b_verb = mdp.state_to_braid(state,false);

  std::cout << "Braid simplified:\n";
  std::cout << b_simp << "\n";
  std::cout << "Braid verbatim:\n";
  std::cout << b_verb << "\n\n";
  
  Plat_braid b = b_verb;

  auto jones = mdp.jones_poly_.quantum_invariant(b);

  int layer_idx = 2;
  auto Mb = mdp.jones_poly_.quantum_morphism(b.braid().begin(), b.braid().begin()+layer_idx);//Mb

  auto Mt = mdp.jones_poly_.quantum_morphism(b.braid().begin()+(layer_idx), b.braid().end());  

  Mt.ltimes_equal(mdp.jones_poly_.h_tensor());//h * Mt

  // Mt.ltimes_equal(Mb);

  std::cout << Mt << "\n\n";


  std::cout << "Traditional Jones =    " << jones << "\n";
  std::cout << "Other Jones       =    " << Mt.trace() << "\n";

  // Mt.



/****************************************************************/
  //train until possibly stuck in a state
  // size_t num_iteration = (size_t)num_iterations;
//
  // auto max_samp = mdp.train(num_iteration);

  // auto imp_samp = mdp.improve_exhaustive(max_samp);

  return 0;
}


