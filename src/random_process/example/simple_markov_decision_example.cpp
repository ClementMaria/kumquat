/*    This file is part of the KUMQUAT Library - https://kumquat.inria.fr/ 
 *    - released under  GNU GENERAL PUBLIC LICENSE v.3
 *    See file LICENSE or go to https://kumquat.inria.fr/licensing/ for full 
 *    license details.
 *    Author(s):       Clément Maria
 *
 *    Copyright (C) 2023 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include <kumquat/Markov_decision.h>

using namespace kumquat;

int main(int argc, char * argv[]) {

  std::vector<size_t> num_states = {10,3,10,5,10,8,10};
  Markov_decision mdp(num_states);

  size_t sum=0;
  for(auto n : num_states) { sum += (n-1); }
  //amplify the sum of indices in the path, min is 0, max is 9+2+9+... = sum
  //return sum_state/sum - 1/2 in  [-0.5;0.5]
  Markov_decision::Proba_float rew_alpha = [&](std::vector<size_t> samp) {
    size_t sum_state=0; for(auto n : samp) { sum_state += (n-1); }
    return ((Proba_float)(sum_state)/(Proba_float)(sum) - (Proba_float)(0.5));
  };

  size_t num_iteration = 1000;
  mdp.train(num_iteration, rew_alpha);

  return 0;
}