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

#include <kumquat/Markov_decision.h>

using namespace kumquat;

typedef Markov_decision::Proba_float P_float;

int main() {
  //num states per layer
  std::vector<size_t> num_states = {20,39,50,150,400,80,10};
  Markov_decision mdp(num_states);

  // size_t sum=0;
  // for(auto n : num_states) { sum += (n-1); }
  // //amplify the sum of indices in the path, min is 0, max is sum
  // //return 0.01*[sum_state/sum - 1/2] in  [-0.005;0.005]
  // auto rew_sum = [&](std::vector<size_t> &samp)->P_float
  // {
  //   size_t sum_state=0; 
  //   for(auto n : samp) { sum_state += n; }
  //   // return 0;
  //   return 0.01*((P_float)(sum_state)/(P_float)(sum) - (P_float)(0.5));
  // };

  //non-linear reward based on a polynomial on the labels of nodes in the Markov process
  auto rew_poly = [&](std::vector<size_t> &samp)->P_float
  {
    P_float s0 = (P_float)samp[0];
    P_float s1 = (P_float)samp[1];
    P_float s2 = (P_float)samp[2];
    P_float s3 = (P_float)samp[3];
    P_float s4 = (P_float)samp[4];
    P_float s5 = (P_float)samp[5];
    P_float s6 = (P_float)samp[6];

    auto p_x = s0*s1*s5 - s0*s3*s4 + s0*s5*s6 - s2*s5*s1 + s3*s4*s2;
     
    return (p_x / (P_float)(150*400*40))-0.1;
  };

  //evaluate randomly the polynomial reward to later normalize it
  std::cout << "Out of 1000 rewards: ";
  P_float minr = 0;
  P_float maxr = 0; 
  for(size_t k=0; k<1000; ++k) {
    auto samp = mdp.sample();
    auto tmpr = rew_poly(samp);
    if(tmpr > maxr) { maxr = tmpr; }
    if(tmpr < minr) { minr = tmpr; }
    // std::cout << tmpr << " ";
  }
  std::cout << "\n" << "Min reward = " << minr << "  Max reward = " << maxr << "\n";

  //normalized polynomial reward
  auto rew_poly_balanced = [&](std::vector<size_t> &samp)->P_float
  {
    auto rew = rew_poly(samp);
    rew -= (maxr+minr)/2.;
    rew *= 0.01 * (1./(maxr-minr));
    return rew;
  };

  //train until possibly stuck in a state
  size_t num_iteration = 30000;//generally converges
  size_t batch_size = 10;
  auto max_samp = mdp.train(num_iteration, rew_poly_balanced, batch_size);

  //sample the state and brute force improve it a few times
  std::vector<size_t> curr_samp;
  auto imp_samp = max_samp;
  do {
    curr_samp = imp_samp;
    imp_samp = mdp.local_improve(curr_samp, rew_poly_balanced);
  } while(curr_samp != imp_samp);

  // auto samp = mdp.sample();
  // auto imp_samp = mdp.local_improve(max_samp, rew_poly);

  std::cout << "\n\n\n";
  std::cout << "Sample          = " << mdp.to_string(max_samp) << "of rew == " << rew_poly(max_samp) <<  "\n";
  std::cout << "Improved sample = " << mdp.to_string(imp_samp) << "of rew == " << rew_poly(imp_samp) <<  "\n";
  std::cout << "\n";

  // //start a new Markow decision process
  // Markov_decision mdp_1(num_states);
  // //and push the imp_sample first to give it a direction
  // mdp_1.amplify(imp_samp, 2.*rew_poly_balanced(imp_samp));
  // //train to try to do better
  // mdp_1.train(num_iteration, rew_poly_balanced);

  return 0;
}