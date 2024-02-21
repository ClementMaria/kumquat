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
#include <kumquat/Braid.h>
#include <kumquat/Jones_polynomial.h>
#include <iostream>
#include <fstream>
#include <ctime>

using namespace kumquat;

typedef Markov_decision::Proba_float P_float;

int main(int argc, char * argv[]) {

  if(argc != 5) { 
    std::cout << "./exe num_strands num_layers max_twists num_iterations\n";
    return 0;
  }

  int max_twists = std::atoi(argv[3]);
//10;//arbitrary
  int num_strands = std::atoi(argv[1]);//must be even
  int num_layers = std::atoi(argv[2]);//must be odd
  int num_iterations = std::atoi(argv[4]);

  int max_num_crossings = 0;
  if(num_layers % 2 == 0) {
    max_num_crossings = max_twists * ( (num_layers/2) * (num_strands/2) + (num_layers/2 + 1)*(num_layers/2-1) );
  }
  else {
    max_num_crossings = max_twists * ( (num_layers/2) * (num_strands/2) + (num_layers/2)*(num_layers/2-1) );
  }

  int sum_n_squared = max_num_crossings * (max_num_crossings-1)/2;

  std::cout << "---------- Nb. strands = " << num_strands << "  Nb. layers = " << num_layers << "  Max nb. twists = " << max_twists << " ----------\n";

  std::cout << "---> max # crossings = " << max_num_crossings << "\n";
  //plat braids on n (even) strands labelled from 0 to n-1, with k layers labelled from 0 to k-1,
  //are represented by non-zero integers:
  //     [a_k-1,1] [a_k-1,3]  [a_k-1,5] ...[a_k-1,n-3]  
  //                      ....
  //      [a_2,1]   [a_2,3]   [a_2,5]  ...  [a_2,n-3]
  //  [a_1,0]   [a_1,2]   [a_1,4]   [a_1,6] ...    [a_1,n-2]
  //      [a_0,1]   [a_0,3]   [a_0,5]  ...  [a_0,n-3]
  //
  // where    [a_i,j] represents (s_j)^{a_i,j} at layer i, a_i,j positive or negative
  // 

  //In the Markov process, each layer correspond to a layer of the plat braid, i.e.
  // with t the max number of twists (0 included), -t <= a_i,j <= +t, there are
  // (2t+1)^{n/2-1} states for even layers 0,2,4 ..., and,
  // (2t+1)^{n/2} states for odd layers 1,3,5 etc 

  //call T the max number of twists; for any twist a_i,j, we have -T <= a_i,j <= T,
  //and (a_i,j + T) \in [0;2T]
  // consider layer 
  //      [a_i,1]   [a_i,3]   [a_i,5]  ...  [a_i,n-3]
  // we give it label:
  // 
  //  (a_i,1 + T) + (a_i,3 + T)*(2T+1) + (a_i,5 + T)*(2T+1)^2 + ...
  //  from 0 to (2T+1)^{n/2-1}
  // 
  // consider layer 
  //  [a_i,0]   [a_i,2]   [a_i,4]   [a_i,6] ...    [a_i,n-2]
  // we give it label:
  //
  //  (a_i,0 + T) + (a_i,2 + T)*(2T+1) + (a_i,4 + T)*(2T+1)^2 + ...
  // from 0 to (2T+1)^{n/2}

  //NB, size_t takes values in 0 ... 4,294,967,295
  // e.g., T=20, n=10 -> 41^5 <        116.000.000

//create the markov process with appropriate number of nodes in each layer
  std::cout << "Initialize the Markov process with layers of size: ";
  std::vector<size_t> num_states(num_layers);
  for(int lay=0; lay<num_layers; ++lay) {
    if(lay % 2) { 
      num_states[lay] = std::pow((2*max_twists+1),num_strands/2);
     }
    else { 
      num_states[lay] = std::pow((2*max_twists+1),num_strands/2-1);
    }
  }
  for(auto s : num_states) { std::cout << s << " "; }
  std::cout << "\n";
  Markov_decision mdp(num_states);
  std::cout << "                done.\n";
  



//compute the braid from a given state
  auto state_to_braid = [&](std::vector<size_t> &s)->Plat_braid {
    
    // std::cout << "   state_to_braid - state = ";
    // for(auto i : s) { 
    //   std::cout << i << " ";
    // }
    // std::cout << "\n";

    Plat_braid b(num_strands);
    for(int layer_idx=0; layer_idx<(int)s.size(); ++layer_idx) {
      int x = s[layer_idx];
      int idx_twist;

      // std::cout << "   ---- layer = " << layer_idx << "\n";
      if(layer_idx % 2 == 0) {//even layer [a_i,1]   [a_i,3]   [a_i,5]  ...  [a_i,n-3]
        idx_twist = 1;
      }
      else {//odd layer [a_i,0]   [a_i,2]   [a_i,4] 
        idx_twist = 0;
      }
      while(x != 0) {
        int rem = x % (2*max_twists+1);//== a_i,j+T

        // std::cout << "    add_twist idx = " << idx_twist << "  num twists = " << rem-max_twists << "\n";
        b.add_twist(idx_twist,rem - max_twists);
        x /= (2*max_twists+1);
        // std::cout << "      after division x == " << x << "\n";
        idx_twist +=2 ;
      }
      // std::cout << "       done layer " << layer_idx << "\n";
    }
    // std::cout << "     done state to braid\n";
    return b;
  };

//compute the corresponding labels for a given braid
  auto braid_to_state = [&](const Plat_braid &b)->std::vector<size_t> {
    std::vector<size_t> state; state.reserve(num_layers);
    int layer_idx = 0;
    std::vector< int > twist_regions(num_strands-1,0);//twist_regions[j] = a_i,j, when considering layer of index i

    for(auto &layer : b.braid()) {
      size_t lab = 0;
      //set all the a_k,i+T
      for(auto pp : layer) {
        twist_regions[pp.first] = pp.second + max_twists;//in [0,2T]
      }
      if(layer_idx % 2) {//odd layer (a_i,0 + T) + (a_i,2 + T)*(2T+1) + ...
        int pow_num = 0;
        for(int j=0; j<num_strands-1; j+=2) {
          lab += twist_regions[j] * (std::pow(2*max_twists+1,pow_num));
          ++pow_num;//==j/2
        }
      }
      else {//even layer
        int pow_num = 0;
        for(int j=1; j<num_strands-2; j+=2) {
          lab += twist_regions[j] * (std::pow(2*max_twists+1,pow_num));
          ++pow_num;//==(j-1)/2
        }
      }
      state.push_back(lab);
      ++layer_idx;//next layer
      //put all a_k,i to 0
      for(size_t i=0; i<twist_regions.size(); ++i) {
        twist_regions[i] = 0;
      }
    }
    return state;
  };


  std::cout << "Initialize the quantum data:\n";

//precompute matrices and quantum data
  Jones_polynomial J(num_strands,num_layers*max_twists);

  std::cout << "                done.\n";

  // //high reward for jones polynomials close to the trivial Jones,
  // auto rew_jones = [&](std::vector<size_t> &samp)->P_float
  // {
  //   Plat_braid b = state_to_braid(samp);
  //   auto num_comp = b.num_components();
  //   if(num_comp > 1) { return 0; }
  //   auto jones_poly = J.quantum_invariant(b);
  //   int norm = jones_poly.sq_norm();
  //   auto num_cross = b.num_crossings();

  //   P_float to_decrease = (P_float)(norm) / (P_float)(sum_n_squared); //approximately between (0;1], closer to 0
  //   P_float to_increase = (P_float)(num_cross)/(P_float)(max_num_crossings);//approx between (0;1] 
  //   P_float rew = to_increase/to_decrease;

  //   // if(num_comp > 1) { rew *= -1; }

  //   std::cout << "------------------------------------------------------------------- #components = " << num_comp << "  #crossings = " << num_cross << "      reward = " << rew << "    sq_norm = " << norm << "\n";
  //   std::cout << b.oriented_Gauss_code() << "\n";
  //   return rew;
  // };


  time_t now = time(0);
  // Convert time to tm structure
  tm *local_time = localtime(&now);
  // Define format string
  char buffer[80];
  strftime(buffer, 80, "%Y-%m-%d_%H:%M:%S", local_time);

  std::string filename = "~/git_repo/kumquat/data/numstrands_" + std::to_string(num_strands) + "_numlayers_ " + std::to_string(num_layers) + "_maxtwists_" + std::to_string(max_twists) + "__" + buffer + ".dat";

  std::ofstream outfile(filename);

  //high reward for jones polynomials close to the trivial Jones,
  auto rew_jones = [&](std::vector<size_t> &samp)->P_float
  {
    Plat_braid b = state_to_braid(samp);
    auto num_comp = b.num_components();
    if(num_comp > 1) { return 0; }
    auto jones_poly = J.quantum_invariant(b);
    int sprd = jones_poly.spread();
    auto num_cross = b.num_crossings();

    std::cout << "------------------------------------------------------- #components = " << num_comp << "  #crossings = " << num_cross << "  spread = " << sprd << "\n";
//"  reward = " << rew << 

    // P_float to_decrease = (P_float)(norm) / (P_float)(sum_n_squared); //approximately between (0;1], closer to 0
    // P_float to_increase = (P_float)(num_cross)/(P_float)(max_num_crossings);//approx between (0;1] 
    P_float rew = (P_float)(1) / (P_float)(sprd);

    // if(num_comp > 1) { rew *= -1; }

    outfile << "--- #components = " << num_comp << "  #crossings = " << num_cross << "  reward = " << rew << "  spread = " << sprd << "\n";
    outfile << b.oriented_Gauss_code() << "\n";

    // std::cout << b.oriented_Gauss_code() << "\n";
    return rew;
  };






  //evaluate randomly the polynomial reward to later normalize it
  // std::cout << "Out of 10 rewards: ----------------------------\n\n";
  // P_float minr = 0;
  // P_float maxr = 0; 
  // for(size_t k=0; k<10; ++k) {
  //   auto samp = mdp.sample();
  //   auto tmpr = rew_jones(samp);
  //   if(tmpr > maxr) { maxr = tmpr; }
  //   if(tmpr < minr) { minr = tmpr; }
  //   std::cout << "       -------------- end of ite " << k << "/100   with rew_jones = " << tmpr << "        min rew = " << minr << "  max rew = " << maxr << "\n";
  // }
  // std::cout << "\n" << "Min reward = " << minr << "  Max reward = " << maxr << "\n";


  //normalized polynomial reward
  // auto rew_jones_balanced = [&](std::vector<size_t> &samp)->P_float
  // {
  //   auto rew = rew_jones(samp);
  //   rew -= (maxr+minr)/2.;
  //   rew *= 0.01 * (1./(maxr-minr));
  //   return rew;
  // };



  //train until possibly stuck in a state
  size_t num_iteration = (size_t)num_iterations;//generally converges
  auto max_samp = mdp.train(num_iteration, rew_jones);

  //sample the state and brute force improve it a few times
  std::vector<size_t> curr_samp;
  auto imp_samp = max_samp;

  size_t count_turns=0;
  do {
    std::cout << "======================================================= Turn #" << count_turns++ << "\n";
    curr_samp = imp_samp;
    imp_samp = mdp.local_improve(curr_samp, rew_jones);
  } while(curr_samp != imp_samp);

  // auto samp = mdp.sample();
  // auto imp_samp = mdp.local_improve(max_samp, rew_poly);

  std::cout << "\n\n\n";
  std::cout << "Sample          = " << mdp.to_string(max_samp) << "of rew == " << rew_jones(max_samp) <<  "\n";
  std::cout << "Improved sample = " << mdp.to_string(imp_samp) << "of rew == " << rew_jones(imp_samp) <<  "\n";
  std::cout << "\n";

  auto imp_braid = state_to_braid(imp_samp);
  auto jon = J.quantum_invariant(imp_braid);

  std::cout << "______________________________________________\n";
  std::cout << "Maximizing braid is: \n" << imp_braid << "\n";

  std::cout << "  with polynomial = " << jon << "      vs J(O) = " << J.quantum_invariant_unknot() << "\n";

  outfile.close();


return 0;




}
















