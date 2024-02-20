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

using namespace kumquat;

typedef Markov_decision::Proba_float P_float;

int main() {

  // std::vector<size_t> max_samp = {21, 1124, 19, 937, 20, 916, 0, 127, 37, 260, 17, 927, 19, 836, 21, 562, 23, 509, 0, 1058};

// std::vector<size_t> max_samp = {19, 1001, 25, 1091, 16, 619, 20, 989, 13, 1038, 32, 1430, 16, 467, 22, 233, 20, 1564, 0, 442};

// std::vector<size_t> max_samp = {13, 1081, 23, 1360, 20, 486, 23, 333, 30, 226, 0, 1498, 33, 1040, 26, 1075, 0, 70, 0, 1212};

// std::vector<size_t> max_samp = {15, 1338, 24, 441, 12, 7, 0, 636, 16, 19, 31, 842, 0, 848, 27, 387, 15, 684, 31, 550};

// std::vector<size_t> max_samp = {11, 1378, 0, 464, 12, 1424, 18, 676, 27, 707, 11, 1172, 17, 1098, 3, 1327, 0, 1409, 20, 187};

// std::vector<size_t> max_samp = {18, 260, 11, 1339, 0, 463, 39, 392, 10, 507, 16, 1041, 39, 687, 27, 50, 20, 1631, 25, 886};

//   int max_twists = 20;
// //10;//arbitrary
//   int num_strands = 4;//must be even
//   int num_layers = 20;//must be odd



// std::vector<size_t> max_samp = {18, 2779, 66, 4316, 3, 5234, 33, 3212, 37, 4777, 41, 6389, 76, 3173, 41, 4909, 47, 444, 19, 3459, 16, 4157, 57, 212, 10, 4131, 59, 2950, 24, 626, 43, 3546, 48, 1891, 52, 5858, 77, 3372, 36, 3353, 79, 2884, 37, 164, 12, 3547, 70, 2400, 36, 2297};
//   int num_strands = 4;//must be even
//   int num_layers = 50;//must be odd
//   int max_twists = 40;



std::vector<size_t> max_samp = {124, 2784, 151, 5543, 365, 4504, 12, 4459, 359, 197, 175, 5581, 141, 5649, 199, 2935, 274, 6927, 330, 3282};
  int num_strands = 6;//must be even
  int num_layers = 20;//must be odd
  int max_twists = 10;




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

  //high reward for jones polynomials close to the trivial Jones,
  auto rew_jones = [&](std::vector<size_t> &samp)->P_float
  {
    Plat_braid b = state_to_braid(samp);
    
    auto num_comp = b.num_components();
    if(num_comp > 1) { return 0; }

    auto jones_poly = J.quantum_invariant(b);

    int norm = jones_poly.sq_norm();

    auto num_cross = b.num_crossings();

    P_float to_decrease = (P_float)(norm) / (P_float)(sum_n_squared); //approximately between (0;1], closer to 0
    P_float to_increase = (P_float)(num_cross)/(P_float)(max_num_crossings);//approx between (0;1] 
    P_float rew = to_increase/to_decrease;

    // if(num_comp > 1) { rew *= -1; }

    std::cout << "------------------------------------------------------------------- #components = " << num_comp << "  #crossings = " << num_cross << "      reward = " << rew << "    sq_norm = " << norm << "\n";
    std::cout << b.oriented_Gauss_code() << "\n";
    return rew;
  };


  //sample the state and brute force improve it a few times
  std::vector<size_t> curr_samp;
  auto imp_samp = max_samp;
  do {
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


return 0;
}