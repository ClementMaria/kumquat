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
  if(argc != 8) { 
    std::cout << "./exe num_strands num_layers max_twists num_iterations batch_size num_power_trainings alpha_reinitialize\n";
    return 0;
  }
  int num_strands = std::atoi(argv[1]);//must be even
  int num_layers = std::atoi(argv[2]);//better if odd
  int max_twists = std::atoi(argv[3]);
  int num_iterations = std::atoi(argv[4]);
  int batch_size = std::atoi(argv[5]);
  int num_pow_train = std::atoi(argv[6]);
  P_float alpha = std::atof(argv[7]);//0.5;

  // std::vector<size_t> starting_state = {7, 86, 0, 2, 6, 42, 6, 129, 7, 98, 1, 7, 7, 148, 6, 9, 12, 37, 6, 65, 8, 86, 0, 30, 0, 26, 0, 124, 0, 154, 3, 74, 6, 83, 4, 167, 0, 21, 7, 127, 6, 70, 5, 49, 0, 24, 6, 107, 5, 71, 11, 94, 5, 70, 11, 97, 7, 80, 0, 113, 6, 110, 0, 1, 6, 23, 6, 15, 2, 66, 0, 141, 0, 39, 0, 127, 7, 34, 0, 59, 0, 94, 8, 58, 11, 6, 8, 49, 6, 106, 5, 45, 1, 9, 7, 98, 7, 52, 11, 132, 4};


  // std::vector<size_t> starting_state = {7, 8, 0, 2, 0, 42, 0, 129, 7, 98, 3, 7, 7, 96, 0, 87, 0, 37, 0, 65, 2, 80, 12, 36, 0, 26, 0, 124, 6, 150, 5, 68, 0, 85, 0, 167, 6, 21, 1, 97, 0, 70, 11, 51, 6, 22, 6, 109, 5, 71, 11, 94, 5, 70, 11, 97, 7, 80, 6, 115, 0, 108, 6, 79, 0, 15, 0, 21, 2, 118, 0, 9, 0, 19, 0, 125, 9, 92, 0, 83, 0, 96, 6, 58, 11, 162, 0, 25, 6, 132, 7, 97, 3, 113, 7, 98, 7, 56, 5, 112, 4 };

  // std::vector<size_t> starting_state = {7, 8, 0, 2, 0, 42, 0, 129, 7, 98, 5, 7, 7, 96, 0, 85, 6, 35, 0, 67, 4, 6, 2, 36, 0, 26, 0, 124, 0, 144, 7, 72, 0, 5, 0, 165, 0, 23, 1, 97, 0, 70, 11, 51, 0, 22, 0, 109, 5, 71, 11, 94, 5, 70, 11, 97, 7, 80, 0, 115, 6, 108, 6, 1, 6, 15, 0, 21, 2, 92, 6, 9, 0, 123, 0, 21, 9, 92, 0, 83, 6, 96, 6, 58, 11, 160, 6, 25, 6, 134, 7, 97, 3, 113, 7, 98, 7, 56, 5, 112, 4 };

  std::vector<size_t> starting_state = {7, 8, 0, 2, 0, 42, 0, 129, 7, 98, 5, 7, 7, 96, 0, 7, 0, 35, 0, 67, 2, 0, 2, 32, 0, 38, 6, 118, 0, 144, 11, 74, 0, 3, 0, 167, 6, 23, 1, 97, 6, 70, 11, 51, 0, 22, 0, 109, 5, 71, 11, 94, 5, 70, 11, 97, 7, 2, 0, 115, 6, 108, 6, 79, 0, 15, 0, 21, 2, 118, 0, 35, 6, 97, 0, 73, 9, 92, 6, 83, 6, 96, 6, 58, 11, 160, 6, 25, 0, 134, 7, 97, 3, 113, 7, 98, 7, 56, 5, 112, 4 };

  std::cout << "Number of strands:          " << num_strands << "\n";
  std::cout << "Number of layers:           " << num_layers << "\n";
  std::cout << "Max number of twists:       " << max_twists << "\n";
  std::cout << "Number iterations (train)   " << num_iterations << "\n";
  std::cout << "Batch size (train):         " << batch_size << "\n";
  std::cout << "Number of power trainings:  " << num_pow_train  << "\n";
  std::cout << "Alpha for reinitialization: " << alpha << "\n";



/****************************************************************/
  auto start_mdp = std::chrono::high_resolution_clock::now();
  std::cout << "*** Initializes Markov_decision_quantum_knot\n";
//
  Markov_decision_quantum_knot mdp(num_strands, num_layers, max_twists);
//
  auto end_mdp = std::chrono::high_resolution_clock::now();
  auto duration_mdp = duration_cast<std::chrono::seconds>(end_mdp-start_mdp);  
  std::cout << "    done in " << duration_mdp.count() << " sec.\n";
/****************************************************************/
  //train until possibly stuck in a state


  std::cout << "*** Power training\n";

  auto best_state = mdp.power_train( num_pow_train
                                   , num_iterations
                                   , batch_size
                                   , alpha
                                   // , starting_state // 
                                   );



  std::cout << "Finishes with best_state = ";
  for(auto i : best_state) { std::cout << i << " "; }
  std::cout << std::endl;

  std::cout << mdp.state_to_braid(best_state).oriented_Gauss_code() << "\n";
  
















  return 0;

//   size_t num_iteration = (size_t)num_iterations;

//   auto start_tr = std::chrono::high_resolution_clock::now();  
//   std::cout << "*** Train Markov process\n";
// //
//   auto max_samp = mdp.train(num_iteration, batch_size);
// //
//   auto end_tr = std::chrono::high_resolution_clock::now();
//   auto duration_tr = duration_cast<std::chrono::seconds>(end_tr-start_tr);  
//   std::cout << "    done in " << duration_tr.count() << " sec.\n";
// /****************************************************************/
//   auto start_ie = std::chrono::high_resolution_clock::now();  
//   std::cout << "*** Improve exhaustive\n";
// //
//   auto imp_samp = mdp.improve_exhaustive(max_samp);
// //
//   auto end_ie = std::chrono::high_resolution_clock::now();
//   auto duration_ie = duration_cast<std::chrono::seconds>(end_ie-start_ie);  
//   std::cout << "    done in " << duration_ie.count() << " sec.\n";
// /****************************************************************/
// //open file for output
//   time_t now = time(0);  // Convert time to tm structure
//   tm *local_time = localtime(&now);
//   char buffer[80];  // Define format string
//   strftime(buffer, 80, "%Y-%m-%d_%H:%M:%S", local_time);
//   std::string filename = "/Users/cmaria/git_repositories/kumquat/data/numstrands_" + std::to_string(num_strands) + "_numlayers_ " + std::to_string(num_layers) + "_maxtwists_" + std::to_string(max_twists) + "__" + buffer + ".dat";
//   std::ofstream outfile(filename);


//   outfile.close();

// //conclusion:
//   std::cout << "\n\n\n";
//   std::cout << "Max samp = " << mdp.to_string(max_samp) << " of rew == " << mdp.rew_jones(max_samp) <<  "\n";
//   std::cout << "Imp samp = " << mdp.to_string(imp_samp) << " of rew == " << mdp.rew_jones(imp_samp) <<  "\n";
//   std::cout << "\n";

//   // auto imp_braid = state_to_braid(imp_samp);
//   // auto jon = J.quantum_invariant(imp_braid);

//   // std::cout << "______________________________________________\n";
//   // std::cout << "Maximizing braid is: \n" << imp_braid << "\n";
//   // std::cout << "\nGauss code:\n";
//   // std::cout << imp_braid.oriented_Gauss_code();
//   // std::cout << "\n\n";
//   // std::cout << "  with polynomial = " << jon << "      vs J(O) = " << mdp.quantum_invariant_unknot() << "\n";

//   return 0;
}



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












