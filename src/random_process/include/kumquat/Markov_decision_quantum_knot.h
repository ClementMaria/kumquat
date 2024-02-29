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

#ifndef KUMQUAT_MARKOV_DECISION_QUANTUM_KNOT 
#define KUMQUAT_MARKOV_DECISION_QUANTUM_KNOT

#include <chrono>
#include <ctime>
#include <chrono>
#include <fstream>
#include <set>
#include <kumquat/R.h>
#include <kumquat/Dense_matrix.h>
#include <kumquat/Markov_decision.h>
#include <kumquat/Plat_braid.h>
#include <kumquat/Jones_polynomial.h>
#include "link/link.h"//regina
#include <tbb/parallel_for.h>

namespace kumquat {

/** A decision matrix for a Markov decision process, where states are in bijections with plat braids of a certain topology (number of stransd, number of layers, maximum number of twists).
 * 
 plat braids on n (even) strands labelled from 0 to n-1, with k layers labelled from 0 to k-1,
are represented by non-zero integers:
    [a_k-1,1] [a_k-1,3]  [a_k-1,5] ...[a_k-1,n-3]  
                     ....
     [a_2,1]   [a_2,3]   [a_2,5]  ...  [a_2,n-3]
 [a_1,0]   [a_1,2]   [a_1,4]   [a_1,6] ...    [a_1,n-2]
     [a_0,1]   [a_0,3]   [a_0,5]  ...  [a_0,n-3]

where    [a_i,j] represents (s_j)^{a_i,j} at layer i, a_i,j positive or negative


In the Markov process, each layer correspond to a layer of the plat braid, i.e.
with t the max number of twists (0 included), -t <= a_i,j <= +t, there are
(2t+1)^{n/2-1} states for even layers 0,2,4 ..., and,
(2t+1)^{n/2} states for odd layers 1,3,5 etc 

call T the max number of twists; for any twist a_i,j, we have -T <= a_i,j <= T,
and (a_i,j + T) \in [0;2T]
consider layer 
     [a_i,1]   [a_i,3]   [a_i,5]  ...  [a_i,n-3]
we give it label:

 (a_i,1 + T) + (a_i,3 + T)*(2T+1) + (a_i,5 + T)*(2T+1)^2 + ...
 from 0 to (2T+1)^{n/2-1}

consider layer 
 [a_i,0]   [a_i,2]   [a_i,4]   [a_i,6] ...    [a_i,n-2]
we give it label:

 (a_i,0 + T) + (a_i,2 + T)*(2T+1) + (a_i,4 + T)*(2T+1)^2 + ...
from 0 to (2T+1)^{n/2}

NB, size_t takes values in 0 ... 4,294,967,295
e.g., T=20, n=10 -> 41^5 <        116.000.000
 * 
 * 
 */
class Markov_decision_quantum_knot {
public:
  /** A float type to represent probabilities.**/
  typedef typename Markov_decision::Proba_float Proba_float;

  /** \brief Initialize the transition matrices with uniform 
   * probabilities. 
   *
   * @param[in] num_states the vector of sizes of the sets of nodes 
   *                       num_states[i] = |Vi| of the multipartite 
   *                       graph.
   **/
  Markov_decision_quantum_knot(int num_strands, int num_layers, int max_twists) : 
  num_strands_(num_strands),
  num_layers_(num_layers),
  max_twists_(max_twists)
  {
    //compute the maximum number of crossings
    int max_num_crossings = 0;
    if(num_layers % 2 == 0) {
      max_num_crossings = max_twists * ( (num_layers/2) * (num_strands/2) + (num_layers/2 + 1)*(num_layers/2-1) );
    }
    else {
      max_num_crossings = max_twists * ( (num_layers/2) * (num_strands/2) + (num_layers/2)*(num_layers/2-1) );
    }
    // int sum_n_squared = max_num_crossings * (max_num_crossings-1)/2;
    std::cout << "---------- Nb. strands = " << num_strands << "  Nb. layers = " << num_layers << "  Max nb. twists = " << max_twists << " ----------\n";
    std::cout << "---> max # crossings = " << max_num_crossings << "\n";
    //create the markov process with appropriate number of nodes in each layer
    auto start_mp = std::chrono::high_resolution_clock::now();
    std::cout << "Initialize the Markov process with layers of size: ";
    //compute the number of nodes per layer
    std::vector<size_t> num_states(num_layers);
    for(int lay=0; lay<num_layers; ++lay) {
      if(lay % 2) { 
        num_states[lay] = std::pow((2*max_twists+1),num_strands/2);
       }
      else { 
        num_states[lay] = std::pow((2*max_twists+1),num_strands/2-1);
      }
    }
    //display the size of each layer
    for(auto s : num_states) { std::cout << s << " "; }
    std::cout << "\n";
    //initialize the transition matrices
    markov_dec_ = Markov_decision(num_states);
    //
    auto end_mp = std::chrono::high_resolution_clock::now();
    auto duration_mp = duration_cast<std::chrono::seconds>(end_mp-start_mp);  
    std::cout << "--- in " << duration_mp.count() << " sec.\n";
    //initialize the quantum data
    auto start_qd = std::chrono::high_resolution_clock::now();
    std::cout << "Initialize the quantum data:\n";
    //precompute matrices and quantum data (jones polynomial)
    jones_poly_ = Jones_polynomial(num_strands,num_layers*max_twists);
    //
    auto end_qd = std::chrono::high_resolution_clock::now();
    auto duration_qd = duration_cast<std::chrono::seconds>(end_qd-start_qd);  
    std::cout << "--- in " << duration_qd.count() << " sec.\n";


  //open file for output
    // time_t now = time(0);  // Convert time to tm structure
    // tm *local_time = localtime(&now);
    // char buffer[80];  // Define format string
    // strftime(buffer, 80, "%Y-%m-%d_%H:%M:%S", local_time);
    // std::string filename = "/Users/cmaria/git_repositories/kumquat/data/numstrands_" + std::to_string(num_strands) + "_numlayers_ " + std::to_string(num_layers) + "_maxtwists_" + std::to_string(max_twists) + "__" + buffer + ".dat";
    // out_ = std::ofstream(filename);

  }

  /** \brief Return a state sampled from the distribution of the Markov process. 
   * @param[out] a vector of integers, representing the label of a node in V0, V1 ... Vt.
   * **/
  std::vector< size_t > sample() {
    return markov_dec_.sample();
  }

  /** \brief Increase or decrease the probabilities of transitions 
   * corresponding to a certfain state.
   * 
   * @param[in] samp the state for which the transitions are amplified,
   * @param[in] alpha an float between -1 < alpha < 1.
   * */
  void amplify(std::vector<size_t> samp, Proba_float alpha) {
    markov_dec_.amplify(samp,alpha);
  }

/** \brief For a given state, try, for each layer, to change a node locally to see if it improves the reward.*/
  std::vector<size_t> local_improve(std::vector<size_t> &state) 
  { //starting state
    std::vector<size_t> curr_state(state);
    //for all layers in the braid/Markov process
    for(int layer_idx=0; layer_idx < state.size(); ++layer_idx) {
      
      std::cout << "   " << layer_idx << "/" << state.size() << "\n";

      Plat_braid b_start = state_to_braid(curr_state,false);//get the braid corresponding to curr_state
      //we are modifying layer 'layer_idx', compute the morphisms for layers [0;layer_idx-1] and [layer_idx+1, top]*h_tensor
      auto Mb = jones_poly_.quantum_morphism(b_start.braid().begin(), b_start.braid().begin()+layer_idx);

      auto Mt = jones_poly_.quantum_morphism(b_start.braid().begin()+(layer_idx+1), b_start.braid().end());

      Mt.ltimes_equal(jones_poly_.h_tensor());
      Mb.rtimes_equal(Mt);//tr(Mt M Mb) = tr(M Mb Mt)

      //in parallel, look for best label in layer 'layer_idx'
      size_t size_layer = markov_dec_.size_layer(layer_idx);
      //thread i computes the reward for label i in layer layer_idx, and writes it in thread_rew[i]
      std::vector< Proba_float > thread_rew(size_layer,-1);

      tbb::parallel_for(size_t(0), size_t(size_layer),
        [&](size_t i) {
          //copy the current state and modify its layer 'layer_idx'
          std::vector<size_t> curr_state_thread(curr_state);
          curr_state_thread[layer_idx] = i;//the state in this thread
          //turn it into a braid (without simplifying)
          Plat_braid b_thread = state_to_braid(curr_state_thread, false);        
          //multiply by hand the matrices to get the Jones polynomial
          Proba_float rew;//rew is set to 0 if b_thread is not a knot
          if(b_thread.num_components() > 1) { rew = 0; }
          else {//we have a knot:
            Jones_polynomial::Morphism M = jones_poly_.quantum_morphism(
                             b_thread.braid().begin()+layer_idx, 
                             b_thread.braid().begin()+(layer_idx+1));
            M.rtimes_equal(Mb);//reuse product Mb h Mt
            auto q_inv = jones_poly_.trace(M);//jones poly of b_thread
            int sprd = q_inv.spread();


            bool flag_unknot = false;
            if(sprd == 2 && (q_inv == jones_poly_.quantum_invariant_unknot())) {

              std::string oGc = b_thread.oriented_Gauss_code();

              //do our best to identify the unknot
              flag_unknot = unknot(oGc); 

              if(!flag_unknot) {
                std::cout << "difficult Jones unknot: " << oGc << "\n";
                candidates_cex_.insert(oGc);
              }

            }

            if(flag_unknot) { rew = 0.; }
            else {
              //simplify with regina to get a better approximation of the number of crossings
              auto gausscode = b_thread.oriented_Gauss_code();
              regina::Link regina_link(gausscode);
              regina_link.intelligentSimplify();
              auto num_cross = regina_link.crossings().size();
              //compute the reward knowing #crossing and spread(jones)
              rew = rew_formula(num_cross,sprd);
            }
          }
          thread_rew[i] = rew;
        });
      //extract the best reward and idx
      size_t best_label = state[layer_idx];
      Proba_float best_reward = -1;
      for(int j=0; j<size_layer; ++j) {
        if(thread_rew[j] > best_reward) {
          best_reward = thread_rew[j];
          best_label = j;
        }
      }
      //upgrade the current state and go to next layer
      curr_state[layer_idx] = best_label;
    }
    return curr_state;
  }

/**
 * @brief From an input state, try to iteratively modify the node layer after layer, and keep greedily the best node in each layer. Work on each layer iteratively, and loop until no imporvement is possible.
 * @details Recycle matrices when working on a given layer.
 * 
 * @param state A state to improve on
 * @return A local maximum, i.e., not local improvement is possible.
 */
  std::vector<size_t> improve_exhaustive(std::vector<size_t> &state)
  {//sample the state and brute force improve it a few times
    std::vector<size_t> curr_samp;
    auto imp_samp = state;
    size_t count_turns=0;//count the number of rounds
    do {
      auto start_li = std::chrono::high_resolution_clock::now();
      std::cout << "   - improve_exhaustive - turn #" << count_turns++ << "\n";
      curr_samp = imp_samp;
      //
      imp_samp = local_improve(curr_samp);
      //
      auto end_li = std::chrono::high_resolution_clock::now();
      auto duration_li = duration_cast<std::chrono::seconds>(end_li-start_li);  
      std::cout << "--- in " << duration_li.count() << " sec. \n";
      std::cout << "{";
      for(auto s : imp_samp) { std::cout << s << ", "; }
      std::cout << "}" << std::endl;

      write_file();

    } while(curr_samp != imp_samp);

    return imp_samp;
  }


/** \brief Turn a state into a string.*/
  std::string to_string(std::vector<size_t> &state) {
    return markov_dec_.to_string(state);
  }

  /** \brief Train the Markov decision process by sampling a state, compute its reward and turn it into an amplifying value alpha, and amplify the probabilities of transtion.
   * 
   * @param[in] num_iterations the number of runs for the training,
   * @param[in] reward a reward function to apply on a state ; it returns values in (-1,1) that can be given to the amplify method.
   **/
  std::vector<size_t> train(size_t num_iterations, size_t batch_size = 100) 
  {
    auto reward = [&](std::vector<size_t> &state)->Proba_float {
      return rew_jones(state);
    };
    return markov_dec_.train(num_iterations,reward, batch_size);
  }

/**
 * @brief Meta training process.
 * @details Train the Markov process, extract the best encountered state, local improve it, reinitialize the Markov chain with the best state, train again and so on.
 * 
 * @param num_power_trainings Number of overall train-local_improve-reinitialize iterations.
 * @param num_ite_training Number of iterations in each training
 * @param batch_size_training Size of batches for parallelized training
 * @return The best state encountered overall.
 */
  std::vector<size_t> power_train(size_t num_power_trainings, size_t num_ite_training, size_t batch_size_training, Proba_float alpha, std::vector<size_t> starting_state = std::vector<size_t>()) {
    
    // //open file for output
    time_t now = time(0);  // Convert time to tm structure
    tm *local_time = localtime(&now);
    char buffer[80];  // Define format string
    strftime(buffer, 80, "%Y-%m-%d_%H:%M:%S", local_time);
    outfilename_ = "/Users/cmaria/git_repositories/kumquat/data/counter_examplesJO__" + std::string(buffer) + ".dat";

    flag_cex_ = false;
    //the reward function
    auto reward = [&](std::vector<size_t> &state)->Proba_float {
      return rew_jones(state);
    };

    std::vector<size_t> best_state_overall;
    Proba_float best_reward_overall = -1;

    std::vector<size_t> best_state_training, best_state_improve;
    
    if(!starting_state.empty()) {
      //improve exhaustive this state
      std::cout << " - improve starting_state\n";
      best_state_improve = improve_exhaustive(starting_state);
      auto rew = reward(best_state_improve);
      if(rew > best_reward_overall) {
        best_reward_overall = rew;
        best_state_overall = best_state_improve;
      }
      std::cout << " - initialize starting_state_improve\n";
      markov_dec_.initialize(best_state_improve, alpha);
    }


    do {
      std::cout << "Power training left: " << num_power_trainings << "        [" << flag_cex_ << "]\n";
      auto start = std::chrono::high_resolution_clock::now();

      //train mdp and keep best state encountered
      
      std::cout << " - train\n";
      best_state_training = markov_dec_.train(num_ite_training, reward, batch_size_training);

      write_file();

      std::cout << " - improve\n";

      //improve exhaustive this state
      best_state_improve = improve_exhaustive(best_state_training);
     
      std::cout << " - reward\n";

      auto rew = reward(best_state_improve);
      if(rew > best_reward_overall) {
        best_reward_overall = rew;
        best_state_overall = best_state_improve;
      }

      // std::cout << " - flag_cex\n";
      // if(flag_cex_) {
      //   for(auto cex : counter_examples_) {
      //     out_ << cex << "\n";
      //   }
      //   out_ << "\n";
      // }
      auto end = std::chrono::high_resolution_clock::now();
      auto duration = duration_cast<std::chrono::seconds>(end-start);  

      std::cout << "    - best rew|state: [" << best_reward_overall << "]    {";
      for(auto s : best_state_overall) { std::cout << s << ", "; }
      std::cout << "}         in " << duration.count() << " sec.\n";

      std::cout << " - initialize\n";
      //reinitialize the mdp with the improved state
      markov_dec_.initialize(best_state_improve, alpha);

    } while(--num_power_trainings != 0);

  // out_.close();

  return best_state_overall;
}


/** \brief Return the number of layers in the Markov decision process.*/
  size_t num_layers() { return markov_dec_.num_layers(); }
/** Return the size of a layer given by its index.*/
  size_t size_layer(size_t i) {
    return markov_dec_.size_layer(i);
  }













  /**
   * @brief Compute a state from a braid.
   * @details 
   * 
   * @param b A plat braid.
   * @return A state as a list of indices for each layer.
   */
  std::vector<size_t> braid_to_state(const Plat_braid &b) 
  {
    std::vector<size_t> state; state.reserve(num_layers_);
    int layer_idx = 0;
    std::vector< int > twist_regions(num_strands_-1,0);//twist_regions[j] = a_i,j, when considering layer of index i

    for(auto &layer : b.braid()) {
      size_t lab = 0;
      //set all the a_k,i+T
      for(auto pp : layer) {
        twist_regions[pp.first] = pp.second + max_twists_;//in [0,2T]
      }
      if(layer_idx % 2) {//odd layer (a_i,0 + T) + (a_i,2 + T)*(2T+1) + ...
        int pow_num = 0;
        for(int j=0; j<num_strands_-1; j+=2) {
          lab += twist_regions[j] * (std::pow(2*max_twists_+1,pow_num));
          ++pow_num;//==j/2
        }
      }
      else {//even layer
        int pow_num = 0;
        for(int j=1; j<num_strands_-2; j+=2) {
          lab += twist_regions[j] * (std::pow(2*max_twists_+1,pow_num));
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
  }


/**
 * @brief Compute a braid from a state.
 * @details 
 * 
 * @param s A state
 * @return The corresponding braid
 */
    Plat_braid state_to_braid(const std::vector<size_t> &s, bool simplify = true) {
    Plat_braid b(num_strands_);
    if(simplify) {
      for(int layer_idx=0; layer_idx<(int)s.size(); ++layer_idx) {
        int x = s[layer_idx];
        int idx_twist;
        //even layer [a_i,1]   [a_i,3]   [a_i,5]  ...  [a_i,n-3]
        if(layer_idx % 2 == 0) { idx_twist = 1; }
        //odd layer [a_i,0]   [a_i,2]   [a_i,4]
        else { idx_twist = 0; }
        while(x != 0) {
          int rem = x % (2*max_twists_+1);//== a_i,j+T
          b.add_twist(idx_twist,rem - max_twists_);
          x /= (2*max_twists_+1);
          idx_twist +=2 ;
        }
      }
      b.greedy_simplify();//simplify the braid with Markov moves
      return b;
    }
    else {//do not simplify to match the structure of the markov process
      for(int layer_idx=0; layer_idx<(int)s.size(); ++layer_idx) {
        b.push_layer();
        int x = s[layer_idx];
        int idx_twist;
        //even layer [a_i,1]   [a_i,3]   [a_i,5]  ...  [a_i,n-3]
        if(layer_idx % 2 == 0) { idx_twist = 1; }
        //odd layer [a_i,0]   [a_i,2]   [a_i,4]
        else { idx_twist = 0; }
        while(x != 0) {
          int rem = x % (2*max_twists_+1);//== a_i,j+T
          b.add_twist_idx(idx_twist,rem - max_twists_, layer_idx);
          x /= (2*max_twists_+1);
          idx_twist +=2 ;
        }
      }
      return b;
    }
  }

  // std::map<int,int> state_to_braid_layer(size_t label, int layer_idx) {
  //   std::map<int,int> layer;

  //   int idx_twist;
  //   //even layer [a_i,1]   [a_i,3]   [a_i,5]  ...  [a_i,n-3]
  //   if(layer_idx % 2 == 0) { idx_twist = 1; }
  //   //odd layer [a_i,0]   [a_i,2]   [a_i,4]
  //   else { idx_twist = 0; }
  //   while(label != 0) {
  //     int rem = label % (2*max_twists_+1);//== a_i,j+T
  //     if(rem - max_twists_ != 0 ) {
  //       layer[idx_twist] = rem - max_twists_;
  //     }
  //     label /= (2*max_twists_+1);
  //     idx_twist +=2 ;
  //   }
  //   return layer;
  // }



  Proba_float rew_formula(int num_crossings, int spread_jones) {
    return (Proba_float)(num_crossings) / (Proba_float)(spread_jones*spread_jones);
  }
/**
 * @brief Compute the reward for an input state, using Kumquat's own implementation.
 * @details Turn the state into a braid, and compute its Jones polynomial. The reward is the tradeoff spread(Jones poly)/number crossings after simplification (convert to Regina, simplify, then output number of crossings).
 * 
 * @param samp A state
 * @return A reward value
 */
   Proba_float rew_jones(std::vector<size_t> &samp)
  {
    Plat_braid b = state_to_braid(samp);
    b.greedy_simplify();
    auto num_comp = b.num_components();
    if(num_comp > 1) { return 0.; }//not connected => 0.
    //compute the jones polynomial
    auto jones_poly = jones_poly_.quantum_invariant(b);
    int sprd = jones_poly.spread();
    //check if jones poly is trivial
    if(sprd == 2 && (jones_poly == jones_poly_.quantum_invariant_unknot())) {

      std::string oGc = b.oriented_Gauss_code();

      //do our best to identify the unknot
      bool flag_unknot = unknot(oGc); 

      if(flag_unknot) { return 0.; }//identified the unknot
      //otherwise, candidate counter example
      candidates_cex_.insert(oGc);
      std::cout << "difficult Jones unknot: " << oGc << "\n";
    }
    //compute number of crossings after simplification
    auto num_cross_pbraid = b.num_crossings(); //<- naive
    auto gausscode = b.oriented_Gauss_code();
    regina::Link L(gausscode);
    L.intelligentSimplify();
    // auto jones_poly = L.jones();
    // int sprd = jones_poly.maxExp() - jones_poly.minExp();
    auto num_cross = L.crossings().size();
    // std::cout << "------------- #components = " << num_comp << "  #crossings = " << num_cross_pbraid << "|" << num_cross << "  spread = " << sprd << "      spread/num_crossings = " << (double)(sprd) / (double)(num_cross) << "\n";
    Proba_float rew = rew_formula(num_cross,sprd);
    return rew;
  };





/**
 * @brief Compute the reward for an input state, using Regina.
 * @details Turn the state into a braid, convert it into a Regina knot, simplify it and compute its Jones polynomial. The reward is the tradeoff spread(Jones poly)/number crossings after simplification.
 * 
 * @param samp A state
 * @return A reward value
 */
  // Proba_float rew_jones_regina(std::vector<size_t> &samp)
  // {
  //   Plat_braid b = state_to_braid(samp);
  //   auto num_comp = b.num_components();
  //   if(num_comp > 1) { return 0; }

  //   auto gausscode = b.oriented_Gauss_code();
  //   regina::Link L(gausscode);

  //   auto jones_poly = L.jones();
  //   int sprd = jones_poly.maxExp() - jones_poly.minExp();
  //   auto num_cross = L.crossings().size();

  //   std::cout << "------------------------------------------------------- #components = " << num_comp << "  #crossings = " << num_cross << "  spread = " << sprd << "      spread/num_crossings = " << (double)(sprd) / (double)(num_cross) << "\n";

  //   Proba_float rew = 0.01*(Proba_float)(num_cross) / (Proba_float)(sprd);

  //   // outfile << "--- #components = " << num_comp << "  #crossings = " << num_cross << "  reward = " << rew << "  spread = " << sprd << "      spread/num_crossings = " << (double)(sprd) / (double)(num_cross) << "\n";
    
  //   return rew;
  // }


  /** \brief Return true if this represents the unknot. If false, this may or may not be the unknot.*/
  bool unknot(std::string oGc) {
    regina::Link K(oGc);
    
    std::cout << "    - K.intelligentSimplify()\n";
    K.intelligentSimplify();
    //do we have the trivial diagram?
    if(K.crossings().size() == 0) { return true; }
    
// bool  simplifyExhaustive (int height=1, unsigned threads=1
      // std::cout << "    - K.simplifyExhaustive(1,16)\n";
      // K.simplifyExhaustive(1,16);

    if(K.crossings().size() == 0) { return true; }

    //check if knot group is Z after simplifications
    std::cout << "    - grpK.intelligentSimplify()\n";
    auto grpK = K.group(true);
    grpK.intelligentSimplify();

    std::cout << "    - grpK.proliferateRelators()\n";
    grpK.proliferateRelators();//1

    std::cout << "    - grpK.identifyAbelian()\n";
    if(grpK.identifyAbelian()) { return true; }
    
    //check is the pi_1 of the complement (==knot group) is Z
    std::cout << "    - T.intelligentSimplify()\n";    
    auto T = K.complement();
    T.intelligentSimplify();

    // std::cout << "    - T.simplifyExhaustive(1,16)\n";    
    // T.simplifyExhaustive(1,16);

    if(T.size() <= 6) { 
      std::cout << "    - T.isSolidTorus()\n";
      return T.isSolidTorus(); 
    }

    std::cout << "    - grpT.intelligentSimplify()\n";    
    auto grpT = T.group();
    grpT.intelligentSimplify();
    
    std::cout << "    - grpT.proliferateRelators()\n";
    grpT.proliferateRelators();//1

    std::cout << "    - grpT.identifyAbelian()\n";
    if(grpT.identifyAbelian()) { return true; }


    std::cout << "candidate: " << oGc << "\n";

    //true unknot recognition is too slow
    // return T.isSolidTorus();
    return false;    

  }

  void write_file() {
    
    if( counter_examples_.empty() 
        && candidates_cex_.empty()) {
      return;
    }

    // std::ofstream out_(outfilename_);
    // out_ << "#non-trivial knots with trivial Jones polynomial.\n";
    
    std::cout << "--------------------------WRITE FILE:  #knots =  [" << candidates_cex_.size() << "]\n";
    std::cout << "#non-trivial knots with trivial Jones polynomial.\n";

    for(auto str : counter_examples_) {
      // out_ << str << "\n";
      std::cout << str << "\n\n";
    }

    // out_ << "#possibly-trivial knots with trivial Jones polynomial.\n";
    std::cout << "#possibly-trivial knots with trivial Jones polynomial.\n";

    for(auto str : candidates_cex_) {
      // out_ << str << "\n";
      std::cout << str << "\n\n";
    }
    // out_.close();

    std::cout << "--------------------------end WRITE FILE\n";

  }

  // //high reward for jones polynomials close to the trivial Jones,
  // auto rew_jones = [&](std::vector<size_t> &samp)->P_float
  // {
  //   Plat_braid b = state_to_braid(samp);
  //   auto num_comp = b.num_components();
  //   if(num_comp > 1) { return 0; }
  //   auto jones_poly = J.quantum_invariant(b);
  //   int norm = jones_poly.sq_norm();
  //   auto num_cross = b.num_crossings();

  //   Proba_float to_decrease = (Proba_float)(norm) / (Proba_float)(sum_n_squared); //approximately between (0;1], closer to 0
  //   Proba_float to_increase = (Proba_float)(num_cross)/(Proba_float)(max_num_crossings);//approx between (0;1] 
  //   Proba_float rew = to_increase/to_decrease;

  //   // if(num_comp > 1) { rew *= -1; }

  //   std::cout << "------------------------------------------------------------------- #components = " << num_comp << "  #crossings = " << num_cross << "      reward = " << rew << "    sq_norm = " << norm << "\n";
  //   std::cout << b.oriented_Gauss_code() << "\n";
  //   return rew;
  // };





// private:
public:
  Markov_decision markov_dec_;
  Jones_polynomial jones_poly_;

  int max_twists_;
  int num_strands_;
  int num_layers_;
  int num_iterations_;

  bool flag_cex_;
  std::string outfilename_;
  std::set< std::string > counter_examples_;
  std::set< std::string > candidates_cex_;
};

} //namespace kumquat

#endif // KUMQUAT_MARKOV_DECISION_QUANTUM_KNOT












