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

#ifndef KUMQUAT_MARKOV_DECISION_H_ 
#define KUMQUAT_MARKOV_DECISION_H_

#include <kumquat/R.h>
#include <kumquat/Dense_matrix.h>

namespace kumquat {

/** A decision matrix for a Markov decision process. The process is 
 * a complete 
 * multipartite graph V0 V1 ... Vt, represented by a set of transition matrices 
 * M0 ... M_t-1. 
 * 
 * Matrix Mi is a (Vi x Vi+1) dense matrix, representing the 
 * transition of local state at time i to time i+1. Precisely, 
 * Mi[s] is a row vector of length Vi+1 containing the probability 
 * distribution of transitioning to the nodes of Vi+1, knowing that 
 * we are in the s-th node of Vi. 
 * 
 * Mi[p][q] is the probability to jumping to q \in Vi+1, knowing 
 * that we are in p at time i.
 * 
 * One state of the universe is exactly a left to right path in the 
 * multipartite graph.
 * 
 * A state is a left to right path in the multipartite graph (i.e., a sequence v0 v1 ... vt of nodes in V0, V1, ... Vt respectively). The set of all states in the universe.
 **/
// template<typename NodeType>
class Markov_decision {
public:
  /** A float type to represent probabilities.**/
  typedef R::Element Proba_float;

  // typedef NodeType Node_type;

  /** \brief Initialize the transition matrices with uniform 
   * probabilities. 
   *
   * @param[in] num_states the vector of sizes of the sets of nodes 
   *                       num_states[i] = |Vi| of the multipartite 
   *                       graph.
   **/
  Markov_decision(std::vector<size_t> num_states) 
  {
    size_t num_layers = num_states.size(); 
    //idd distribution for starting point
    start_distribution_.reserve(num_states[0]);
    Proba_float unif = R_.element(1)/R_.element(num_states[0]);
    for(size_t i=0; i<num_states[0]; ++i) {
      start_distribution_.push_back(unif);
    }
    //vector containing the transitions matrices
    // decision_mat_.reserve( num_layers-1 ); 
  
    for(size_t t=0; t<num_layers-1; ++t) {
      size_t num_rows = num_states[t];
      size_t num_cols = num_states[t+1];
      decision_mat_.emplace_back(num_rows,num_cols,R_); 
      //fill the matrices with appropriately normalize uniform proba
      auto f = [&]([[maybe_unused]] size_t i, [[maybe_unused]] size_t j)->Proba_float { return (Proba_float)(1)/(Proba_float)num_cols; };   
      decision_mat_[t].fill(f);
    }
  }

  /** \brief Return a state sampled from the distribution of the Markov process. 
   * 
   * Transverse the mutipartite graph V0 ... Vt from left to right, 
   * and transition by sampling the starting and transition 
   * matrices.
   * 
   * @param[out] a vector of integers, representing the label of a node in V0, V1 ... Vt.
   * **/
  std::vector< size_t > sample() {
    std::random_device rd;//initialize random generator
    std::mt19937 gen(rd());

    std::vector<size_t> samp(num_layers());//initialize the state
    //sample a node of V0 using the starting distribution
    std::discrete_distribution<> dist_start(
             start_distribution_.begin(),start_distribution_.end());
    samp[0] = dist_start(gen);

    size_t curr_node = samp[0];
    for(size_t t=0; t<num_layers()-1; ++t) {
      //folow the transition matrices
      std::discrete_distribution<> dist_t_to_next(
                        decision_mat_[t][curr_node].begin(), 
                        decision_mat_[t][curr_node].end());
      curr_node = dist_t_to_next(gen);
      samp[t+1] = curr_node;
    }
    return samp;
  }

  /** \brief Increase or decrease the probabilities of transitions 
   * corresponding to a certain state.
   * 
   * The input state give the transitions to modify. The float alpha, such that -1< alpha <1 gives the change to apply to the transitions probabilities.
   * 
   * If the input state is v0 v1 ... vt, the starting probability to start with v0, and the transition probabilities to go from v0 to v1, to v1 to v2, etc, from vt-1 to vt are modified as follows. Let p be such starting or transition probability:
   * - if alpha < 0, the probability is decreased to become 
   * p-|alpha|p,
   * - if alpha > 0, the probability is increased to become p+|alpha|(1-p). 
   * 
   * Each transtition vector is then normalized to add up to 1.
   * 
   * @param[in] samp the state for which the transitions are amplified,
   * @param[in] alpha an float between -1 < alpha < 1.
   * 
   * */
  void amplify(std::vector<size_t> samp, Proba_float alpha) {
    //increase starting proba
    auto ps = start_distribution_[samp[0]];
    if(ps + alpha > 0) { start_distribution_[samp[0]] += alpha; }
    else { start_distribution_[samp[0]] /= 1.01; } 
    normalize_proba(start_distribution_.begin(),start_distribution_.end());
    //increase transition proba
    for(size_t t=0; t<num_layers()-1; ++t) {
      auto p = decision_mat_[t][ samp[t] ][ samp[t+1] ];
      if(p + alpha > 0) { 
        decision_mat_[t][ samp[t] ][ samp[t+1] ] += alpha;
      }
      else { decision_mat_[t][ samp[t] ][ samp[t+1] ] /= 1.01; } 
      //ensure the entries of the proba row sum to 1
      normalize_proba(decision_mat_[t][ samp[t] ].begin(), 
                      decision_mat_[t][ samp[t] ].end() );
    }

    // if(alpha > 0) {
    //   //increase starting proba
    //   auto p_start = start_distribution_[samp[0]];
    //   start_distribution_[samp[0]] += alpha*(1-p_start); 
    //   normalize_proba(start_distribution_.begin(),start_distribution_.end());
    //   //increase transition proba
    //   for(size_t t=0; t<num_layers()-1; ++t) {
    //     auto p = decision_mat_[t][ samp[t] ][ samp[t+1] ];//transition proba
    //     p += alpha * (1-p);//increase the transition proba
    //     decision_mat_[t][ samp[t] ][ samp[t+1] ] = p;//
    //     //ensure the entries of the proba row sum to 1
    //     normalize_proba(decision_mat_[t][ samp[t] ].begin(), 
    //                     decision_mat_[t][ samp[t] ].end() );
    //   }
    // }
    // if(alpha < 0 ) {
    //   //decrease starting proba
    //   auto p_start = start_distribution_[samp[0]];
    //   start_distribution_[samp[0]] += alpha*p_start;
    //   normalize_proba(start_distribution_.begin(),start_distribution_.end()); 
    //   //decrease transition proba
    //   for(size_t t=0; t<num_layers()-1; ++t) {
    //     auto p = decision_mat_[t][ samp[t] ][ samp[t+1] ];
    //     p += alpha * p;//decrease the proba
    //     decision_mat_[t][ samp[t] ][ samp[t+1] ] = p;
    //     //ensure the entries of the proba row sum to 1
    //     normalize_proba(decision_mat_[t][ samp[t] ].begin(), 
    //                     decision_mat_[t][ samp[t] ].end() );
    //   }
    // }
  }

//train the decision process
  // template< typename Function >
  // void train(size_t num_iterations, size_t size_batch, Function &reward) {
  //   for(size_t ite = 0; ite < num_iterations; ++ite) {
  //     std::vector< Proba_float > rewards(size_batch);
  //     std::vector< std::vector<size_t> > samples(size_batch);
  //     for(size_t ite_batch=0; ite_batch<size_batch; ++ite_batch) {
  //       samples[ite_batch] = sample();
  //       rewards[ite_batch] = reward(samples[ite_batch]); 
  //     }

  //     //push good behaviors / don't penalize bad ones
  //     auto cp_rewards(reward);
  //     sort(cp_rewards.begin(),cp_rewards.end());


  //   }
  // }

/** \brief For a given state, try, for each layer, to change a node locally to see if it improves the reward.*/
  template< typename Function >
  std::vector<size_t> local_improve(std::vector<size_t> &state, Function &reward) 
  {
    auto max_rew = reward(state);
    std::vector<size_t> imp_state(state);
    //for all layers
    for(size_t i = 0; i<imp_state.size(); ++i) {
      size_t new_s = state[i];
      //try to change state locally to improve reward
      for(size_t j=0; j<size_layer(i); ++j) {
        imp_state[i] = j;
        auto curr_rew = reward(imp_state);
        if(curr_rew > max_rew) {
          max_rew = curr_rew;
          new_s = j;
        }
      }
      imp_state[i] = new_s;
    }
    return imp_state;
  }

  std::string to_string(std::vector<size_t> &state) {
    std::string str = "";
    for(auto i : state) { str = str + std::to_string(i) + " "; }
    return str;
  }

  /** \brief Train the Markov decision process by sampling a state, compute its reward and turn it into an amplifying value alpha, and amplify the probabilities of transtion.
   * 
   * @param[in] num_iterations the number of runs for the training,
   * @param[in] reward a reward function to apply on a state ; it returns values in (-1,1) that can be given to the amplify method.
   **/
  template< typename Function >
  std::vector<size_t> train(size_t num_iterations, Function &reward) {
    
    Proba_float max_rew = -1;
    std::vector<size_t> max_state;

    for(size_t ite = 0; ite < num_iterations; ++ite) {
      Proba_float rew;
      auto samp = sample();
      rew = reward(samp); //between -1 and 1
      if(rew > max_rew) { max_rew = rew; max_state = samp; }
      /////
      std::cout << "#it: " << ite << " -- Reward = " << rew << "  ---  ";
      for(auto x : samp) { std::cout << x << " "; }
      std::cout << "    -- Max = " << max_rew << " [";
      for(auto x : max_state) { std::cout << x << " "; }
      std::cout << "] ";
      std::cout << std::endl;
      

      /////
      // Proba_float factor = 10.*(Proba_float)(ite+1)/(Proba_float)num_iterations;
      amplify(samp,rew);

      for(auto p : start_distribution_) { std::cout << p << " "; }
        std::cout << "\n\n";
    }
    return max_state;
  }

/** \brief Return the number of layers in the Markov decision process.*/
  size_t num_layers() { return decision_mat_.size()+1; }
/** Return the size of a layer given by its index.*/
  size_t size_layer(size_t i) {
    if(i<decision_mat_.size()) {
      if(i==0) { return start_distribution_.size(); }
      return decision_mat_[i-1].num_columns();
    }
    return 0;
  }
private:
  /* Normalize the float numbers in a range such that they sum to 1.
   */
  template<typename FloatIterator >
  void normalize_proba(FloatIterator beg, FloatIterator end) {
    auto sum = Proba_float(0);
    for(auto it=beg; it!=end; ++it) { sum += *it; }
    for(auto it=beg; it!=end; ++it) { *it /= sum; }
  }

private:
  //nodes[l][i] contains the data type of a node in the layer l, l=0...t-1, with label i
  // std::vector< std::vector<Node_type> > nodes_;
  //The Markov decision process is encoded as a multipartite graph 
  //V0 ... Vt with transition |Vi| x |Vi+1| matrices Mi, 
  //i=0...t-1, for 
  //transitioning from the nodes of Vi to the nodes of Vi+1
  /* probability distribution over the nodes of V0 for start*/
  std::vector< Proba_float > start_distribution_;
  //the matrices decision_mat_[i] = Mi
  std::vector< Dense_matrix< R > > decision_mat_;
  //the field a real numbers used to work with multiprecision float for the probabilities
  R R_;
};

} //namespace kumquat

#endif // KUMQUAT_MARKOV_DECISION_H_












