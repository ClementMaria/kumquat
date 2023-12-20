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

/** A decision matrix for a Markov decision process. The process is a complete 
 * multipartite graph V0 V1 ... Vt, represented by a set of transition matrices 
 * M0 ... M_t-1. 
 * 
 * Matrix Mi is a (Vi x Vi+1) dense matrix, representing the transition of local state at time i to time i+1. Precisely, Mi[p] is a row vector of length Vi+1 containing the probability distribution of transitioning to the nodes of Vi+1, knowing that we are in the p-th node of Vi. 
 * 
 * Mi[p][q] is the probability to jumping to q \in Vi+1, knowing that we are in p at time i.
 * 
 **/
class Markov_decision {

  typedef R::Element Proba_float;

  //balanced probabilities 
  //num_states[i] contains the size of |Vi|, i=0...t
  Markov_decision(std::vector<size_t> num_states) 
  : num_layers_(num_states.size()) 
  {
    //idd starting point
    start_distribution_(num_states[0], R_.element(1)/R_.element(num_states[0]));    
    decision_mat_( num_layers_-1 ); 
  
    for(size_t t=0; t<num_layers_-1; ++t) {
      size_t num_rows = num_states[t];
      size_t num_cols = num_states[t+1];
      decision_mat_[t] = Dense_matrix<R>(num_rows,num_cols,R_); 
      Proba_t f = [&](size_t i, size_t j) { return 1/num_cols; };    
      decision_mat_.[t].fill(f);
    }
  }

//return a sample using the decision process
  std::vector< size_t > sample() {
    std::random_device rd;
    std::mt19937 gen(rd());

    std::vector<size_t> samp(num_layers_);

    std::discrete_distribution<> dist_start(start_distribution_.begin(),start_distribution_.end());
    samp[0] = dist_start(gen);

    size_t curr_node = samp[0];
    for(size_t t=0; t<num_layers_-1; ++t) {
      std::discrete_distribution<> dist_t_to_next(
                        decision_mat_[t][curr_node].begin(), 
                        decision_mat_[t][curr_node].end());
      curr_node = dist_t_to_next(gen);
      samp[t+1] = curr_node;
    }
    return samp;
  }

// if alpha > 0, a transition probability p is amplified to become p + |alpha|(1-p). If alpha < 0, p becomes p - |alpha| p
// -1 < alpha < 1 
  void amplify(std::vector<size_t> samp, Proba_float alpha) {
    if(alpha > 0) {
      for(size_t t=0; t<num_layers_-1; ++t) {
        auto p = decision_mat_[t][ samp[t] ][ samp[t+1] ];//transition proba
        p += alpha * (1-p);//increase the transition proba
        decision_mat_[t][ samp[t] ][ samp[t+1] ] = p;//
        //ensure the entries of the proba row sum to 1
        normalize_proba(decision_mat_[t][ samp[t] ].begin(), 
                        decision_mat_[t][ samp[t] ].end() );
      }
    }
    if(alpha < 0 ) {
      for(size_t t=0; t<num_layers_-1; ++t) {
        auto p = decision_mat_[t][ samp[t] ][ samp[t+1] ];
        p += alpha * p;//decrease the proba
        decision_mat_[t][ samp[t] ][ samp[t+1] ] = p;
        //ensure the entries of the proba row sum to 1
        normalize_proba(decision_mat_[t][ samp[t] ].begin(), 
                        decision_mat_[t][ samp[t] ].end() );
      }
    }
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

  template< typename Function >
  void train(size_t num_iterations, Function &reward) {
    for(size_t ite = 0; ite < num_iterations; ++ite) {
      Proba_float rew;
      auto samp = sample();
      rew = reward(samp); //between -1 and 1
      /////
      std::cout << "Reward = " << rew << "  ---  ";
      for(auto x : samp) { std::cout << x << " "; }
      std::cout << std::endl;
      /////
      amplify(samp,rew);
    }
  }

private:
  template<typename FloatIterator >
  void normalize_proba(FloatIterator beg, FloatIterator end) {
    auto sum = Proba_float(0);
    for(auto it=beg; it!=end; ++it) { sum += *it; }
    for(auto it=beg; it!=end; ++it) { *it /= sum; }
  }

private:
  //M[i][j] contains a float that gives the probability of going to state M[i+1][j] when we are in state M[i][j]. The overall decision is the set of col idx visited over time (time = row index)
  std::vector< Proba_float > start_distribution_;
  std::vector< Dense_matrix< Proba_float > > decision_mat_;
  size_t num_steps_;

  R R_;
};

} //namespace kumquat

#endif // KUMQUAT_MARKOV_DECISION_H_
