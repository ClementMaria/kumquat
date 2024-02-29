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

#ifndef KUMQUAT_MARKOV_DECISION 
#define KUMQUAT_MARKOV_DECISION

#include <kumquat/R.h>
#include <kumquat/Dense_matrix.h>

namespace kumquat {

/** A decision matrix for a Markov decision process. The process is 
 * a complete 
 * multipartite graph V0 V1 ... Vt, represented by a set of transition matrices 
 * M0 ... M_t-1. 
 * 
 * Matrix Mi, i=0...t-1. is a (Vi x Vi+1) dense matrix, representing the 
 * transition of local state at time i to time i+1. Precisely, 
 * Mi[s][...] is a row vector of length |Vi+1| containing the probability 
 * distribution of transitioning to the nodes of Vi+1, knowing that 
 * we are in the s-th node of Vi. 
 * 
 * In other terms, Mi[p][q] is the probability to jumping to q \in Vi+1, knowing 
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
  //default constructor with empty states and matrices
  Markov_decision() {}

/** \brief Copy constructor.*/
  Markov_decision(const Markov_decision& other) {
    start_distribution_ = other.start_distribution_;
    decision_mat_ = other.decision_mat_;
    R_ = other.R_;
  }
/** \brief Move constructor.*/
  Markov_decision(Markov_decision&& other) noexcept {
    start_distribution_ = std::move(other.start_distribution_);
    decision_mat_ = std::move(other.decision_mat_);
    R_ = std::move(other.R_);
  }
/** brief Destructor.*/
  ~Markov_decision() {};
/**
 * @brief Copy assignment.
 * @details Copy assignment on all members.
 * 
 * @param other Original to copy.
 * @return *this.
 */
  Markov_decision& operator=(const Markov_decision& other) {
    if(this != &other) { 
      start_distribution_ = other.start_distribution_;
      decision_mat_ = other.decision_mat_;
      R_ = other.R_;
    }
    return *this;
  }
/**
 * @brief Move assignment.
 * @details Move assignment on all members.
 * 
 * @param other Original to move.
 * @return *this.
 */
  Markov_decision& operator=(Markov_decision&& other) noexcept
  {
    if(this != &other) { 
      start_distribution_ = std::move(other.start_distribution_);
      decision_mat_ = std::move(other.decision_mat_);
      R_ = std::move(other.R_);
    }
    return *this;
  }


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
    decision_mat_.reserve( num_layers-1 ); 
  
    for(size_t t=0; t<num_layers-1; ++t) {
      size_t num_rows = num_states[t];
      size_t num_cols = num_states[t+1];
      decision_mat_.emplace_back(num_rows,num_cols,R_); 
      //fill the matrices with appropriately normalize uniform proba
      auto f = [&]([[maybe_unused]] size_t i, [[maybe_unused]] size_t j)->Proba_float { return (Proba_float)(1)/(Proba_float)num_cols; };   
      decision_mat_[t].fill(f);
    }
  }

/**
 * @brief Reinitializes the transition matrices such that the probability of transition from state[i] to state[i+1] is exactly alpha, and all other transition are uniform. 
 * @details Useful to explore the surroundings of a state.
 * 
 * @param state The state to favor.
 * @param alpha The fraction of probability.
 */
  void initialize(std::vector<size_t> state = std::vector<size_t>(), Proba_float alpha = 0.5)
  {
    if(state.empty()) 
    {
      //idd distribution for starting point
      Proba_float unif = R_.element(1)/R_.element(size_layer(0));
      for(size_t i=0; i<size_layer(0); ++i) {
        start_distribution_[i] = unif;
      }
      //for each transition matrix    
      for(size_t t=0; t<num_layers()-1; ++t) {
        size_t num_rows = size_layer(t);
        size_t num_cols = size_layer(t+1);
        //fill the matrices with appropriately normalize uniform proba
        auto f = [&]([[maybe_unused]] size_t i, 
                     [[maybe_unused]] size_t j)->Proba_float { 
          return (Proba_float)(1)/(Proba_float)num_cols; 
        };
        decision_mat_[t].fill(f);
      }
    }
    else {//state != empty
      //idd distribution for starting point
      

      // std::cout << "+++ size_layer(0) = " << size_layer(0) << "\n";
      Proba_float unif_alpha = (Proba_float(1) - alpha)/Proba_float(size_layer(0));
      for(size_t i=0; i<size_layer(0); ++i) {
        if(i == state[0]) { start_distribution_[i] = alpha; }
        else { start_distribution_[i] = unif_alpha; }
      }
      //for each transition matrix    
      for(size_t t=0; t<num_layers()-1; ++t) {
        
        // std::cout << "+++ size_layer(" << t << ") = " << size_layer(t) << "\n";
        // std::cout << "----- size_layer(t+1) = " << size_layer(t+1) << "\n";

        size_t num_rows = size_layer(t);
        size_t num_cols = size_layer(t+1);
        //fill the matrices with appropriately normalize uniform proba
        auto f = [&]([[maybe_unused]] size_t i, 
                     [[maybe_unused]] size_t j)->Proba_float { 
          if(i == state[t]) {
            if(j == state[t+1]) { return alpha; }
            else { 
              return (Proba_float(1) - alpha)/Proba_float(num_cols);
            }
          }
          else {//not in row state[t]
            return (Proba_float)(1)/(Proba_float)num_cols;  
          }
        };
        decision_mat_[t].fill(f);
      }
    }
  }
  

 
  /** \brief Return a state sampled from the distribution of the Markov process. 
   * 
   * Transverse the multipartite graph V0 ... Vt from left to right, 
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
      // std::cout << "     - Layer " << i << "/" << imp_state.size()-1 << "\n";
      //prepare batches of size size_batch
      int size_batch = 20;
      int num_batches = size_layer(i)/size_batch +1;
      //local_best[b] == (new_s, rew) means that when trying to locally improve within layer i, the batch b (improve locally for states of labels in [b*size_batch, (b+1)*size_batch-1]), the best valiue attained is rew, when locally setting to state new_s
      std::vector< std::pair<size_t, Proba_float> > local_best(num_batches, std::make_pair(state[i], max_rew));

      // for(int batch=0; batch < num_batches; ++batch) {
      tbb::parallel_for(size_t(0), size_t(num_batches), 
        [&](size_t batch) {
          
          std::vector<size_t> tmp_imp_state(imp_state);
          //try to change state locally to improve reward
          for(size_t j = batch * size_batch ; 
                     j < std::min(size_layer(i), (batch+1)*size_batch); ++j) 
          {
            // std::cout << "          - state " << j << "/" << size_layer(i) << "\n";

            tmp_imp_state[i] = j;
            auto curr_rew = reward(tmp_imp_state);
            if(curr_rew > local_best[batch].second) {
              local_best[batch].second = curr_rew;
              local_best[batch].first = j;
            }
          }
        }
      );

      //now extract the best option among all parallel computations
      auto max_rew_batch = max_rew;
      auto max_state_batch = state[i];
      for(auto pp : local_best) {
        if(pp.second > max_rew_batch) {
          max_rew_batch = pp.second;
          max_state_batch = pp.first;
        }
      }

      max_rew = max_rew_batch;
      imp_state[i] = max_state_batch;

    }

  return imp_state;
    // auto max_rew = reward(state);
    // std::vector<size_t> imp_state(state);
    // //for all layers
    // for(size_t i = 0; i<imp_state.size(); ++i) {
    //   std::cout << "     - Layer " << i << "/" << imp_state.size()-1 << "\n";
    //   size_t new_s = state[i];
    //   //try to change state locally to improve reward
    //   for(size_t j=0; j<size_layer(i); ++j) {
    //     std::cout << "              - state " << j << "/" << size_layer(i) << "\n";

    //     imp_state[i] = j;
    //     auto curr_rew = reward(imp_state);
    //     if(curr_rew > max_rew) {
    //       max_rew = curr_rew;
    //       new_s = j;
    //     }
    //   }
    //   imp_state[i] = new_s;
    // }
    // return imp_state;
  }

  std::string to_string(std::vector<size_t> &state) {
    std::string str = "";
    for(auto i : state) { str = str + std::to_string(i) + " "; }
    return str;
  }

  /** \brief Train the Markov decision process by sampling a batch of states in parallel, compute their rewards, and amplify the probabilities of transitions according to the rewards.
   * 
   * @param[in] num_iterations the number of runs for the training,
   * @param[in] reward a reward function to apply on a state ; it returns any value, that can be given to the amplify method.
   * @param[in] batch_size number of samples taken in parallel, for which rewards are computed. The amplifying is then done on the samples, one after the other. batch_size should divide num_iterations
   **/
  template< typename Function >
  std::vector<size_t> train(size_t num_iterations, Function &reward, size_t batch_size) {

    Proba_float max_rew = -1;
    std::vector<size_t> max_state;

    // std::vector<size_t> prev_samp(0,num_layers());
    // int num_repetitions = 0;
    for(size_t ite = 0; ite < num_iterations; ite += batch_size) {
      
      std::cout << "   " << ite << "/" << num_iterations << "\n";
      // if(ite % (batch_size) == 0) {  std::cout << "-----------iteration #" << ite << "\n"; }
      //after the parallel_for, samples_rew[i] contains a pair made of a sample and its reward value, as computed by thread #i
      std::vector< std::pair< std::vector<size_t>, Proba_float > > samples_rew(batch_size);
      tbb::parallel_for(size_t(0), size_t(batch_size),
        [&](size_t i) {
          std::vector<size_t> samp = sample();
          auto rew = reward(samp); //between -1 and 1
          samples_rew[i] = std::make_pair(samp, rew);
        });
      // if(prev_samp == samp) {
      //   ++num_repetitions;
      //   if(num_repetitions > 9) { return max_state; }
      // }
      // else {
      //   prev_samp = samp;
      //   num_repetitions = 1;
      // }
      // rew = reward(samp); //between -1 and 1
      size_t thread_idx=0;
      for(auto s_r : samples_rew) {
        if(s_r.second > max_rew) {
          max_state = s_r.first;
          max_rew = s_r.second;
          //////////////////////////////////new max
          // std::cout << "train ite " << ite+thread_idx << "  max_rew [" << max_rew << "]   max_state: ";
          // for(auto i : max_state) { std::cout << i << " "; }
          // std::cout << std::endl;
          //////////////////////////////////

        }
        ++thread_idx;
        if(s_r.second != 0) { amplify(s_r.first,s_r.second); }
      }


      /////
      // if(rew > 0) {
      //   std::cout << "#it: " << ite << " -- Reward = " << rew << "  ---  ";
      //   for(auto x : samp) { std::cout << x << " "; }
      //   std::cout << "    -- Max = " << max_rew << " [";
      //   for(auto x : max_state) { std::cout << x << " "; }
      //   std::cout << "] ";
      //   std::cout << std::endl;
      // }

      /////
      // Proba_float factor = 10.*(Proba_float)(ite+1)/(Proba_float)num_iterations;
      // if(rew != 0) { amplify(samp,rew); }

      //max value of start distrib
      // Proba_float max_p_start = 0;
      // for(auto p : start_distribution_) { 
      //   if(p > max_p_start) { max_p_start = p; }
      //   // std::cout << p << " "; 
      // }
      // std::cout << "max start proba = " << max_p_start << "   vs uniform = " << (Proba_float)(1)/(Proba_float)(size_layer(0)) << "\n";
        // std::cout << "\n\n";
    }
    std::cout << "\n";
    return max_state;
  }

/** \brief Return the number of layers in the Markov decision process.*/
  size_t num_layers() { return decision_mat_.size()+1; }
/** Return the size of a layer given by its index.*/
  size_t size_layer(size_t i) {
    // std::cout << "decision_mat_.size() = " << decision_mat_.size() << "   i = " << i << "\n";
    if(i<=decision_mat_.size()) {
      if(i==0) { return start_distribution_.size(); }
      return decision_mat_[i-1].num_columns();
    }
    // std::cout << "X\n";
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












