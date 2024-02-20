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

#ifndef KUMQUAT_PLAT_BRAID_H_ 
#define KUMQUAT_PLAT_BRAID_H_

#include <boost/multiprecision/gmp.hpp>
#include <set>
#include <map>

namespace kumquat {

/**
 * @brief      Data structure for a plat braid.
 */
class Plat_braid {
public:
  //must be signed
  typedef int Crossing_index;

  // Plat_braid() : num_strands_(-1) {};

  // Plat_braid(size_t num_strands, std::vector< std::map<Crossing_index,int> >& braid) : num_strands_(num_strands), braid_(braid) {}
  Plat_braid(size_t num_strands) : num_strands_(num_strands), braid_() {}

  int num_crossings() {
    int num_cr = 0;
    for(auto &layer : braid_) {
      for(auto twist_reg : layer) {
        num_cr += std::abs(twist_reg.second);
      }
    }
    return num_cr;
  }

  int num_strands() { return num_strands_; }

  const std::vector< std::map< Crossing_index, int > >& braid() const { return braid_; }

/**
 * @brief      Adds a twist region on top of the braid.
 *
 * @param[in]  i            Index of the crossing. Must be positive.
 * @param[in]  num_twists  The number twists (positive or negative).
 */
  void add_twist(Crossing_index i, int num_twists) {
    if(i > (int)num_strands_-2) {
      std::cerr << "Invalid crossing index i = " << i <<".\n";
      return;
    }
    if(i < 0) {
      std::cerr << "Invalid crossing index i = " << i <<".\n";
      return;
    }
    if(num_twists == 0) { return; }
    if(braid_.empty()) {
      braid_.emplace_back();//-> create new layer
      braid_.begin()->emplace(i,num_twists);
      return;      
    }
    auto top_layer = braid_.rbegin();
    auto find_i = top_layer->find(i);
    if(find_i != top_layer->end()) {//sigma_i alreayd in top layer
      find_i->second += num_twists;//add the twists
      if(find_i->second == 0) {//if twist region cancelled
        top_layer->erase(find_i);//remove sigma_i
        if(top_layer->empty()) { braid_.pop_back(); }//empty layer?
      }
    }
    else {//sigma_i not in top layer
      //check for absence of conflict
      if(top_layer->find(i-1) != top_layer->end()
         || top_layer->find(i+1) != top_layer->end() ) {//conflict
        braid_.emplace_back();//-> create new layer
        braid_.rbegin()->emplace(i,num_twists);
      }
      else {//no conflict, insert the twist region in top layer
        top_layer->emplace(i,num_twists);
      }
    }
  }


/**
 * \brief Compute the number of compoennts of the standard closure of the braid.*/
  size_t num_components() {
    std::set<int> not_visited;
    for(int i=1; i<(int)num_strands_; ++i) {
      not_visited.insert(i);
    }
    size_t num_comp = 0;
    int curr_idx = 0;
    while(true) {
      //image of curr_idx under the permutation
      curr_idx = bottom_to_top(curr_idx);
      auto find_idx = not_visited.find(curr_idx);
      if(find_idx == not_visited.end()) {//finished reading a component
        ++num_comp;
        
        if(not_visited.empty()) { return num_comp; }
        curr_idx = *(not_visited.begin());//pick a new index
        not_visited.erase(not_visited.begin());
      }
      else {
        not_visited.erase(find_idx);
      }
    }
    return num_comp;
  }

/** Return the oriented Gauss code for the braid.
 * 
 * crossings in a same layer: [a0][a2] ... [a_n/2] are 
 * reorganised as
 *                   [a_n/2]
 *             ...
 *     [a2]
 * [a0]
 * and crossings are labelled (strating at 1) according to there 
 * height, from bottom to top.
 * 
 * Must be a knot.
 */
  std::string oriented_Gauss_code() {
    //map[k,j] == l means that the labels of the consecutive twists
    // on strand j|j+1 in layer k start at label l (included).
    std::map< std::pair<int,int>, int > layer_strand_to_label;

    int curr_label = 1;
    int layer_idx = 0;
    for(auto &layer : braid_) {
      for(auto twist_reg : layer) {
        layer_strand_to_label.emplace(std::make_pair(layer_idx,twist_reg.first),curr_label);
        curr_label += std::abs(twist_reg.second);//+= nb of cross in region
      }
      ++layer_idx;
    }

    // //////
    // for(auto ppp : layer_strand_to_label) {
    //   std::cout << ppp.first.first << "," << ppp.first.second << "|" << ppp.second << "\n";
    // }
    // //////

    std::string ogc = "";
    int curr_idx = 0;
    do {
      //image of curr_idx under the permutation, track the crossings
      curr_idx = bottom_to_top_str(curr_idx, ogc, layer_strand_to_label);
    } while(curr_idx != 0);

    return ogc;
  }

private:
  /** Returns the index of the top of strand of index idx.*/
  int bottom_to_top(int idx) {
    int curr_idx = idx;
    for(auto & layer : braid_) {
      auto find_idx = layer.find(curr_idx-1);
      if(find_idx != layer.end()) {
        if(find_idx->second % 2) {//odd number of twists
          --curr_idx;
        }
      }
      else {
        find_idx = layer.find(curr_idx);
        if(find_idx != layer.end()) {
          if(find_idx->second % 2) {//odd number of twists
            ++curr_idx;
          }
        }
      }
    }
    return curr_idx;
  }

  /** Returns the index of the top of strand of index idx.*/
  int bottom_to_top_str(int idx, std::string &ogc, std::map< std::pair<int,int>, int > l_to_lab) {

    int curr_idx = idx;
    int layer_idx = 0;
    for(auto & layer : braid_) {
      auto find_idx = layer.find(curr_idx-1);
      if(find_idx != layer.end()) {
        //find l such that all crossings of twist reg index from [l; l+T-1]
        int lab_twist = 
                    l_to_lab[std::make_pair(layer_idx,curr_idx-1)];

        std::string sign, direction;
        if(find_idx->second > 0) { sign = "-"; direction = ">"; }
        else { sign = "+"; direction = ">"; }

        for(int i=lab_twist; i<lab_twist+ std::abs(find_idx->second); ++i) 
        {
          ogc += " " + sign + direction + std::to_string(i);

          if(sign == "-") { sign = "+"; }
          else { sign = "-"; }
          if(direction == "<") { direction = ">"; }
          else { direction = "<"; }
        }
        //odd number of twists
        if(find_idx->second % 2) { --curr_idx; }
      }
      else {
        find_idx = layer.find(curr_idx);
        if(find_idx != layer.end()) {
          //find l such that all crossings of twist reg index from [l; l+T-1]
          int lab_twist = 
                      l_to_lab[std::make_pair(layer_idx,curr_idx)];

          std::string sign, direction;

          if(find_idx->second > 0) { sign = "+"; direction = "<"; }
          else { sign = "-"; direction = "<"; }

          for(int i=lab_twist; i<lab_twist+ std::abs(find_idx->second); ++i) 
          {
            ogc += " " + sign + direction + std::to_string(i);

            if(sign == "-") { sign = "+"; }
            else { sign = "-"; }
            if(direction == "<") { direction = ">"; }
            else { direction = "<"; }
          }
          //odd number of twists
          if(find_idx->second % 2) { ++curr_idx; }
        }
      }
      ++layer_idx;
    }
    return curr_idx;
  }


private:
  /**
   * The number of strands in the braid. Strands are indexed from 0 to num_strands-1 included.
   */
  size_t num_strands_;
  /**
   * A plat braid is a braid represented in layers, starting from the bottom towards the top. Each layer consists of twists regions in parallel that can commute in the braid group. In particular, no consecutive braid generators \f$\sigma_i\f$ and \f$\sigma_{i+1}\f$ can be both involved in the same layer.
   * 
   * Each layer is represented by and std::map, containing pairs (i,n) such that i is the index of the generator \f$\sigma_i\f$, and \f$n \neq 0\f$ is the length of the twist region. The convention is the following:
   * - strands are indexed from 0 to num_strands-1 included,
   * - a pair (i,n) satisfies , 0 <= i < num_strands_, and n is positive for |n| consecutive positive crossings \f$\sigma_i\f$, and n is negative for |n| consecutive negative crossings \f$\sigma_i^{-1}\f$
   * 
   * braid_[j] is associated to layer j, starting from j=0.
   */
  std::vector< std::map< Crossing_index, int > > braid_;

};

std::ostream& operator<<(std::ostream& os, const Plat_braid& b) 
{
  for(auto layer_it=b.braid().rbegin(); layer_it!=b.braid().rend(); ++layer_it) {
    for(auto i_tw : *layer_it) {
      os << "(s_" << i_tw.first << ")^" << i_tw.second << " ";
    }
    os << "\n";
  }
  return os;
}

} //namespace kumquat

#endif // KUMQUAT_PLAT_BRAID_H_
