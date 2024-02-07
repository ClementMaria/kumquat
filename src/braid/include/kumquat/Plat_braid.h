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

namespace kumquat {

/**
 * @brief      Data structure for a plat braid.
 */
class Plat_braid {
public:
  //must be signed
  typedef int Crossing_index;

  Plat_braid(size_t num_strands, std::vector< std::map<Crossing_index,int> >& braid) : num_strands_(num_strands), braid_(braid) {}
  Plat_braid(size_t num_strands) : num_strands_(num_strands), braid_() {}

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
      std::cerr << "Invalid crossing index.\n";
      return;
    }
    if(i < 0) {
      std::cerr << "Invalid crossing index.\n";
      return;
    }
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
