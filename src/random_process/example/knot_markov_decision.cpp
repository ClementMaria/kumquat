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

#include <kumquat/Markov_decision.h>
#include <kumquat/Braid.h>

using namespace kumquat;

typedef Markov_decision::Proba_float P_float;

int main() {

  //A decision process where states v0 ... vt represent braids with the following encoding:
  //we work with braids on k strands, k even, where strands are labeled from 1 to k included.
  //vi = (i,w) where 1 <= |i| < k encode the location of the crossing \sigma_i. If i>0 it is a positive crossing (strand i above strand i+1), and i negative is a negative crossing
  //and w is the length of the twist region.
  int num_strands = 6;//|i| = 1 ... num_strands-1
  int max_twist = 5;//w = 1 ... max_twist

  //each layer has (k-1)*(max_twist)*2 nodes (all sigma_i times number of possible twists times the two possible signs). 
  /* A node (i,w) with 1<= i < k-1, and 1 <= w <= max_twist is 
   * given label: 
   * 
   *  (i-1) * (max_twist) + (w-1) \in [0; (k-1)*max_w - 1]
   * if i>0, and label  
   * (k-1)*max_w + (|i|-1) * (max_twist) + (w-1) \in [(k-1)*max_w; 2*(k-1)*max_w -1] if i<0.  
   */

  //from a label as above, return the pair (i,w)
  auto label_to_cross = [&](size_t &s)->std::pair<int,int> {
    int x = (int)s;
    if(x < (num_strands-1)*max_twist) {//positive crossing
      //x == (i-1) * (max_twist) + (w-1)
      int i = x/max_twist; ++i; //label of crossing
      int w = x - (i-1)*max_twist; ++w; //number of twists
      return std::make_pair(i,w);
    }
    else {//negative crossing
      //x==(k-1)*max_w + (|i|-1) * (max_twist) + (w-1)
      x -= (num_strands-1)*max_twist; //now x==(|i|-1) * (max_twist) + (w-1)
      int i = x/max_twist; ++i; //label of crossing
      int w = x - (i-1)*max_twist; ++w; //number of twists
      return std::make_pair((-1)*i,w);
    } 
  }

  //from a state return a braid
  auto state_to_braid = [&](std::vector<size_t> &samp)->std::vector< std::pair<int,int> > {
    std::vector< std::pair<int,int> > b;
    for(auto s : samp) {
      b.push_back(label_to_cross(s));
    }
    return b;
  };

  //high reward for jones polynomials close to the trivial Jones,
  //and high reward for knots (and not braids)
  auto rew_jones [&](std::vector<size_t> &samp)->P_float
  {


    return ...;
  }

}
















