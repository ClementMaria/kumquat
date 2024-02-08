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

using namespace kumquat;

typedef Markov_decision::Proba_float P_float;

int main() {

  int max_twists = 10;//arbitrary
  int num_strands = 6;//must be even
  int num_layers = 11;//must be odd
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


  to do


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
















