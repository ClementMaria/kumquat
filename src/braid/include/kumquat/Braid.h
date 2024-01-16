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

#ifndef KUMQUAT_BRAID_H_ 
#define KUMQUAT_BRAID_H_

#include <boost/multiprecision/gmp.hpp>

namespace kumquat {

/** Finite cyclic group \bigoplus_{k,n} \left(Z/p^k Z\right)^n, of elements of order a power of a fixed prime p.
 * 
 **/
class Braid {
  //must be signed
  typedef int Crossing_index;

  Braid(size_t num_strands, std::vector<Crossing_index> braid) : num_strands_(num_strands_), braid_(braid) {}

private:
  //a braid is represented by a number of strands, and a list of crossings read from bottom to top. Strands are labeled from 1 to k (included). braid_[i] = j means that the i-th generator of the braid is sigma_|j| (involves strands j and j+1), and if j > 0, strands j goes over strand j+1, and if j < 0, strands j goes under strand j+1. 
  size_t num_strands_;
  std::vector< Crossing_index > braid_;

};

} //namespace kumquat

#endif // KUMQUAT_X_XX_H_
