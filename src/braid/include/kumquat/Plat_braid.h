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

  int num_strands() { return num_strands_; }

  std::vector< std::map< Crossing_index, int > >& braid() { return braid_; }
private:
  /**
   * The number of strands in the braid. Strands are indexed from 1 to num_strands included.
   */
  size_t num_strands_;
  /**
   * A plat braid is a braid represented in layers, starting from the bottom towards the top. Each layer consists of twists regions in parallel that can commute in the braid group. In particular, no consecutive braid generators \f$\sigma_i\f$ and \f$\sigma_{i+1}\f$ can be both involved in the same layer.
   * 
   * Each layer is represented by and std::map, containg pairs (i,n) such that i is the index of the generator \f$\sigma_i\f$, and \f$n\f$ is the length of the twist region. The convention is the following:
   * - strands are indexed from 1 to |num_strands| included,
   * - a pair (i,n) satisfied n > 0, and if i > 0, in represents \f$(\sigma_{|i|})^n\f$, and if \f$i<0\f$, the pair represents \f$(\sigma_{|i|})^{-n}\f$.
   * 
   * braid_[i] is associated to layer i, starting from i=0.
   */
  std::vector< std::map< Crossing_index, int > > braid_;

};

} //namespace kumquat

#endif // KUMQUAT_PLAT_BRAID_H_
