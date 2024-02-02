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

/**
 * @brief      This class representing a braid on a fixed number of strands.
 */
class Braid {
public:
/**
 * Type for the indices of the strands/crossings. Must be a signed integer ; strands are labelled from 1 to num_strands, and crossings are labelled from 1 to num_strands-1.
 * 
 * Must be sign to represent \f$\sigma_i\f$ by index i>0 and \f$\sigma_i^{-1}\f$ by index -i < 0.
 */
  typedef int Crossing_index;
/**
 * @brief      Creates a braid with input sequence of twist regions.
 *
 * @param[in]  num_strands  The number strands in the braid.
 * @param      braid        Creates the braid.
 */
  Braid(size_t num_strands, std::vector< std::pair<Crossing_index,int > >& braid) : num_strands_(num_strands), braid_(braid) {}
/**
 * @brief      Return the number of strands in the braid.
 *
 * @return     Number of strands.
 */
  int num_strands() { return num_strands_; }
/**
 * @brief      Return the seuqnece of crossings.
 *
 * @return     Vector of pairs indicating crossing index and number of twists.
 */
  std::vector< std::pair<Crossing_index, int> >& braid() { return braid_; }
/**
 * @brief      Adds a single crossing.
 *
 * @param[in]  i     Index of the crossing. Positive for positive crossing at index i, negative for negative crossing at index |i|.
 */
  void add_crossing(Crossing_index i) {
    add_twist(i,1);
  }
/**
 * @brief      Adds a twist region on top of the braid.
 *
 * @param[in]  i            Index of the crossing. Positive for positive crossing at index i, negative for negative crossing at index |i|.
 * @param[in]  num_twists  The number twists.
 */
  void add_twist(Crossing_index i, int num_twists) {
    if((int)(std::abs(i)) > (int)num_strands_-1) {
      std::cerr << "Invalid crossing index.\n";
      return;
    }
    if(num_twists <= 0) {
      std::cerr << "Invalid number of twists.\n";
      return;
    }
    auto last = braid_.rbegin();
    if(last->first == i) { last->second += num_twists; return; }
    if(last->first == (-1)*i) { 
      last->second -= num_twists; 
      if(last->second == 0) { braid_.pop_back(); return; }
      if(last->second < 0) {
        last->first *= -1; last->second *= -1;
      }
      return;
    }
    braid_.emplace_back(i,1);
  }

private:
  /**
   * The number of strands in the braid,
   */
  size_t num_strands_;
  /**
   * A braid is represented by a number of strands, and a list of crossings read from bottom to top. Strands are labeled from 1 to k (included). braid_[i] = j means that the i-th generator of the braid is sigma_|j| (involves strands j and j+1), and if j > 0, strands j goes over strand j+1, and if j < 0, strands j goes under strand j+1.
   * 
   * This is a compact representation of the braid, where braid_[i] = (j,n) significate that we have n consecutive sigma_|j|, with appropriate sign.
   */
  std::vector< std::pair<Crossing_index, int> > braid_;

};

std::ostream& operator<<(std::ostream& os, const Braid& b) const
{
  for(auto ct : b.braid()) {
    os << ct.first << "^" << ct.second << " ";
  }
  return os;
}

} //namespace kumquat

#endif // KUMQUAT_X_XX_H_
