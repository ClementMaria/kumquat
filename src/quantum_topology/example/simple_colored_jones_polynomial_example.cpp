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

#include <iostream>
#include <kumquat/Quantum_group.h>

using namespace kumquat;

int main() {
  int N = 2;
  Quantum_group_Uqsl2_gen_q<N> uq;
  
  int num_strands = 4;
  int max_twists = 3;
  std::vector< std::pair<int,int> > braid; braid.reserve(5);
  braid.push_back(std::make_pair(2,2));
  braid.push_back(std::make_pair(-1,1));
  braid.push_back(std::make_pair(-3,1));
  braid.push_back(std::make_pair(2,1));
  Braid
  
  return 0;
}