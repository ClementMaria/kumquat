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
#include <kumquat/Jones_polynomial.h>

using namespace kumquat;

int main() {

  int num_strands = 2;
  int max_twist = 10;
  std::cout << "Initialization of Uq(sl2(C) Jones with num_strands = " << num_strands << " and max_twist = " << max_twist << "\n";
  Jones_polynomial J(num_strands,max_twist);
  

  J.display();



  Plat_braid pbu(num_strands);
  pbu.add_twist(0,1);
  // pbu.add_twist(1,1);
  // pbu.add_twist(2,1);

  auto q_pbu = J.quantum_invariant(pbu);

  std::cout << "Quantum invariant of the unknot: " << J.quantum_invariant_unknot().to_string() << "\n";

  std::cout << "Quantum invariant of unknot with non-standard diagram = " << q_pbu << "\n";

  return 0;
}