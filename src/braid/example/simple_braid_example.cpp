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

#include <kumquat/Braid.h>

using namespace kumquat;

int main() {//(int argc, char * argv[]) {

  std::vector< std::pair<int,int> > cross;
  cross.emplace_back(1,2);   cross.emplace_back(-2,4);
  Braid b(3,cross);
  cross.add_twist(2,5);

  std::cout << b << "\n";

  return 0;
}