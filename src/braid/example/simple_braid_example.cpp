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
#include <kumquat/Plat_braid.h>

using namespace kumquat;

int main() {//(int argc, char * argv[]) {

  std::vector< std::pair<int,int> > cross;
  cross.emplace_back(1,2);   cross.emplace_back(2,4);
  Braid b(3,cross);
  b.add_twist(2,-5);
  b.add_twist(0,-1);

  std::cout << b << "\n\n";


  Plat_braid pb(4);//plat braid on 4 strands
  pb.add_twist(0,-2);
  pb.add_twist(1,1);//new layer
  pb.add_twist(1,-1);//cancels the s_1 twist
  pb.add_twist(0,-2);//more twisting on 0
  pb.add_twist(2,3);//same bottom layer
  pb.add_twist(1,4);//new layer
  pb.add_twist(0,-2);//new layer
  pb.add_twist(2,3);//same top layer
  pb.add_twist(0,1);//same top layer
  pb.add_twist(1,1);//same top layer
  pb.add_twist(2,1);//same top layer

  //is the braid:
  //
  //        [ 1 ] 
  //     [  1]
  //  [  1]
  //   | |  | |
  //  [ -2][ 3 ]
  //  | [ 4 ] |
  // [ -4] [  3]
  //  | |   | |


  std::cout << pb << "\n";
  std::cout << "#crossings = " << pb.num_crossings() << "\n";
  std::cout << "#strands = " << pb.num_strands() << "\n";
  std::cout << "#components = " << pb.num_components() << "\n";

  std::cout << pb.oriented_Gauss_code() << "\n";

  return 0;
}