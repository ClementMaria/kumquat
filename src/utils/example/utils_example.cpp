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

#include <kumquat/utils.h>
#include <boost/multiprecision/gmp.hpp>

using namespace kumquat;

// using mp = boost::multiprecision::mpz_int;
// using Qmp = Q<mp>;
// using Q_U1mp = Q_U1<mp>;
// using Q_mod_Zmp = Q_mod_Z<mp>;
// using Zmp = Z<mp>;
// using Z_mod_nZmp = Z_mod_nZ<mp>;

// using Qint = Q<int>;
// using Q_U1int = Q_U1<int>;
// using Q_mod_Zint = Q_mod_Z<int>;
// using Zint = Z<int>;
// using Z_mod_nZint = Z_mod_nZ<int>;


int main(int argc, char * argv[]) {

  // int x = 17;
  // int y = 22;
  // int z = 4;
  // int t = 3;

  // Qmp qmp;
  // std::cout << kumquat::to_string<Qmp>( qmp.element(x,y) ) << "\n";
  // Qint qint;
  // std::cout << kumquat::to_string<Qint>( qint.element(x,y) ) << "\n";

  // Q_U1mp qu1mp;
  // std::cout << kumquat::to_string<Q_U1mp>( qu1mp.element(x,y,z,t) ) << "\n";
  // Q_U1int qu1int;
  // std::cout << kumquat::to_string<Q_U1int>( qu1int.element(x,y,z,t) ) << "\n";

  // Q_mod_Zmp qmodzmp;
  // std::cout << kumquat::to_string<Qmp>( qmodzmp.element(x,y) ) << "\n";
  // Q_mod_Zint qmodzint;
  // std::cout << kumquat::to_string<Q_mod_Zint>( qmodzint.element(x,y) ) << "\n";

  // Zmp zmp;
  // std::cout << kumquat::to_string<Zmp>( zmp.element(x) ) << "\n";
  // Zint zint;
  // std::cout << kumquat::to_string<Zint>( zint.element(x) ) << "\n";

  // Z_mod_nZmp zmodnzmp;
  // std::cout << kumquat::to_string<Z_mod_nZmp>( zmodnzmp.element(x) ) << "\n";
  // Zint zmodnzint;
  // std::cout << kumquat::to_string<Z_mod_nZint>( zmodnzint.element(x) ) << "\n";

  return 0;
}