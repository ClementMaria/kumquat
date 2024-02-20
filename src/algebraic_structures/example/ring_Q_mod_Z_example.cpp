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
#include <kumquat/Q_mod_Z.h>

using namespace kumquat;

int main() {
  Q_mod_Z_mp qmz;

  Q_mod_Z_mp::Element x = qmz.element(150,65);
  Q_mod_Z_mp::Element y = qmz.element(12,22);

  auto zero = qmz.additive_identity();
  std::cout << "Additive identity = " << zero.first << "/" << zero.second << "\n";
  auto x_plus_y = qmz.plus(x,y);
  std::cout << "150/65 + 12/22 = " << x_plus_y << "\n";
  std::cout << "150/65 = " << x << "\n";
  std::cout << "12/22 = " << y << "\n";
  std::cout << "Order of 150/65 mod Z = " << qmz.order(x) << "\n";

  return 0;
}