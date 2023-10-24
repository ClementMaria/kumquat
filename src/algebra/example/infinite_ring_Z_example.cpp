/*    This file is part of the KUMQUAT Library - https://kumquat.inria.fr/ 
 *    - released under  GNU GENERAL PUBLIC LICENSE v.3
 *    See file LICENSE or go to https://kumquat.inria.fr/licensing/ for full 
 *    license details.
 *    Author(s):       Clément Maria
 *
 *    Copyright (C) 2023 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include <iostream>
#include <kumquat/Z.h>

using namespace kumquat;

int main() {
  Z_int intZ;

  for(int i=0; i<10; ++i) {
    auto a = intZ.element(i);
    auto inv_a = intZ.multiplicative_inverse(a);
    std::cout << i << " " << a << " " << inv_a << "\n";
  }

  auto unit = intZ.multiplicative_identity();
  auto zero = intZ.additive_identity();
  auto x = intZ.plus(intZ.element(4),intZ.element(5));//4+5 mod 6 == 3
  auto y = intZ.times(intZ.element(3),intZ.element(4));//3*4 mod 6 = 0
  std::cout << "Multiplicative identity = " << unit << "\n";
  std::cout << "Additive identity = " << zero << "\n";
  std::cout << "4+5 = " << x << "\n";
  std::cout << "3*4 = " << y << "\n";
  std::cout << "Order of 4 in Z = " << intZ.order(intZ.element(4)) << "\n";
  std::cout << "4 to the power of 5 = " << intZ.pow(4,5) << "\n";

  for(int i=1; i < 200; ++i) {
    std::cout << "2 to the power of " << i << " = " << intZ.pow(2,i) << "\n";
  }

  return 0;
}