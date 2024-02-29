/*    This file is part of the KUMQUAT Library - https://kumquat.inria.fr/ 
 *    - released under  GNU GENERAL PUBLIC LICENSE v.3
 *    See file LICENSE or go to https://kumquat.inria.fr/licensing/ for full 
 *    license details.
 *    Author(s):       Clément Maria
 *
 *    Copyright (C) 2024 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include <iostream>
#include <kumquat/Z_mod_nZ.h>

using namespace kumquat;

typedef ABELIAN_GROUP AGroup;

int main() {
  AGroup G(...);

  for(int i=0; i<10; ++i) {
    auto a = Z6.element(i);
    auto inv_a = Z6.multiplicative_inverse(a);
    std::cout << i << " " << a << " " << inv_a << "\n";
  }

  auto unit = Z6.multiplicative_identity();
  auto zero = Z6.additive_identity();
  auto x = Z6.plus(Z6.element(4),Z6.element(5));//4+5 mod 6 == 3
  auto y = Z6.times(Z6.element(3),Z6.element(4));//3*4 mod 6 = 0
  std::cout << "Multiplicative identity = " << unit << "\n";
  std::cout << "Additive identity = " << zero << "\n";
  std::cout << "4+5 = " << x << "\n";
  std::cout << "3*4 = " << y << "\n";
  std::cout << "Order of 4 in Z/6Z = " << Z6.order(Z6.element(4)) << "\n";
  std::cout << "4 to the power of 5 = " << Z6.pow(4,5) << "\n";

  return 0;
}