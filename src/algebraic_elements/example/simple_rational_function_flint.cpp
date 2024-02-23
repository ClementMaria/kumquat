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
#include <fstream>      // std::ofstream
#include <ctime>
#include <string>
#include <utility>  // for std::pair
#include <numeric> //std::gcd
#include <chrono>

#include <kumquat/Rational_function_flint.h> 

  using Rational_function_flint Rf;
  using boost::multiprecision::mpz_int Intmp;

int main (int argc, char* const argv[])
{ 
  Rf a; //default
  Rf b(5); //5 X^0 / 1
  Rf c(std::make_pair(2,14));//14 X^2 / 1

  std::vector< std::pair<int, Intmp> > num;
  num.emplace_back(1,2);
  num.emplace_back(1,-5);
  num.emplace_back(3,-9); //-3X - 9X^3 = -3X(1 + 3X^2)

  std::vector< std::pair<int, Intmp> > den;
  num.emplace_back(0,6);
  num.emplace_back(1,6);
  num.emplace_back(3,18);
  num.emplace_back(2,18);//6+6X+18X^2+18X^3 = 2*3*(1+3X^2)(1+X)

  Rf d(num,den);
  Rf e(d);
  Rf f = e;

  
  return 0;
}


