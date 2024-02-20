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

#ifndef KUMQUAT_Z_ZZ_H_ 
#define KUMQUAT_Z_ZZ_H_

#include <vector>
#include <numeric>
#include <kumquat/number_theory.h>
#include <boost/multiprecision/gmp.hpp>//boost multiprecision wrap over gmp
// #include <boost/integer/extended_euclidean.hpp>//boost extended gcd
#include <boost/math/tools/polynomial.hpp>

namespace kumquat {

/** \brief A data type for a scalar element of type...
 * 
 * The scalar are represented by...
 * 
 * \implements ScalarType
 */
class Z_zz {
public:

/** \brief Default initialization to 0/1.*/
  Z_zz()  {}
/** \brief Initialization with an integer z, z/1.*/
  Z_zz(int z) {}

// Copy assignment
Z_zz& operator=(const Z_zz& other) 
{
  return *this;
}
// Move assignment
Z_zz& operator=(Z_zz&& other) noexcept
{
  numerator_ = other.numerator_;
  denominator_ = other.denominator_;  
  return *this;
}
//operator+=
Z_zz& operator+=(const Z_zz& rhs) 
{ 
  ...           
  return *this;
}
//operator-=
Z_zz& operator-=(const Z_zz& rhs) 
{             
  ...
  return *this;
}
//operator*=
Z_zz& operator*=(const Z_zz& rhs) 
{                 
  ...
  return *this; 
}
//operator/= if applies
Z_zz& operator/=(const Z_zz& rhs) 
{                     
  ...
  return *this;
}
//operator==
inline bool operator==(const Z_zz& lhs, const Z_zz& rhs) { 
  return ...;
}

Z_zz& operator%=(const Z_zz& rhs) {
  ...
  return *this;
}
//other %= etc

//operator+ deduced from operator+=
friend Z_zz operator+(Z_zz lhs, const Z_zz& rhs) 
{
  lhs += rhs;
  return lhs;
}
//operator- deduced from operator-=
friend Z_zz operator-(Z_zz lhs, const Z_zz& rhs) 
{
  lhs -= rhs;
  return lhs;
}
//operator* deduced from operator*=
friend Z_zz operator*(Z_zz lhs, const Z_zz& rhs) 
{
  lhs *= rhs;
  return lhs;
}
//operator/ deduced from operator/=
friend Z_zz operator/(Z_zz lhs,  const Z_zz& rhs) 
{
  lhs /= rhs;
  return lhs;
}
//operator/ deduced from operator/=
friend Z_zz operator%(Z_zz lhs,  const Z_zz& rhs) 
{
  lhs %= rhs;
  return lhs;
}
//operator!= defined as the negation of operator==
inline bool operator!=(const Z_zz& lhs, const Z_zz& rhs) { return !(lhs == rhs); }

};

}  //namespace kumquat

#endif //KUMQUAT_RATIONAL_FUNCTION_INTEGRAL_MP_H_