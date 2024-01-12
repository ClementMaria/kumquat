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

#ifndef KUMQUAT_RATIONAL_FUNCTION_H_ 
#define KUMQUAT_RATIONAL_FUNCTION_H_

#include <vector>
#include <numeric>
#include <kumquat/number_theory.h>
#include <boost/multiprecision/gmp.hpp>//boost multiprecision wrap over gmp
// #include <boost/integer/extended_euclidean.hpp>//boost extended gcd
#include <boost/math/tools/polynomial.hpp>

namespace kumquat {

/** \brief A data type for multi-precision rational function represented by two polynomials (P,Q) such that gcd(P,Q) = cste and P and Q have multiprecision integer coefficients.
 */
class Rational_function_integral_mp {
public:
/** \brief An integer type for the coefficients.*/
  typedef boost::multiprecision::mpz_int Coefficient;//type of coefficients
/** The boost polynomial type.*/
  typedef boost::math::tools::polynomial<Coefficient> Polynomial;

/** \brief Default initialization to 0/1.*/
  Rational_function_integral_mp() : numerator_({{0}}), denominator_({{1}}) {}
/** \brief Initialization with an integer z, z/1.*/
  Rational_function_integral_mp(int z) : numerator_({{(Coefficient)z}}), denominator_({{1}}) {}
/** \brief Initialization with an integer z, z/1.*/
  template<typename RangeType>
  Rational_function_integral_mp(RangeType &num, RangeType &den) 
  : numerator_(num.begin(),num.end()), denominator_(den.begin(),den.end()) {
    normalize();
  }

// copy assignment
Rational_function_integral_mp& operator=(const Rational_function_integral_mp& other) 
{
  numerator_ = other.numerator_;
  denominator_ = other.denominator_;  
  return *this;
}
// move assignment
Rational_function_integral_mp& operator=(Rational_function_integral_mp&& other) noexcept
{
  numerator_ = other.numerator_;
  denominator_ = other.denominator_;  
  return *this;
}

Rational_function_integral_mp& operator+=(const Rational_function_integral_mp& rhs) 
{ 
  numerator_ = numerator_*rhs.denominator_ + denominator_*rhs.numerator_;
  denominator_ *= rhs.denominator_;
  normalize();             
  return *this;
}

friend Rational_function_integral_mp operator+(Rational_function_integral_mp lhs,        
                   const Rational_function_integral_mp& rhs) 
{
  lhs += rhs;
  return lhs;
}

Rational_function_integral_mp& operator-=(const Rational_function_integral_mp& rhs) 
{             
  numerator_ = numerator_*rhs.denominator_ - denominator_*rhs.numerator_;
  denominator_ *= rhs.denominator_;
  normalize();             
  return *this;
}

friend Rational_function_integral_mp operator-(Rational_function_integral_mp lhs, const Rational_function_integral_mp& rhs) 
{
  lhs -= rhs;
  return lhs;
}

Rational_function_integral_mp& operator*=(const Rational_function_integral_mp& rhs) 
{                 
  numerator_ *= rhs.numerator_;
  denominator_ *= rhs.denominator_;
  normalize();             
  return *this; 
}

friend Rational_function_integral_mp operator*(Rational_function_integral_mp lhs, const Rational_function_integral_mp& rhs) 
{
  lhs *= rhs;
  return lhs;
}

Rational_function_integral_mp& operator/=(const Rational_function_integral_mp& rhs) 
{                     
  numerator_ *= rhs.denominator_;
  denominator_ *= rhs.numerator_;
  normalize();             
  return *this;
}

friend Rational_function_integral_mp operator/(Rational_function_integral_mp lhs,  const Rational_function_integral_mp& rhs) 
{
  lhs /= rhs;
  return lhs;
}

inline bool operator==(const Rational_function_integral_mp& lhs, const Rational_function_integral_mp& rhs) { 
  return (lhs.numerator_ == rhs.numerator_) && (lhs.denominator_ == rhs.denominator_);
 }
inline bool operator!=(const Rational_function_integral_mp& lhs, const Rational_function_integral_mp& rhs) { return !(lhs == rhs); }


private:
//make sure the numerator and denominator are coprime polynomials, i.e., gcd = integer
  void normalize(Element &rf) {
    auto gcd_poly = boost::math::tools::gcd(numerator_,denominator_);
    numerator_ /= gcd_poly;
    denominator_ /= gcd_poly; 
  }


private:
  Polynomial numerator_;
  Polynomial denominator_;
};

}  //namespace kumquat

#endif //KUMQUAT_RATIONAL_FUNCTION_INTEGRAL_MP_H_