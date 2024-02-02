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
 * 
 * Is model of concept ScalarFieldOperations
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
  template<typename IntegerType>
  Rational_function_integral_mp(IntegerType z) : numerator_({{(Coefficient)z}}), denominator_({{1}}) {}
/** \brief Initialization with a range.*/
  template<typename RangeType>
  Rational_function_integral_mp(RangeType &num, RangeType &den) 
  : numerator_(num.begin(),num.end()), denominator_(den.begin(),den.end()) {
    normalize();
  }
/** \name Model of ScalarSetOperations
 * 
 * @{ */
/** \brief Copy constructor.*/
  Rational_function_integral_mp(const Rational_function_integral_mp& other) {
    numerator_ = other.numerator_;
    denominator_ = other.denominator_;
  }
/** \brief Move constructor.*/
  Rational_function_integral_mp(Rational_function_integral_mp&& other) noexcept {
    numerator_ = std::move(other.numerator_);
    denominator_ = std::move(other.denominator_);
  }
/** brief Destructor.*/  
  ~Rational_function_integral_mp() {}
/** \brief Copy assignment.*/
  Rational_function_integral_mp& operator=(const Rational_function_integral_mp& other) 
  {
    numerator_ = other.numerator_;
    denominator_ = other.denominator_;  
    return *this;
  }
/** \brief Move assignment.*/
  Rational_function_integral_mp& operator=(Rational_function_integral_mp&& other) noexcept
  {
    numerator_ = other.numerator_;
    denominator_ = other.denominator_;  
    return *this;
  }

  inline bool operator==(const Rational_function_integral_mp& rhs) { 
    return (numerator_ == rhs.numerator_) && (denominator_ == rhs.denominator_);
   }
  inline bool operator!=(const Rational_function_integral_mp& rhs) { return !((*this) == rhs); }

  inline bool operator==(const int& rhs) { 
    Rational_function_integral_mp rf_rhs(rhs);
    return (*this) == rf_rhs;
   }
  inline bool operator!=(const int& rhs) { return !((*this) == rhs); }
/* @} */

/** \name Model of ScalarGroupOperations
 * 
 * @{ */
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

  /**
   * @brief      Multiplication by an integer (Z-module structure).
   *
   * @param[in]  rhs   The right hand side
   *
   * @return     The result of the multiplication assignment
   */
  Rational_function_integral_mp& operator*=(Coefficient rhs)
  {                 
    numerator_ *= rhs;
    // denominator_ *= rhs.denominator_;
    normalize();             
    return *this; 
  }
/* @} */

/** \name Model of ScalarRingOperations
 * 
 * @{ */
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
/* @} */

/** \name Model of ScalarFieldOperations
 * 
 * @{ */
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
/* @} */

  Polynomial numerator() const { return numerator_; }
  Polynomial denominator() const { return denominator_; }

  std::string to_string() {
    std::stringstream ss;
    ss << numerator_ << "/" << denominator_;
    return ss.str();
  }

private:
//make sure the numerator and denominator are coprime polynomials, i.e., gcd = integer
  void normalize() {
    auto gcd_poly = boost::math::tools::gcd(numerator_,denominator_);
    numerator_ /= gcd_poly;
    denominator_ /= gcd_poly; 
  }


private:
  Polynomial numerator_;
  Polynomial denominator_;
};

  

std::ostream& operator<<(std::ostream& os, const Rational_function_integral_mp& x)
{
  auto num = x.numerator();
  auto den = x.denominator();
  for(int n=(int)num.size()-1; n>-1; --n) {
    if(num[n] != 0) { 
      if(n == 0) {
        os << num[n];
      }
      else {
        if(num[n] == 1) { os << "X^" << n ; }
        else { os << num[n] << " X^" << n ; }
        os << " + ";
        // if(n != (int)num.size()-1) { os << " + "; } 
      }
    }
  }
  os << " / ";
  for(int n=(int)den.size()-1; n>-1; --n) {
    if(den[n] != 0) {
      if(n == 0) { 
        if(den[n] < 0) { os << " " << den[n]; }
        else { os << " + " << den[n]; }
      }
      else { //n>0
        if(den[n] == 1) { 
          if(n == (int)den.size()-1) { os << "X^" << n ; }
          else { os << " + X^" << n ; }
        }
        else { //den[n] != 0
          if(n == (int)den.size()-1) { 
            os << den[n] << " X^" << n ; 
          }
          else { 
            if(den[n] < 0) { os << den[n] << " X^" << n ; }
            else { os << " + " << den[n] << " X^" << n ; }
          }
        }
      }
    }
  }
  // os << x.numerator() << "/" << x.denominator();
  return os;
}

}  //namespace kumquat

#endif //KUMQUAT_RATIONAL_FUNCTION_INTEGRAL_MP_H_