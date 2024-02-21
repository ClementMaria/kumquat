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

#ifndef KUMQUAT_RATIONAL_FUNCTION_FLINT
#define KUMQUAT_RATIONAL_FUNCTION_FLINT

#include "fmpz_poly_qxx.h"

namespace kumquat {

/** \class Rational_function_flint Rational_function_flint.h kumquat/Rational_function_flint.h 
  * \brief .
  *
  * \implements AlgebraicElementSet
  * \implements AlgebraicElementGroup
  * \implements AlgebraicElementRing
  * \implements AlgebraicElementPID
  * \implements AlgebraicElementPoset
  */
class Rational_function_flint {
private:
  typedef flint::fmpz_poly_qxx Rational_f;

public:
  /** \brief Creator from an integer value.*/
  template<typename IntegerType>
  Rational_function_flint(IntegerType& z) {
    rf_ = z;
  }
  /** \brief Creator from a pair (degree,coefficient).
   *
   * MonomialRange are ranges whose iterators have value type std::pair<int, Integer>
   */
  Rational_function_flint(std::pair<int,Integer> monom)
  {
    Rational_f num(0);
    Rational_f x("2  0 1");//monom X / 1
    //prepare the numerator num
    num += monom.second * pow(x, monom.first);
    rf_ = num;
  }
  /** \brief Creator from two ranges of monomials.
   *
   * MonomialRange are ranges whose iterators have value type std::pair<int, Integer>
   */
  template<typename MonomialRange>
  Rational_function_flint( MonomialRange& lnum, MonomialRange & lden)
  {
    Rational_f num(0);
    Rational_f den(0);
    Rational_f x("2  0 1");//monom X / 1
    //prepare the numerator num
    if(lnum.empty()) { num = 1; }
    else {
      for(auto monom : lnum) {//monom of type std::pair<int,Integer>, (d,a) -> a*X^d
        num += monom.second * pow(x, monom.first);
      }
    }
    //prepare the denominator den
    if(lden.empty()) { den = 1; }
    else {
      for(auto monom : lden) {//(d,a) -> a*x^d
        den += monom.second * pow(x, monom.first);
      }
    }
    rf_ = num/den;
  }

/** \name Implementation of AlgebraicElementSet.
 * @{ */

/** \brief Copy constructor.*/
  Rational_function_flint(const Rational_function_flint& other) {
    rf_ = other.rf_;
  }
/** \brief Move constructor.*/
  Rational_function_flint(Rational_function_flint&& other) noexcept {
    rf_ = std::move(other.rf_);
  }
/** brief Destructor.*/
  ~Rational_function_flint();
/** \brief Copy assignment.
 * 
 * return *this;
 * */
  Rational_function_flint& operator=(const Rational_function_flint& other) {
    if(this != &other) { rf_ = other.rf_; }
    return *this;
  }
/** \brief Move assignment.
 * 
 * return *this;
 * */
  Rational_function_flint& operator=(Rational_function_flint&& other) noexcept
  {
    if(this != &other) { rf_ = std::move(other.rf_); }
    return *this;
  }
/** \brief Test for equality.*/
  inline bool operator==(const Rational_function_flint& rhs) const { 
    return rf_ == rhs.rf_;
  }
/** \brief Test whether two element are different.
 * 
 * Based on ==.
 * */
  inline bool operator!=(const Rational_function_flint& rhs) const { 
    return !(lhs == rhs); 
  }
/* @} */  //end implementation AlgebraicElementSet

/** \name Implementation of AlgebraicElementGroup.
 * @{ */

/** \brief An integer type for the \f$\mathbb{Z}\f$-module structure of the group.
 * 
 * Must be convertible to int.
 */
  typedef boost::multiprecision::mpz_int Integer;
/** \brief Multiplication by an integer.
 * 
 * Set *this <- *this + rhs. Return *this;
 * */
  Rational_function_flint& operator*=(const Integer& rhs) 
  { 
    rf_ *= rhs;
    return *this;
  }
/** \brief Addition to the right.
 * 
 * Set *this <- *this + rhs. Return *this;
 * */
  Rational_function_flint& operator+=(const Rational_function_flint& rhs) 
  {
    rf_ += rhs;
    return *this;
  }
/** Return the addition of two group elements.
 * 
 * Return (lhs + rhs), based on operator +=.
 * */
  friend Rational_function_flint operator+(Rational_function_flint lhs,        
       const Rational_function_flint& rhs) 
  {
    lhs += rhs;
    return lhs;
  }
/** \brief Subtraction to the right.
 * 
 * Set *this <- *this - rhs.
 * return *this;
 * */
  Rational_function_flint& operator-=(const Rational_function_flint& rhs) 
  {
    rf_ -= rhs;
    return *this;
  }
/** Return the subtraction of an element to the other.
 * 
 * Return (lhs - rhs), based on operator -=.
 * */
  friend Rational_function_flint operator-(Rational_function_flint lhs, const Rational_function_flint& rhs) 
  {
    lhs -= rhs;
    return lhs;
  }
/* @} */  //end implementation AlgebraicElementGroup


/** \name Implementation of AlgebraicElementRing.
 * @{ */
/** \brief Multiplication to the right.
* 
* Set *this <- *this * rhs.
* return *this;
* */
  Rational_function_flint& operator*=(const Rational_function_flint& rhs) 
  {
    rf_ *= rhs.rf_;
    return *this;
  }
/** Return the multiplication of two ring elements.
 * 
 * Return (lhs * rhs), based on operator *=.
 * */
  friend Rational_function_flint operator*(Rational_function_flint lhs, const Rational_function_flint& rhs) 
  {
    lhs *= rhs;
    return lhs;
  }
/* @} */  //end implementation AlgebraicElementRing

/** \name Implementation of AlgebraicElementPID.
 * @{ */
  /** \brief Division to the right.
   * 
   * Set *this <- *this / rhs.
   * return *this;
   * */
  Rational_function_flint& operator/=(const Rational_function_flint& rhs) 
  {
    if(rhs_.rf_.is_zero()) { std::cout << "Division by zero.\n"; return *this; }
    rf_ /= rhs.rf_;
    return *this;
  }
  /** \brief Return the division of one field element to the other.
   * 
   * Return (lhs / rhs), based on operator /=.
   * */
  friend Rational_function_flint operator/(Rational_function_flint lhs, const Rational_function_flint& rhs) 
  {
    lhs /= rhs;
    return lhs;
  }
  /** \brief Remainder of the division.
   * 
   * Set *this <- *this % rhs.
   * return *this;
   * */
  Rational_function_flint& operator%=(const Rational_function_flint& rhs) {
    rf_ = 0;
    return (*this);
  }
  /** \brief Return the remainder of the division of a scalar by another.
   * 
   * Return (lhs % rhs), based on operator %=.
   * */
  friend Rational_function_flint operator%(Rational_function_flint lhs, const Rational_function_flint& rhs) {
    lhs %= rhs;
    return lhs;
  }
  /** \brief Return the value of an Euclidean function compatible with the division with remainder.
  */
  Integer Euclidean_function() {
    return 0;
  }
/* @} */  //end implementation AlgebraicElementPID

/** \name Implementation of AlgebraicElementPoset.
 * @{ */
  /** \brief Return true iff the two elements are comparable.*/
  bool comparable(const Rational_function_flint& rhs) { 
    return false;
  }
  /** \brief Strictly less than.*/
  inline bool operator<(const Rational_function_flint& lhs, const Rational_function_flint& rhs) {
    return false;
  }
  /** \brief Less or equal than.*/
  inline bool operator<=(const Rational_function_flint& lhs, const Rational_function_flint& rhs) {
    return false;
  }
  /** \brief Strictly greater than.*/
  inline bool operator>(const Rational_function_flint& lhs, const Rational_function_flint& rhs) {
    return false;
  }
  /** \brief Greater or equal than.*/
  inline bool operator>=(const Rational_function_flint& lhs, const Rational_function_flint& rhs) {
    return false;
  }
/* @} */  //end implementation AlgebraicElementPoset

  // std::string to_string() {}


//Return true iff the denominator is a monomial X^d
  bool is_laurent() const {
    auto deg = rf_.den().degree();//degree of denominator
    for(int i=0; i<deg; ++i) {
      if(val_.den().get_coeff(i) != 0) { return false; }
    }
    if(val_.den().get_coeff(deg) != 1) { return false; }
    return true;
  }

private:
  Rational_f rf_;
};

} // namespace kumquat

#endif // KUMQUAT_RATIONAL_FUNCTION_FLINT

