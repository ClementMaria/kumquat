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

namespace kumquat {

/** \brief The concept for a scalar with standard ring operations, where the unerlying ring is denoted (R,+,*).
  * */
class AlgebraicElementPID : AlgebraicElementRing {
  /** \brief Return the value of an Euclidean function compatible with the division with remainder.
  */
  Integer Euclidean_function() {}
  /** \brief Division to the right.
   * 
   * Set *this <- *this / rhs.
   * return *this;
   * */
  AlgebraicElementPID& operator/=(const AlgebraicElementPID& rhs);
  /** Return the division of one field element to the other.
   * 
   * Return (lhs / rhs), based on operator /=.
   * */
  friend AlgebraicElementPID operator/(AlgebraicElementPID lhs, const AlgebraicElementPID& rhs);  
  /** \brief Remainder of the division.
   * 
   * Set *this <- *this % rhs.
   * return *this;
   * */
  AlgebraicElementPID& operator%=(const AlgebraicElementPID& rhs);
  /** Return the remainder of the division of a scalar by another.
   * 
   * Return (lhs % rhs), based on operator %=.
   * */
  friend AlgebraicElementPID operator%(AlgebraicElementPID lhs, const AlgebraicElementPID& rhs);
};

} // namespace kumquat