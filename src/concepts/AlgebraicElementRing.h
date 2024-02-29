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
class AlgebraicElementRing : AlgebraicElementGroup {
/** \brief Multiplication to the right.
* 
* Set *this <- *this * rhs.
* return *this;
* */
  AlgebraicElementRing& operator*=(const AlgebraicElementRing& rhs);
/** Return the multiplication of two ring elements.
 * 
 * Return (lhs * rhs), based on operator *=.
 * */
  friend AlgebraicElementRing operator*(AlgebraicElementRing lhs, const AlgebraicElementRing& rhs);

};

} // namespace kumquat