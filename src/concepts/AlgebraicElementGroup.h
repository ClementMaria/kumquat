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

/** \brief The concept for a scalar with standard group operations, where the underlying group is denoted (G,+).
  * */
struct AlgebraicElementGroup : public AlgebraicElementSet {
/** \brief An integer type.
*/
  typedef unspecified Integer;
/** \brief Addition to the right.
 * 
 * Set *this <- *this + rhs.
 * return *this;
 * */
  AlgebraicElementGroup& operator+=(const AlgebraicElementGroup& rhs); 
/** Return the addition of two group elements.
 * 
 * Return (lhs + rhs), based on operator +=.
 * */
  friend AlgebraicElementGroup operator+(AlgebraicElementGroup lhs, const AlgebraicElementGroup& rhs);
/** \brief Subtraction to the right.
 * 
 * Set *this <- *this - rhs.
 * return *this;
 * */
  AlgebraicElementGroup& operator-=(const AlgebraicElementGroup& rhs);
/** Return the subtraction of an element to the other.
 * 
 * Return (lhs - rhs), based on operator -=.
 * */
  friend AlgebraicElementGroup operator-(AlgebraicElementGroup lhs, const AlgebraicElementGroup& rhs); 
};

} //namespace kumquat