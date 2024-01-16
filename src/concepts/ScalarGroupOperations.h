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
struct ScalarGroupOperations : public ScalarSetOperations {
  /** \brief Addition to the right.
   * 
   * Set *this <- *this + rhs.
   * return *this;
   * */
  ScalarGroupOperations& operator+=(const ScalarGroupOperations& rhs); 
  /** Return the addition of two group elements.
   * 
   * Return (lhs + rhs), based on operator +=.
   * */
  friend ScalarGroupOperations operator+(ScalarGroupOperations lhs, const ScalarGroupOperations& rhs);
  /** \brief Substraction to the right.
   * 
   * Set *this <- *this - rhs.
   * return *this;
   * */
  ScalarGroupOperations& operator-=(const ScalarGroupOperations& rhs);
  /** Return the substraction of an element to the other.
   * 
   * Return (lhs - rhs), based on operator -=.
   * */
  friend ScalarGroupOperations operator-(ScalarGroupOperations lhs, const ScalarGroupOperations& rhs); 
};

} //namespace kumquat