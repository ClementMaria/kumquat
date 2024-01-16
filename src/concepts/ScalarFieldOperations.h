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

/** \brief The concept for a scalar with standard field operations, where the underlying ring is denoted (R,+,*).
  * */
class ScalarFieldOperations : ScalarRingOperations {
  /** \brief Division to the right.
   * 
   * Set *this <- *this / rhs.
   * return *this;
   * */
  ScalarFieldOperations& operator/=(const ScalarFieldOperations& rhs);
  /** Return the division of one field element to the other.
   * 
   * Return (lhs / rhs), based on operator /=.
   * */
  friend ScalarFieldOperations operator/(ScalarFieldOperations lhs, const ScalarFieldOperations& rhs);
};

} // namespace kumquat