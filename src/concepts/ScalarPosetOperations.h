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

/** \brief The concept for a scalar with standard poset operations, where the poset is denoted by (P,<).
 * 
 * The operators must satisfy the following:
 * if (a<b) and (b>a) iff a and b are not comparable,
 * if !(a<b) and !(b<a) iff (a==b) (set operation).
 */ 
class ScalarPosetOperations : public ScalarSetOperations{
  /** \brief Return true iff the two elements are comparable.*/
  bool comparable(const ScalarSetOperations& rhs);
  /** \brief Strictly less than.*/
  inline bool operator<(const ScalarSetOperations& lhs, const ScalarSetOperations& rhs);
  /** \brief Less or equal than.*/
  inline bool operator<=(const ScalarSetOperations& lhs, const ScalarSetOperations& rhs);
  /** \brief Strictly greater than.*/
  inline bool operator>(const ScalarSetOperations& lhs, const ScalarSetOperations& rhs);
  /** \brief Greater or equal than.*/
  inline bool operator>=(const ScalarSetOperations& lhs, const ScalarSetOperations& rhs);
};

} // namespace kumquat