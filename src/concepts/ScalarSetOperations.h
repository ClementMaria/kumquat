/*    This file is part of the KUMQUAT Library - https://kumquat.inria.fr/ 
 *    - which is a licence protected library. 
 *    See file LICENSE or go to https://kumquat.inria.fr/licensing/ for full 
 *    license details.
 *    Author(s):       Clément Maria
 *
 *    Copyright (C) 2024 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

namespace kumquat {

/** \brief The concept for a scalar in a set.
 * 
 * A scalar in \kumquat is a mathematical object on which C++ operators are defined, that correspond to algebraic operations (such as addition, multiplication, etc). Scalars are often described by the type of algebraic structures they belong too (such as a group, an abelian group, a ring, etc).
 * */ 
class ScalarSetOperations {
  /** \brief Copy assignment.
   * 
   * return *this;
   * */
  ScalarSetOperations& operator=(const ScalarSetOperations& other);
  /** \brief Move assignment.
   * 
   * return *this;
   * */
  ScalarSetOperations& operator=(ScalarSetOperations&& other) noexcept;
  /** \brief Test for equality.*/
  inline bool operator==(const ScalarSetOperations& lhs, const ScalarSetOperations& rhs);
  /** \brief Test whether two element are different.
   * 
   * Based on ==.
   * */
  inline bool operator!=(const ScalarSetOperations& lhs, const ScalarSetOperations& rhs);
};

} // namespace kumquat