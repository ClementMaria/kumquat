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

class Set {
/** Constructor.*/
  Set();
/** Destructor.*/
  ~Set();
/** \brief Copy constructor.*/
  Set(const Set& other);
/** \brief Move constructor.*/
  Set(Set&& other) noexcept;
/** brief Destructor.*/
  ~Set();
/** \brief Copy assignment.
  * 
  * return *this;
  * */
  Set& operator=(const Set& other);
/** \brief Move assignment.
  * 
  * return *this;
  * */
  Set& operator=(Set&& other) noexcept;
/** \brief Test for equality.*/
  inline bool operator==(const Set& lhs, const Set& rhs);
/** \brief Test whether two element are different.
  * 
  * Based on ==.
  * */
  inline bool operator!=(const Set& lhs, const Set& rhs);
/** The type of elements in the set. Must be a model of ScalarSetOperations. */
  typedef ScalarSetOperations Element;
/** Return true iff a and b represent the same element of the set.*/
  bool equal(Element a, Element b);
};

} // namespace kumquat