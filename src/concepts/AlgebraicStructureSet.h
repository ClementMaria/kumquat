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

class AlgebraicStructureSet {
/** Constructor.*/
  AlgebraicStructureSet();
/** Destructor.*/
  ~AlgebraicStructureSet();
/** \brief Copy constructor.*/
  AlgebraicStructureSet(const AlgebraicStructureSet& other);
/** \brief Move constructor.*/
  AlgebraicStructureSet(AlgebraicStructureSet&& other) noexcept;
/** \brief Copy assignment.
  * 
  * return *this;
  * */
  AlgebraicStructureSet& operator=(const AlgebraicStructureSet& other);
/** \brief Move assignment.
  * 
  * return *this;
  * */
  AlgebraicStructureSet& operator=(AlgebraicStructureSet&& other) noexcept;
/** \brief Test for equality.*/
  inline bool operator==(const AlgebraicStructureSet& lhs, const AlgebraicStructureSet& rhs);
/** \brief Test whether two element are different.
  * 
  * Based on ==.
  * */
  inline bool operator!=(const AlgebraicStructureSet& lhs, const AlgebraicStructureSet& rhs);
/** The type of elements in the AlgebraicStructureSet. Must be a model of ScalarAlgebraicStructureSetOperations. */
  typedef ScalarAlgebraicStructureSetOperations Element;
/** Return true iff a and b represent the same element of the AlgebraicStructureSet.*/
  bool equal(Element a, Element b);
};

} // namespace kumquat