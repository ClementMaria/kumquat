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

/** \brief The concept for an element in a set.
 * 
 * An algebraic element in \kumquat is a mathematical object on which C++ operators are defined, that correspond to algebraic operations (such as addition, multiplication, etc). Algebraic elements are often described by the type of algebraic structures they belong too (such as a group, an abelian group, a ring, etc). 
 * 
 * The operations (+, -, *, /, etc) intrinsically defined in an algebraic element are the ones of the natural algebraic structure they belong too (an integer in \f$\mathbb{Z}\f$, a float in \f$\mathbb{R}\f$, a rational in \f$\mathbb{Q}\f$, etc).  
 * */ 
class AlgebraicElementSet {
/** \brief Default constructor.*/
  AlgebraicElementSet();
/** \brief Copy constructor.*/
  AlgebraicElementSet(const AlgebraicElementSet& other);
/** \brief Move constructor.*/
  AlgebraicElementSet(AlgebraicElementSet&& other) noexcept;
/** \brief Destructor.*/
  ~AlgebraicElementSet();
/** \brief Copy assignment.
 * 
 * return *this;
 * */
  AlgebraicElementSet& operator=(const AlgebraicElementSet& other);
/** \brief Move assignment.
 * 
 * return *this;
 * */
  AlgebraicElementSet& operator=(AlgebraicElementSet&& other) noexcept;
/** \brief Test for equality.*/
  inline bool operator==(const AlgebraicElementSet& lhs, const AlgebraicElementSet& rhs);
/** \brief Test whether two elements are different.
 * 
 * Based on ==.
 * */
  inline bool operator!=(const AlgebraicElementSet& lhs, const AlgebraicElementSet& rhs);
};

} // namespace kumquat