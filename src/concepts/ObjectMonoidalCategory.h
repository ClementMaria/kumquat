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

/** \brief The concept for an object with monoidal category operations.
  * */
class ObjectMonoidalCategory : ScalarSetOperations {
  /** \brief Tensor product to the right.
   * 
   * Set \f$\mathtt{(*this) <- (*this)} \otimes \mathtt{rhs}\f$.
   * return *this;
   * */
  void rtensor_equal(const ScalarRingOperations& rhs);
  /** \brief Tensor product to the right.
   * 
   * Return \f$\mathtt{(*this)} \otimes \mathtt{rhs}\f$.
   * */
  void rtensor(const ScalarRingOperations& rhs);
  /** \brief Tensor product to the left.
   * 
   * Set \f$\mathtt{(*this) <- rhs} \otimes \mathtt{(*this)}\f$.
   * return *this;
   * */
  void ltensor_equal(const ScalarRingOperations& lhs);
  /** \brief Tensor product to the left.
   * 
   * Return \f$\mathtt{lhs} \otimes \mathtt{(*this)}\f$.
   * */
  void ltensor(const ScalarRingOperations& lhs);
};

} // namespace kumquat