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

/** Concept for a monoidal category (algebra). 
 * 
 * A MonoidalCategory defines the notion of tensor product.
 */
struct StrictMonoidalCategory : public Category {
  /** \brief Tensor product of objects, to the right.
   * 
   * Set \f$\mathtt{lhs <- lhs} \otimes \mathtt{rhs}\f$.
   * */
  void rtensor_equal(Object &lhs, const Object& rhs);
  /** \brief Tensor product of objects.
   * 
   * Return \f$\mathtt{lhs} \otimes \mathtt{rhs}\f$, based on \f$\mathtt{rtensor_equal}\f$.
   * */
  Object tensor(Object lhs, const Object& rhs);
  /** \brief Tensor product of objects, to the left.
   * 
   * Set \f$\mathtt{rhs <- lhs} \otimes \mathtt{rhs}\f$.
   * */
  void ltensor_equal(Object& rhs, const Object& lhs);
  /** \brief Tensor product of morphisms, to the right.
   * 
   * Set \f$\mathtt{lhs <- lhs} \otimes \mathtt{rhs}\f$.
   * */
  void rtensor_equal(Morphism& lhs, const Morphism& rhs);
  /** \brief Tensor product of morphisms.
   * 
   * Return \f$\mathtt{lhs} \otimes \mathtt{rhs}\f$, based on \f$\mathtt{rtensor_equal}\f$.
   * */
  Morphism tensor(Morphism lhs, const Morphism& rhs);
  /** \brief Tensor product of morphisms, to the left.
   * 
   * Set \f$\mathtt{rhs <- lhs} \otimes \mathtt{rhs}\f$.
   * */
  void ltensor_equal(Morphism& rhs, const Morphism& lhs);

};  

} //namespace kumquat