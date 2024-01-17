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
/** \brief Type of the base monoid \f$\operatorname{End}(\mathbbm{1})\f$, which must be a model of concept Monoid.
 * 
 * When the unit object is (the vector space) \f$\mathbb{C}\f$, the base monoid is \f$\operatorname{End}(\mathbb{C})\f$ (where the monoid law + is the composition of morphisms) which is isomorph the \f$\mathbb{C}\f$ (any morphism \f$\mathbb{C}\to\mathbb{C}\f$ is a multiplication by a scalar). We use \f$\mathbb{C}\f$ as base monoid in this case.
 * */ 
  typedef Monoid Base_monoid;
/** \breif The type of elements in the base monoid \f$\operatorname{End}(\mathbbm{1})\f$.*/
  typedef Base_monoid::Element Base_element;
/** \brief Return the unit object \f$\mathbbm{1}\f$ of the strict category.*/
  Object unit();
/** \brief Return the (commutative) monoid \f$\operatorname{End}(\mathbbm{1})\f$.*/
  Base_monoid base_monoid();
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