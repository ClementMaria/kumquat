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

/** \brief Concept for a monoidal category (algebra). 
 * 
 * A MonoidalCategory defines a tensor product.
 */
struct Category : Set {
/** \brief The type of object in the category. Must be copiable. */
  typedef unspecified Element;
/** \brief Return true iff a and b represent the same object of the category.*/
  bool equal(Element a, Element b);
/** \brief The type of arrows in the category.*/
  typedef unspecifed Morphism;
/** \brief Compose a morphism to the right.
 * 
 * Set \f$\phi \rightarrow \phi \circ \psi\f$.
 * */
  void rcompose_equal(Morphism &phi, Morphism psi);
/** \brief Compose a morphism to the left.
 * 
 * Set \f$\phi \rightarrow \psi \circ \phi\f$.
 * */
  void lcompose_equal(Morphism &phi, Morphism psi);
/** \brief Compose two morphisms.
 * 
 * Return \f$\phi \circ \psi\f$.
 * */
  Morphism compose(Morphism phi, Morphism psi);
};  

} //namespace kumquat