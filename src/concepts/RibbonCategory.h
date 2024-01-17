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

/** Concept for a ribbon category (algebra). 
 * 
 * Introduces the notion of handle for objects and morphism to name specific objects and morphisms, that may be pre-computed in practice.
 */
struct RibbonCategory : public MonoidalCategory {
/** \brief A handle type to designate an object in the category.*/
  typedef unspecified Object_handle;
/** \brief A handle type to designate a morphism in the category.*/
  typedef unspecified Morphism_handle;
/** \brief Returns the braiding \f$c_{V,W}: V\otimes W \to W \otimes V\f$.
 * 
 * Input the handle for objects \f$V\f$ and \f$W\f$.*/  
  Morphism_handle braiding(Object_handle v_h, Object_handle w_h);
/** \brief Returns the inverse braiding \f$c^{-1}_{V,W}: W \otimes V \to V \otimes W\f$.
 * 
 * Input the handle for objects \f$V\f$ and \f$W\f$.*/
  Morphism_handle braiding_inv(Object_handle v_h, Object_handle w_h);
/** \brief Returns the twist morphism \f$\theta_V: V \to V\f$.
 * 
 * Input the handle for object \f$V\f$.*/
  Morphism_handle twist(Object_handle v_h);
/** \brief Returns the inverse of the twist morphism \f$\theta^{-1}_V: V \to V\f$.
 * 
 * Input the handle for object \f$V\f$.*/
  Morphism_handle twist_inv(Object_handle v_h);
/** \brief Return the dual of an object \fV^*\f$.
 * 
 * Input the handle for object \f$V\f$.*/
  Object_handle dual(Object_handle v_h);
/** \brief Return the pairing morphism \f$d_V: V^* \otimes V \to \mathbbm{1}\f$.
 * 
 * Input the handle for object \f$V\f$.*/
  Morphism_handle pairing(Object_handle v_h);
/** \brief Return the copairing morphism \f$b_V: \mathbbm{1} \to V \otimes V^* \to \f$.
 * 
 * Input the handle for object \f$V\f$.*/
  Morphism_handle copairing(Object_handle v_h);
/** \brief Return the dimension of an object.
 *
 * Input the handle for object \f$V\f$, return an element of \f$\operatorname{End}(\mathbbm{1})\f$.*/
  Base_element dim(Object_handle v_h);
/** Return the trace of a morphism.*/
  Base_element trace(const Morphism& phi);
};  

} //namespace kumquat