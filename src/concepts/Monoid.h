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

/** Concept for an abelian group (algebra), with additive notation (A,+). 
 * 
 * As a group with multiplication, the concept Ring inherits all types and 
 * methods from the concept AbelianGroup. */
struct Monoid : Set {
/** Return the additive identity 0.*/
  Element additive_identity();
/** Set a <- (a+b). */  
  void plus_equal(Element & a, Element b);
/** Return a+b.*/
  Element plus(Element a, Element b);
};

} //namespace kumquat
