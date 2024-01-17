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

/** Concept for a monoid (algebra), with additive notation (M,+). 
 * 
 * A monoid is a set with a (not necessarily commutative) associative binary 
 * operation +, and a unit 0.*/
struct Monoid : Set {
/** Return the additive identity 0.*/
  Element additive_identity();
/** Set a <- (a+b). */  
  void plus_equal(Element & a, Element b);
/** Return a+b.*/
  Element plus(Element a, Element b);
};

} //namespace kumquat
