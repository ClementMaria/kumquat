/*    This file is part of the KUMQUAT Library - https://kumquat.inria.fr/ 
 *    - released under  GNU GENERAL PUBLIC LICENSE v.3
 *    See file LICENSE or go to https://kumquat.inria.fr/licensing/ for full 
 *    license details.
 *    Author(s):       Clément Maria
 *
 *    Copyright (C) 2023 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

namespace kumquat {

/** Concept for a ring (algebra). 
 * 
 * As a group with multiplication, the concept Ring inherits all types and 
 * methods from the concept AbelianGroup. */
struct Ring : public AbelianGroup {
/** Return the multiplicative identity 1.*/
    Element multiplicative_identity();
/** Set a <- (a*b). */
    void times_equal(Element & a, Element b);
/** Return a*b.*/
    Element times(Element a, Element b);
/** Set a <- a^p. */
    void pow_equal(Element & a, Integer p);
/** Return a <- a^p for a positive integer p. */
    Element pow(Element a, Integer p);
/** Return the multiplicative inverse of a Element a if it exists, and return the additive identity 0 otherwise.*/
    Element multiplicative_inverse(Element a);
};  

} //namespace kumquat