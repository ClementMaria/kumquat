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

/** Concept for a pseudo ring (algebra), i.e., a ring without unit. 
 * 
 * As a group with multiplication, the concept Ring inherits all types and 
 * methods from the concept AbelianGroup. */
struct PseudoRing : public AbelianGroup {
/** Set a <- (a*b). */
    void times_equal(Element & a, Element b);
/** Return a*b.*/
    Element times(Element a, Element b);
/** Set a <- a^p. */
    void pow_equal(Element & a, Integer p);
/** Return a <- a^p for a positive integer p. */
    Element pow(Element a, Integer p);
};  

} //namespace kumquat