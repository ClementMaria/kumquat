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
 * As a pseudo ring with unit, the concept Ring inherits all types and 
 * methods from the concept PseudoRing. */
struct Ring : public PseudoRing {
/** Return the multiplicative identity 1.*/
    Element multiplicative_identity();
/** Return the multiplicative inverse of a Element a if it exists, and return the additive identity 0 otherwise.*/
    Element multiplicative_inverse(Element a);
};  

} //namespace kumquat