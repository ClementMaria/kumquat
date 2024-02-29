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

/** Concept for a field (algebra). 
 * 
 * A Field is a Ring where all elements have multiplicative inverse, except 0. The Ring method multiplicative_inverse must always return a non zero value for a non-zero input.
 */
struct AlgebraicStructureField : public AlgebraicStructureRing {
};  

} //namespace kumquat