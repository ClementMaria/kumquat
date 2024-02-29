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

/** \brief The concept for a scalar with standard field operations, where the underlying ring is denoted (R,+,*).
 * 
 * Same as AlgebraicElementPrincipalIdealDomain, with the garantee that operator%= always returns 0, the additive identity of the field. Inverses can be computed with 1/x.
  * */
class AlgebraicElementField : AlgebraicElementPrincipalIdealDomain {
};

} // namespace kumquat