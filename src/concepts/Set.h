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

class Set {
/** The type of elements in the set. Must be copiable. */
  typedef unspecified Element;
/** Return true iff a and b represent the same element of the set.*/
  bool equal(Element a, Element b);
};

} // namespace kumquat