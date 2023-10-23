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

/** Concept for a principal ideal domain (algebra). 
 * 
 * As ring with extra structure, the concept PrincipalIdealDomain inherits all 
 * types and methods from the concept Ring. */
struct PrincipalIdealDomain : public Ring {
  /** \brief Compute the extended greatest common divisor of two elements of 
   * the ring. 
   * 
   * Return a triple (u,v,gcd) opf ring elements such that gcd is the greatest 
   * common divisor of x and y, and (u,v) satiosfies the Bezout identity:
   * u*x + v*y = gcd, for + the ring addition and * the ring multiplication. 
   * */ 
  std::tuple<Element,Element,Element> extended_gcd(Element x, Element y);
};

} //namespace kumquat
