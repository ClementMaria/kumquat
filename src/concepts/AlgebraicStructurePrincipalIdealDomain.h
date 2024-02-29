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

/** Concept for a principal ideal domain (algebra). 
 * 
 * As ring with extra structure, the concept PrincipalIdealDomain inherits all 
 * types and methods from the concept Ring. */
struct AlgebraicStructurePrincipalIdealDomain : public AlgebraicStructureRing {
/** \brief Return the absolute value of an element (if ring of integers), or 
 * more generally a Euclidean function for Euclidean algorithm.
 */
  Integer abs(Element a); 
/** \brief Compute the extended greatest common divisor of two elements of 
 * the ring. 
 * 
 * Return a triple (u,v,gcd) opf ring elements such that gcd is the greatest 
 * common divisor of x and y, and (u,v) satiosfies the Bezout identity:
 * u*x + v*y = gcd, for + the ring addition and * the ring multiplication. 
 * */ 
  std::tuple<Element,Element,Element> extended_gcd(Element x, Element y);
/** \brief Compute the division x/y in the PID.
 * 
 * Return the value q such that x = q*y + r, with 0 \leq f(r) < f(y) for a 
 * Euclidean function on the Euclidean domain (for exemple, the absolute value 
 * for the integers).*/ 
  Element division(Element x, Element y) {}
/** \brief Compute the remainder of the division x/y in the PID.
 * 
 * Return the value r such that x = q*y + r, with 0 \leq r < |y|, with 
 * 0 \leq f(r) < f(y) for a 
 * Euclidean function on the Euclidean domain (for exemple, the absolute value 
 * for the integers). 
 * */ 
  Element remainder(Element x, Element y) {}
};

} //namespace kumquat