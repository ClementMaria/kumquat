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

/** Concept for an integer type. 
 * 
 * IntegerNumber must be signed.
 * 
 * IntegerNumbers must implement all usual operators +, -, *, /, +=, -=, *=, /=, %, <, >, <=, >=, ==.
 *
 * (int) must be convertible to (IntegerNumber).
 * 
 * IntegerNumber must be a valid instanciation of the template of methods 
 * boost::integer::extended_euclidean<NumberType>(NumberType, NumberType)
 * boost::integer::gcd<NumberType>(NumberType, NumberType)
. */
struct IntegerNumber {
};

} //namespace kumquat
