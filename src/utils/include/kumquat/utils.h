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

#ifndef KUMQUAT_UTILS_H_ 
#define KUMQUAT_UTILS_H_

// #include <kumquat/Q.h>
// #include <kumquat/Q_U1.h>
// #include <kumquat/Q_mod_Z.h>
// #include <kumquat/Z.h>
// #include <kumquat/Z_mod_nZ.h>

namespace kumquat {

// /** Finite cyclic group \bigoplus_{k,n} \left(Z/p^k Z\right)^n, of elements of order a power of a fixed prime p.
//  * 
//  **/
// template< typename T>
// std::string to_string(typename T::Element x) {
//   return std::to_string(x);
// } 

// template< typename IntegerNumber > 
// std::string to_string< Q<IntegerNumber> >(typename Q<IntegerNumber>::Element q) {
//   return to_string(q.first) + "/" + to_string(q.second);
// }

// template< typename IntegerNumber > 
// std::string to_string<Q_mod_Z<IntegerNumber> >(typename Q_mod_Z<IntegerNumber>::Element q) {
//   return to_string(q.first) + "/" + to_string(q.second);
// }

// template< typename IntegerNumber > 
// std::string to_string< Q_U1<IntegerNumber> >(typename Q_U1<IntegerNumber>::Element q) {
//   return to_string< Q<IntegerNumber> >(q.first) + ". e(" + to_string< Q_mod_Z<IntegerNumber> >(q.second) + ")";
// }

// template<> 
// std::string to_string< Z_mod_nZ >(typename Z_mod_nZ::Element m) {
//   return to_string(m);
// }

// template< typename IntegerNumber > 
// std::string to_string< Z<IntegerNumber> >(typename Z<IntegerNumber>::Element z) {
//   return to_string(z);
// }



} //namespace kumquat

#endif // KUMQUAT_X_XX_H_
