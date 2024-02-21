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

#ifndef KUMQUAT_COMPLEX_NUMBER_H_ 
#define KUMQUAT_COMPLEX_NUMBER_H_ 

#include <boost/multiprecision/mpc.hpp>

namespace kumquat {

 typedef boost::multiprecision::mpc_complex Complex_number;

} //namespace kumquat

#endif //KUMQUAT_COMPLEX_NUMBER_H_
