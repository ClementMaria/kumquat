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

#ifndef KUMQUAT_FREE_H_ 
#define KUMQUAT_DENSE_MATRIX_H_

#include <iomanip>
#include <vector>
#include <map>
#include <kumquat/number_theory.h>
#include <kumquat/Q_U1.h>

namespace kumquat {

/** \class Free_module_finite_rank Free_module_finite_rank.h kumquat/Free_module_finite_rank.h 
 * \brief A class for finite rank free modules, or finite dimensional vector spaces.
 * 
 * The template type represents the ring of coefficients for the module, or the field of coefficients for vector spaces. CoefficientStructure is a model of concept Ring.
 * 
 * Morphisms are represented by dense matrices.
 * 
 */
template< class CoefficientStructure >
class Free_module_finite_rank {

  typedef std::vector<Scalar> Vector;
  typedef Dense_matrix<CoefficientStructure> Matrix; 

  int rank() { return rank_; }
private:
  int rank_;//rank or dimension of the space
};