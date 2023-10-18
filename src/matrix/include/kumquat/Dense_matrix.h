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

#ifndef KUMQUAT_DENSE_MATRIX_H_ 
#define KUMQUAT_DENSE_MATRIX_H_

#include <iomanip>
#include <vector>
#include <map>

namespace kumquat {

/** \class Dense_matrix Dense_matrix.h kumquat/Dense_matrix.h 
 * \brief A dense matrix type, suited for normalizing small matrices.
 * 
 * The template type represents the group of coefficients for the matrix entires, 
 * and must be a model of the concept AbelianGroup. 
 * 
 * All the values of a same row are represented contiguously in memory ; hence 
 * row operations are more efficient than column operations.
 */
template< class CoefficientStructure >
class Dense_matrix {
public:
/** The algebraic structure containing the coefficients for the matrix entires.
 * 
 * Must be a model of AbelianGroup.*/
  typedef CoefficientStructure Coeff_struct;
/** \brief The type of coefficients for the matrix entries.*/
  typedef typename Coeff_struct::Element Coefficient;
/** \brief Signed integer type. Must be a model of SignedInteger.*/
  typedef typename Coeff_struct::Integer Integer;
/** \brief Creates an n by m matrix with uninitialized coefficients. */
  Dense_matrix(size_t n, size_t m, CoefficientStructure G) 
  : n_(n), m_(m), G_(G), dim_ker_(-1), dim_im_(-1) {
    mat_ = std::vector< std::vector< Coefficient > >
                                  (n, std::vector< Coefficient >(m));
  }
/** \brief Access the vector encoding the row at index idx.*/
  std::vector< Coefficient > & operator[](size_t idx) {
    return mat_[idx];
  };
/** \brief Fills up the matrix from the values of a bi-variate function f.
 *  
 * The type Function must implement the operator:
 * Coefficient operator(size_t i, size_t j) 
 * that takes as input two indices and returns a scalar f(i,j) convertible 
 * to the type Coefficient of the matrix.
 **/
  template< typename Function >
  void fill(Function &f) {
    for(size_t i=0; i<n_; ++i) {
      for(size_t j=0; j<m_; ++j) {
        mat_[i][j] = f(i,j);
      }
    }
  }
/** \brief Set \f$\col_i \leftarrow z \times \col_i\f$. */
  void times_equal_col(size_t i, Coefficient z) {
    for(size_t k=0; k<n_; ++k) {
      G_.times_equal(mat_[k][i],z);
    }
  }
/** \brief Set row_i <- z * row_i. */
  void times_equal_row(size_t i, Coefficient z) {
    for(size_t k=0; k<m_; ++k) {
      G_.times_equal(mat_[i][k],z);
    }
  }
/** \brief Set col_i <- col_i + col_j. */
  void plus_equal_column(size_t i, size_t j) {
    for(size_t k=0; k<n_; ++k) {
      G_.plus_equal(mat_[k][i],mat_[k][j]);
    }
  }
/** \brief Set col_i <- col_i + z * col_j. */
  void plus_equal_column(size_t i, size_t j, Coefficient z) {
    for(size_t k=0; k<n_; ++k) {
      G_.plus_equal(mat_[k][i], G_.times(mat_[k][j],z) );
    }    
  }
/** \brief Set row_i <- row_i + row_j. */
  void plus_equal_row(size_t i, size_t j) {
    for(size_t k=0; k<m_; ++k) {
      G_.plus_equal(mat_[i][k],mat_[j][k]);
    }
  }
/** \brief Set row_i <- row_i + z * row_j. */
  void plus_equal_row(size_t i, size_t j, Coefficient z) {
    for(size_t k=0; k<m_; ++k) {
      G_.plus_equal(mat_[i][k], G_.times(mat_[j][k],z) );
    }    
  }
/** \brief Exchange the columns of index i and j.*/
  void exchange_col(size_t i, size_t j) {
    if(i==j) { return; }
    for(size_t k=0; k<n_; ++k) {
      std::swap(mat_[k][i],mat_[k][j]);
    }
  }
/** \brief Exchange the rows of index i and j.*/
  void exchange_row(size_t i, size_t j) {
    if(i==j) { return; }
    for(size_t k=0; k<m_; ++k) {
      std::swap(mat_[i][k],mat_[j][k]);
    }
  }

/** \brief Return the total number of rows in the matrix.*/
  size_t num_rows() const { return n_; }
/** \brief Return the total number of columns in the matrix.*/
  size_t num_columns() const { return m_; }
/** \brief Return a reference to the algebraic structure of the entires of 
 * the matrix.*/
  Coeff_struct & coefficient_structure() { return G_; }
/** \brief Reduce the matrix to column echelon form.
 * 
 * Keep trakc of the sequence of operation as a vector of 
 * tuple (i,j,z), meaning operation col_i <- col_i + z * col_j, 
 * pushed into the input col_ops.
 * 
 * The Dense_matrix is in column echelon form after the operation
 * 
 * The coefficient structure must be a field.
 * */
  void column_echelon_form(std::vector< 
                std::tuple<size_t,size_t,Coefficient> > & col_ops) {
    dim_ker_ = 0;
    //low_to_col_idx[low] = i indicates that non-trivial reduced col_i has lowest non-zero element at index low
    std::map<size_t, size_t> low_to_col_idx;
    //track the column operations in the appropriate order
    for(size_t i=0; i<num_columns(); ++i) {
      int curr_low = low_col(i);//lowest non-trivial index
      
      std::cout << "curr_low = " << curr_low << "\n";

      auto it_conflict = low_to_col_idx.find(curr_low);
      while(it_conflict != low_to_col_idx.end()) {
        //found a column with index < i and with same lowest index
        // z = - col_i[low] / col_j[low]
        Coefficient z = G_.additive_inverse(
                          G_.times(mat_[curr_low][i],
                              G_.multiplicative_inverse(
                                mat_[curr_low][it_conflict->second]
                                                             )));
        //col_i <_ col_i - z * col_j
        plus_equal_column(i,it_conflict->second,z);
        col_ops.emplace_back(i,it_conflict->second,z);//mem operation
        curr_low = low_col(i,curr_low);
        it_conflict = low_to_col_idx.find(curr_low);
      }
      if(curr_low != -1) {//col_i != 0 after reduction
        low_to_col_idx[curr_low] = i;//low_to_col_idx does not contain -1
      }
      else { ++dim_ker_; }
    }
    dim_im_ = num_columns() - dim_ker_; 
  }
/** \brief Return the dimension of the kernel of the matrix, seen 
 * as a morphism.
 * 
 * Return -1 if the kernel has not been computed. The matrix must 
 * have coefficient in a field.
 * */
  int dim_kernel() { return dim_ker_; }
/** \brief Return the dimension of the image of the matrix, seen 
 * as a morphism between vector space.
 * 
 * Return -1 if the kernel has not been computed. The matrix must 
 * have coefficient in a field.
 * */
  int dim_image() { return dim_im_; }

private:
/** Return the lowest index of a non-zero entry to the i-th column, and return -1 if the column is all 0.
 *   
 * @param[in] i the index of the column to work on,
 * @param[in] hint an optional upper bound on the indices of the column to check (everything below hint, index hint included, is known to be 0).
 * */  
  int low_col(size_t i, size_t hint 
                         = std::numeric_limits<size_t>::max()) {
    if(hint == 0) { return -1; }
    size_t start = num_rows();
    if(hint < start) { start = hint; }
    for(size_t j = start-1 ; j >= 0; --j) {
      if(!G_.trivial(mat_[j][i])) { return j; }
    }
    return -1;
  }
/** Return the lowest (i.e. rightmost) index of a non-zero entry to the i-th row, and return -1 if the row is all 0.
 * 
 * @param[in] i the index of the row to work on,
 * @param[in] hint an optional upper bound on the indices of the row to check (everything right of hint, index hint included, is known to be 0) .
 * */  
  int low_row(size_t i, size_t hint 
                         = std::numeric_limits<size_t>::infinity()) {
    size_t start = num_columns();
    if(hint < start) { start = hint; }
    for(size_t j = start-1 ; j >= 0; --j) {
      if(!G_.trivial(mat_[i][j])) { return j; }
    }
    return -1;
  }

  //number of rows of the matrix
  size_t n_;
  //number of columns of the matrix
  size_t m_;
  //encoding of the matrix as a bi-dimensional array. The element at row i and col j is in mat_[i][j]
  std::vector< std::vector< Coefficient > > mat_;
  //the group to which the coefficients of the matrix belong.
  CoefficientStructure G_;
  //dim of the kernel, -1 if not initialized
  int dim_ker_;
  //dim of the image, -1 if not initialized
  int dim_im_;
};

/** \brief Write the dense matrix in os.*/
template<class CoefficientStructure >
std::ostream & operator<<(std::ostream & os, Dense_matrix<CoefficientStructure> & mat) {
  for (int i = 0; i < mat.num_rows(); i++) {
    for(int j = 0; j < mat.num_columns(); ++j) {
      os << std::setw(7) << std::left << mat[i][j] << " ";
    }
    os << "\n";
  }
  return os;
}

}  //namespace kumquat

#endif //KUMQUAT_DENSE_MATRIX_H_