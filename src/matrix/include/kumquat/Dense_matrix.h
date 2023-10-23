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
/** \brief Copy constructor. */
  Dense_matrix(Dense_matrix & other) : n_(other.n_), m_(other.m_), mat_(other.mat_), G_(other.G_), dim_ker_(other.dim_ker_), dim_im_(other.dim_im_) {}

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
  void times_equal_column(size_t i, Coefficient z) {
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
      auto it_conflict = low_to_col_idx.find(curr_low);//is it already taken?
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
        col_ops.emplace_back(i,it_conflict->second,z);//memorize operation
        curr_low = low_col(i,curr_low);//new value for lowest index
        it_conflict = low_to_col_idx.find(curr_low);////is there a conflict?
      }
      
      if(curr_low != -1) {//col_i != 0 after reduction
        low_to_col_idx[curr_low] = i;//low_to_col_idx does not contain -1
      }
      else { ++dim_ker_; }
    }
    dim_im_ = num_columns() - dim_ker_; 
  }
/** \brief Reduce the matrix to row echelon form.
 * 
 * Keep track of the sequence of operation as a vector of 
 * tuple (i,j,z), meaning operation row_i <- row_i + z * row_j, 
 * pushed into the input row_ops.
 * 
 * The Dense_matrix is in row echelon form after the operation
 * 
 * The coefficient structure must be a field.
 * */
  void row_echelon_form(std::vector< 
                              std::tuple<size_t,size_t,Coefficient> > & row_ops) {
    // dim_ker_ = 0;
    dim_im_ = 0;
    //low_to_row_idx[low] = i indicates that non-trivial reduced row_i has rightmost non-zero element at index low
    std::map<size_t, size_t> low_to_row_idx;
    //track the column operations in the appropriate order
    for(size_t i=0; i<num_rows(); ++i) {
      int curr_low = low_row(i);//rightmost non-trivial index for row_i
      
      std::cout << " low_row = " << curr_low << "\n";

      auto it_conflict = low_to_row_idx.find(curr_low);//is it already taken?
      while(it_conflict != low_to_row_idx.end()) {
        //found a row with index < i and with same rightmost index
        // z = - row_i[low] / row_j[low]
        Coefficient z = G_.additive_inverse(
                          G_.times(mat_[i][curr_low],
                              G_.multiplicative_inverse(
                                mat_[it_conflict->second][curr_low]
                                                             )));
        //row_i <_ row_i - z * row_j
        plus_equal_row(i,it_conflict->second,z);
        row_ops.emplace_back(i,it_conflict->second,z);//memorize operation
        curr_low = low_row(i,curr_low);//new value for lowest index
        it_conflict = low_to_row_idx.find(curr_low);////is there a conflict?
      }
      
      if(curr_low != -1) {//col_i != 0 after reduction
        low_to_row_idx[curr_low] = i;//low_to_col_idx does not contain -1
        ++dim_im_;
      }
      // else { ++dim_ker_; }
    }
    dim_ker_ = num_columns() - dim_im_; 
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
    if(hint < 1) { return -1; }
    size_t start = num_rows();
    if(hint < start) { start = hint; }
    for(int j = start-1 ; j >= 0; --j) {
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
                         = std::numeric_limits<size_t>::max()) {
    std::cout << "   *** low_row_" << i << "   hint=" << hint << "\n";
    if(hint < 1) { return -1; }
    size_t start = num_columns();
    if(hint < start) { start = hint; }
    for(int j = start-1 ; j >= 0; --j) {
      std::cout << "   " << j << "\n";
      if(!G_.trivial(mat_[i][j])) { return j; }
    }
    return -1;
  }
public:
/** \brief Diagonalize the Gram matrix of a bilinear form.
 * 
 * The matrix mat_ must contain the values b(x_i,x_j) of a bilinear form from a finite abelian group G = <x_1, ... x_n> to the abelian group \f$\mathbb{Q}_{(p)}/\mathbb{Z}\f$, for a prime number \f$p \geq 2\f$. In particular, the matrix is square symmetric.
 * */
  void diagonalize_gram_matrix() {
    //check Matrix_reduction
  }

public:
  enum Operation_type {
    col_plus_eq,//(i,j,z)     col_i <- col_i + z*col_j 
    row_plus_eq,//(i,j,z)     row_i <- row_i + z*row_j
    col_exch,//(i,j,1)        col_i <-> col_j 
    row_exch,//(i,j,1)        row_i <-> row_j 
    col_times_eq,//(i,i,z)    col_i <- z*col_i 
    row_times_eq };//(i,i,z)  row_i <- z*row_i
/** \brief Computes the Smith normal form of a matrix with coefficients in a PID.
 * 
 * After computai
 * */
  void smith_normal_form() {
  // % std::vector< std::tuple<size_t,size_t,Coefficient,Operation_type> > & ops
                          // {
    for(size_t i=0; i < num_columns(); ++i) {
  //at this stage, the top left corner up to i-1 is diagonal
      int c_idx = -1;//idx for a non-trivial column, 0 if not found
      int r_idx = -1;//idx of any non-trivial elem of col_{c_idx}
      for(size_t j=i; j<num_columns(); ++j) {
  //check whether col_j is non-trivial
        r_idx = trivial_column(j);
        if( r_idx != -1 ) { //found a non-trivial column
          c_idx = j; break;
        }
      }
  //if no non-trivial column has been found, we are done
      if(c_idx == -1) { return; }
  //o.w., col_i to col_{c_idx-1} are trivial, and col_{c_idx}[0...r_idx-1] is trivial
  //enforce col_{c_idx}[r_idx] divides all other non-trivial elements of its column and its row
  //first the column
      for(size_t k = r_idx+1; k<num_rows(); ++k) {
  //call x=mat_[r_idx][c_idx] and y=mat_[k][c_idx]
        if(!G_.trivial(mat_[k][c_idx])) {//if y!=0
  //write u*x + v*y = gcd(x,y), with gcd > 0
          auto u_v_gcd = G_.extended_gcd(mat_[r_idx][c_idx], 
                                         mat_[k][c_idx]   );
  //if x does not divide y, i.e., |gcd| < |x|
          if( std::get<2>(u_v_gcd) < G_.abs(mat_[r_idx][c_idx]) ) {
  //perform row_{r_idx} <- u * row_{r_idx} + v * row_k
            times_equal_row(r_idx, std::get<0>(u_v_gcd));
            plus_equal_row(r_idx, k, std::get<1>(u_v_gcd));
          }//now x <- gcd(x,y) and new_x divides y
        }
      }
  //then the row
      for(size_t k = c_idx+1; k<num_columns(); ++k) {
  //call x=mat_[r_idx][c_idx] and y=mat_[r_idx][k]
        if(!G_.trivial(mat_[r_idx][k])) {//if y!=0
  //write u*x + v*y = gcd(x,y), with gcd > 0
          auto u_v_gcd = G_.extended_gcd(mat_[r_idx][c_idx], 
                                         mat_[r_idx][k]   );
  //if x does not divide y, i.e., |gcd| < |x|
          if( std::get<2>(u_v_gcd) < G_.abs(mat_[r_idx][c_idx]) ) {
  //perform col_{c_idx} <- u * col_{c_idx} + v * col_k
            times_equal_column(c_idx, std::get<0>(u_v_gcd));
            plus_equal_column(c_idx, k, std::get<1>(u_v_gcd));
          }//now x <- gcd(x,y) and new_x divides y
        }
      }
  //now mat_[r_idx][c_idx] divides all entries in its column and row -> cancel the column then the row
      for(size_t k = r_idx+1; k<num_rows(); ++k) {
  //call x=mat_[r_idx][c_idx] and y=mat_[k][c_idx]
        if(!G_.trivial(mat_[k][c_idx])) {//if y!=0
  //perform row_k <- y/x * row_{r_idx}
          plus_equal_row(k, r_idx, G_.times(G_.division(mat_[k][c_idx],mat_[r_idx][c_idx]),-1));
        }
      }
  //we can directly trivialize the row r_idx now
      for(size_t k = c_idx+1; k<num_columns(); ++k) {
        mat_[r_idx][k] = G_.additive_identity();
      }
  //put the new divisor on the diagonal
      exchange_row(r_idx,i);
      exchange_col(r_idx,i);
    }
  }
private:
  /** \brief Check whether a column is trivial.*
   * 
   * Return the first row index of a non-zero element if the column is not trivial, and return -1 otherwise. hint is an optional starting index; if given, the procedure checks only whether col[hint ...] is non-trivial, assuming that col[0...hint-1] is known to be 0.
   */
  int trivial_column(size_t idx, int hint = -1) {
    size_t start;
    if(hint < 0) { start = 0; }
    else { start = hint; }
    for(size_t i=start; i<num_rows(); ++i) {
      if(!G_.trivial(mat_[i][idx])) { return i; }
    }
    return -1;
  }

private:
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
  for (size_t i = 0; i < mat.num_rows(); i++) {
    for(size_t j = 0; j < mat.num_columns(); ++j) {
      os << std::setw(7) << std::left << mat[i][j] << " ";
    }
    os << "\n";
  }
  return os;
}

}  //namespace kumquat

#endif //KUMQUAT_DENSE_MATRIX_H_