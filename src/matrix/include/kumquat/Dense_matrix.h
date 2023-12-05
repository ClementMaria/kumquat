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
#include <kumquat/number_theory.h>

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
/** \brief Set \f$\col_i \leftarrow z \times \col_i\f$. 
 * 
 * If CoefficientStructure implements AbelianGroup, WeightType can be of type CoefficientStructure::Integer.
 * 
 * If CoefficientStructure also implements PseudoRing, WeightType can also be of type CoefficientStructure::Element (i.e., Coefficient).
 * */
  template<typename WeightType>
  void times_equal_column(size_t i, WeightType z) {
    for(size_t k=0; k<n_; ++k) {
      G_.times_equal(mat_[k][i],z);
    }
  }
/** \brief Set row_i <- z * row_i. 
 * 
 * If CoefficientStructure implements AbelianGroup, WeightType can be of type CoefficientStructure::Integer.
 * 
 * If CoefficientStructure also implements PseudoRing, WeightType can also be of type CoefficientStructure::Element (i.e., Coefficient).
 * */
  template<typename WeightType>
  void times_equal_row(size_t i, WeightType z) {
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
/** \brief Set col_i <- col_i + z * col_j. 
 * 
 * If CoefficientStructure implements AbelianGroup, WeightType can be of type CoefficientStructure::Integer.
 * 
 * If CoefficientStructure also implements PseudoRing, WeightType can also be of type CoefficientStructure::Element (i.e., Coefficient).
 * */
  template<typename WeightType>
  void plus_equal_column(size_t i, size_t j, WeightType z) {
    if(G_.trivial(z)) return;
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
/** \brief Set row_i <- row_i + z * row_j. 
 * 
 * If CoefficientStructure implements AbelianGroup, WeightType can be of type CoefficientStructure::Integer.
 * 
 * If CoefficientStructure also implements PseudoRing, WeightType can also be of type CoefficientStructure::Element (i.e., Coefficient).
 * */
  template<typename WeightType>
  void plus_equal_row(size_t i, size_t j, WeightType z) {
    if(G_.trivial(z)) return;
    for(size_t k=0; k<m_; ++k) {
      G_.plus_equal(mat_[i][k], G_.times(mat_[j][k], z) );
    }    
  }
/** \brief Exchange the columns of index i and j.*/
  void exchange_column(size_t i, size_t j) {
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
 * Keep track of the sequence of operation as a vector of 
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
    if(hint < 1) { return -1; }
    size_t start = num_columns();
    if(hint < start) { start = hint; }
    for(int j = start-1 ; j >= 0; --j) {
      if(!G_.trivial(mat_[i][j])) { return j; }
    }
    return -1;
  }

public:
  enum Operation_type {// ? means undefined
    plus_equal_column_t,//(i,j,z)  col_i <- col_i + z*col_j 
    plus_equal_row_t,//(i,j,z)     row_i <- row_i + z*row_j
    exchange_column_t,//(i,j,?)    col_i <-> col_j 
    exchange_row_t,//(i,j,?)       row_i <-> row_j 
    times_equal_column_t,//(i,?,z) col_i <- z*col_i 
    times_equal_row_t };//(i,?,z)  row_i <- z*row_i

    struct Elementary_matrix_operation {
    public:
      Elementary_matrix_operation(size_t i, size_t j, Coefficient z, Operation_type op) : i_(i), j_(j), z_(z), op_(op) {}

      size_t lhs() { return i_; }
      size_t rhs() { return j_; }
      Coefficient coefficient() { return z_; }
      Operation_type type() { return op_;}
      
      std::string to_string() {
        switch(type()) {
          case plus_equal_column_t: 
            return "col_" + std::to_string(lhs()) + " <- col_" + std::to_string(lhs()) + " + " + std::to_string((long long)coefficient()) + " * " + "col_" + std::to_string(rhs());
          case plus_equal_row_t: 
            return "row_" + std::to_string(lhs()) + " <- row_" + std::to_string(lhs()) + " + " + std::to_string((long long)coefficient()) + " * row_" + std::to_string(rhs());
            break;
          case exchange_column_t:
            return "col_" + std::to_string(lhs()) + " <-> col_" + std::to_string(rhs());
          case exchange_row_t:
            return "row_" + std::to_string(lhs()) + " <-> row_" + std::to_string(rhs());
          case times_equal_column_t:
            return "col_" + std::to_string(lhs()) + " <- " + std::to_string((long long)coefficient()) + " * col_" + std::to_string(lhs());
          case times_equal_row_t:
            return "row_" + std::to_string(lhs()) + " <- " + std::to_string((long long)coefficient()) + " * row_" + std::to_string(lhs());
          default: return "unknown";     
        }
      }

    private:
      size_t i_;
      size_t j_;
      Coefficient z_;
      Operation_type op_;
    };
/** \brief Computes the Smith normal form of a matrix with coefficients in a PID.
 * 
 * The coefficients of the matrix must be part of a PID (i.e., Coeff_struct 
 * implements PrincipalIdealDomain).
 * */
  void smith_normal_form(std::vector< Elementary_matrix_operation > & ops) {
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
            ops.emplace_back(r_idx, r_idx, std::get<0>(u_v_gcd), times_equal_row_t);
            plus_equal_row(r_idx, k, std::get<1>(u_v_gcd));
            ops.emplace_back(r_idx,k,std::get<1>(u_v_gcd),plus_equal_row_t);
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
            ops.emplace_back(c_idx,c_idx,std::get<0>(u_v_gcd),times_equal_column_t);
            plus_equal_column(c_idx, k, std::get<1>(u_v_gcd));
            ops.emplace_back(c_idx, k, std::get<1>(u_v_gcd),plus_equal_column_t);
          }//now x <- gcd(x,y) and new_x divides y
        }
      }
  //now mat_[r_idx][c_idx] divides all entries in its column and row -> cancel the column then the row
      for(size_t k = r_idx+1; k<num_rows(); ++k) {
  //call x=mat_[r_idx][c_idx] and y=mat_[k][c_idx]
        if(!G_.trivial(mat_[k][c_idx])) {//if y!=0
  //perform row_k <- y/x * row_{r_idx}
          auto z = G_.times(G_.division(mat_[k][c_idx],mat_[r_idx][c_idx]),-1);
          plus_equal_row(k, r_idx, z);
          ops.emplace_back(k, r_idx, z, plus_equal_row_t);
        }
      }
  //we can directly trivialize the row r_idx now
      for(size_t k = c_idx+1; k<num_columns(); ++k) {
        if(!G_.trivial(mat_[r_idx][k])) {//record the operation
          auto z = G_.times(G_.division(mat_[r_idx][k],mat_[r_idx][c_idx]),-1);
          ops.emplace_back(k, c_idx, z, plus_equal_column_t);
          mat_[r_idx][k] = G_.additive_identity();//set directly to 0
        }
      }
  //put the new divisor on the diagonal
      if(r_idx != (int)i) {
        exchange_row(r_idx,i);
        ops.emplace_back(r_idx, i, G_.additive_identity(), exchange_row_t);
        exchange_column(r_idx,i);
        ops.emplace_back(r_idx, i, G_.additive_identity(), exchange_column_t);
      }
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

public:
/** \brief Diagonalize the Gram matrix of a bilinear form.
 * 
 * The matrix mat_ must contain the values b(x_i,x_j) of a bilinear form from a finite abelian group G = <x_1, ... x_n> to the abelian group \f$\mathbb{Q}_{(p)}/\mathbb{Z}\f$, for a prime number \f$p \geq 2\f$. In particular, the matrix is square symmetric.
 * */
  void diagonalize_gram_matrix_Qp_mod_Z() {
    auto n = num_rows();
    auto m = num_columns();
    if(n != m) { 
      std::cerr << "Gram matrix is not square.\n"; return; 
    }
    auto p = G_.p();

    if((p % 2) == 0) { diagonalize_gram_matrix_Qp_mod_Z_p_even(); }
    else { diagonalize_gram_matrix_Qp_mod_Z_p_odd(); }

  }

private:
//diagonalize the gram matrix in Qp_mod_Z with p odd 
  void diagonalize_gram_matrix_Qp_mod_Z_p_odd() {
    auto p = G_.p(); auto n = num_rows();
    //after q iterations, M[0..q][0..q] is block diagonal
    for(size_t num_iteration = 0; num_iteration < n; ++num_iteration) {
      //check if the bottom right block B is uniformly 0, and if not compute the minimum r such that p^r B = 0;
      int pivot_i, pivot_j;//find element of maximal order
      find_pivot_Qp_mod_Z(pivot_i,pivot_j,num_iteration);

      if(pivot_i == -1) { return; }//the remaining matrix is trivial
      //now, the matrix is non trivial, element B[i][j] has minimal order, and if there is a diagonal coefficient of minimal order, then i==j
      if(pivot_i == pivot_j) {//put on top left corner
          exchange_row(num_iteration,pivot_i);
          exchange_column(num_iteration,pivot_i);
      }
      else {//pivot_i != pivot_j and B[i][j]==B[j][i] has strictly larger order than B[i][i] and B[j][j]
        //put element b_ii + b_jj +2*b_ij in B[i][i]
        plus_equal_row(pivot_i,pivot_j);  
        plus_equal_column(pivot_i,pivot_j);
        //now, if p odd, B[i][i] has order the minimal order, and is in the diagonal
        exchange_row(num_iteration,pivot_i);
        exchange_column(num_iteration,pivot_i);
      }
      //use the pivot to cancel row[num_iteration...m] and col[num_iteration...n] 
      cancel_row_column_Qp_mod_Z(num_iteration);
    }
  }

//diagonalize the gram matrix in Qp_mod_Z with p == 2 
  void diagonalize_gram_matrix_Qp_mod_Z_p_even() {
    auto p = G_.p(); auto n = num_rows();
    //after q iterations, M[0..q][0..q] is block diagonal
    for(size_t num_iteration = 0; num_iteration < n; ++num_iteration) {
      std::cout << " ### iteration " << num_iteration << "\n";
      //check if the bottom right block B is uniformly 0, and if not compute the minimum r such that p^r B = 0;
      int pivot_i, pivot_j;
      find_pivot_Qp_mod_Z(pivot_i, pivot_j, num_iteration);

      if(pivot_i == -1) { return; }//the remaining matrix is trivial
      
      //now, the matrix is non trivial, element B[i][j] has minimal order, and if there is a diagonal coefficient of minimal order, then i==j
      if(pivot_i == pivot_j) {//put on top left corner
        exchange_row(num_iteration,pivot_i);
        exchange_column(num_iteration,pivot_i);
        std::cout << "     pivot in diagonal " << pivot_i << " " << pivot_j << "\n";
        //use the pivot to cancel row[num_iteration...m] and col[num_iteration...n] 
        cancel_row_column_Qp_mod_Z(num_iteration);
      }
      else {//pivot_i != pivot_j and B[i][j]==B[j][i] has strictly larger order than B[i][i] and B[j][j]
        std::cout << "     pivot NOT in diagonal. Pivot = " << pivot_i << " " << pivot_j << "\n";
        // std::cout << "       row_" << pivot_i << " += row_" << pivot_j << "\n"; 
        // std::cout << "       col_" << pivot_i << " += col_" << pivot_j << "\n"; 
      
        //put the future block in top left corner
        exchange_row(num_iteration,pivot_i);
        exchange_column(num_iteration,pivot_i);
        exchange_row(num_iteration+1,pivot_j);
        exchange_column(num_iteration+1,pivot_j);

        std::cout << "put block on top left\n";
        std::cout << "\n--------------\n" << *this << "\n--------------\n";

        cancel_2_2_block_Q2_mod_Z(num_iteration);
      }
    }
  }

private:
/* Find the element in the submatrix M[idx...n][idx...n] with maximal order in Qp_mod_Z.
If such element is in the diagonal, always return a diagonal element. If the matrix is uniformly 0, the pivot indices are returned as -1.
*/
  void find_pivot_Qp_mod_Z(int & pivot_i, int & pivot_j, int idx=0) {
    Integer max_order = 0;
    pivot_i = -1; pivot_j = -1;//find element of minimal order
    for(size_t i = idx; i < num_rows(); ++i) {//try to find a pivot in the diagonal
      for(size_t j = idx; j < num_columns(); ++j) {
        if(!G_.trivial(mat_[i][j])) {//!= 0
          auto ord = G_.order(mat_[i][j]);
          if(ord > max_order) { 
            max_order = ord; pivot_i = i; pivot_j = j; 
          }
          else{//prioritize diagonal element of min order
            if((ord == max_order) && (i==j) ) { 
              pivot_i = i; pivot_j = j; 
            }
          }
        }
      }
    }
    if(max_order == 0) { pivot_i = -1; pivot_j = -1; }//matrix is 0
  }
  /* Cancel the row and column of given index idx assuming that the element M[idx][idx] has maximal order in Qp_mod_Z over all other eleemnts in M[idx...n][idx...n]. Subroutine of diagonalization of Gram matrices.*/ 
  void cancel_row_column_Qp_mod_Z(size_t idx) {
        auto p = G_.p(); auto n = num_rows();
    //top left element M[idx][idx] has minimal order
    // and is written as top_left == a p^{-r} with gcd(a,p)==1 and 0 < a < p
    Coefficient top_left = mat_[idx][idx];
    //now, enforce top left corner B[idx][idx]== eps/p^k, for 
    //eps = 1 or eps quadratic non-residue mod p. Compute also eps^{-1} mod p^r
    Integer eps, eps_inv;
    
    std::cout << "now, enforce top left corner == eps/p^k\n";

    //if a quadratic residue mod p, i.e., a = x^2 mod p -> set eps = 1
    if(quadratic_residue(top_left.first,top_left.second,p)) {//compute a solution x
      std::cout << " a = " << top_left.first << " is quadquad res mod p^k=" << top_left.second <<"\n";

      auto x = solve_quadratic_residue(top_left.first,top_left.second,p);
     
      std::cout << "   x=" << x << "  s.t. x*x cong a mod p\n";

      //compute s = x^{-1} mod p^k, which exists because gcd(a,p)=1 => gcd(x,p)=1
      auto s = inverse(x,top_left.second);

      std::cout << "   s=" << s << "   s.t. s=x^{-1} mod p^k (=" << top_left.second << "\n";

      //set row_idx <- s*row_idx and col_idx <- s*col_idx such that B[idx][idx] = 1/p^r mod Z
      times_equal_row(idx,s);
      times_equal_column(idx,s);
      eps = 1; eps_inv = 1;

      std::cout << " now top_left = 1/p^k  -> " << mat_[idx][idx] <<"\n";
    }
    else {//else quadratic non-residue -> top_left is already a/p^k, set eps=a

      std::cout << " a = " << top_left.first << " is quad NON - res mod p=" << p <<"\n";

      eps = top_left.first;
      eps_inv = kumquat::inverse(eps, top_left.second);//in number_theory.h

      std::cout << "   eps=" << eps << "   and eps_inv = " << eps_inv << "(mod p^k=" << top_left.second << "\n";
    }

    top_left = mat_[idx][idx];
//we do have M[idx][idx]== eps/p^k, for eps = 1 or 
//eps quadratic non-residue mod p, and eps_inv = eps^-1 mod p^r
    for(size_t i=idx+1; i<n; ++i) {
//compute beta such that B[idx][i]== beta * b / p^r with beta > 1
// G_.p_normalize(mat_[idx][i]);
      Integer beta = kumquat::division(top_left.second, mat_[idx][i].second);//in number_theory.h
                              
      std::cout << " at B_" << idx <<","<<i<<"   beta=" << beta << "  s.t. " << mat_[idx][i] << "== beta * b / p^r " << "\n";

      //do row_i <- row_i - beta * b * eps^-1 * row_idx and
      //   col_i <- col_i - beta * eps^-1 * col_idx
      auto z = kumquat::times(
                 kumquat::times( 
                   kumquat::times( beta, 
                                   mat_[idx][i].first ) 
                   , eps_inv )
                 , (Integer(-1)));

      std::cout << "  set row_" << i << " += " << z << " * row_" << idx << "\n";
      std::cout << "  set col_" << i << " += " << z << " * col_" << idx << "\n";

      std::cout << "   cancel with: " << G_.times(mat_[idx][idx],z).first << "/" << G_.times(mat_[idx][idx],z).second << "\n";

      plus_equal_row(i,idx,z);
      plus_equal_column(i,idx,z);
    
      std::cout << "   result:\n";
      std::cout << *this << "\n";
    }
  }    


/* Subroutine of the diagonalization of the Gram matrix with coefficients in Q_2 / Z, where the block B[idx,idx+1][idx,idx+1] is of the form:
*  |2a/2^m   b/2^m |   where b is odd and a,c are arbitary
*  |               |
*  |b/2^m    2c/2^m|
*/
  void cancel_2_2_block_Q2_mod_Z(size_t idx) {
    //extract integres a,b,c and 2^m
    Integer two_to_m = mat_[idx][idx+1].second;
    Integer a = mat_[idx][idx].first * (two_to_m / (2*mat_[idx][idx].second)) ;
    Integer b = mat_[idx][idx+1].first;
    Integer c = mat_[idx+1][idx+1].first * 
                              (two_to_m / ((Integer)(2)*mat_[idx+1][idx+1].second));
    
    std::cout << "  alpha = " << a << "\n";
    std::cout << "  beta  = " << b << "\n";
    std::cout << "  gamma = " << c << "\n";
    std::cout << "  2^m   = " << two_to_m << "\n";

    //d is the inverse of 4ac-b^2 modulo 2^m
    Integer d_inv = ( (Integer)4 * a * c - b * b ) % two_to_m;
    if(d_inv < 0) { d_inv += two_to_m; }

    std::cout << "   d^-1    = " << d_inv << "\n";

    Integer d = kumquat::inverse(d_inv, two_to_m);

    std::cout << "    d = " << d << "\n";

    //for all i >idx+1, cancel M[i][idx,idx+1] and M[idx,idx+1][i]
    for(size_t i=idx+2; i<num_rows(); ++i) {
    
      Coefficient u = mat_[i][idx]; Coefficient v = mat_[i][idx+1];

      std::cout << "++++++++++++++++++ i=" << i << "   (u,v)=" << u << "," << v << "\n";
      //prepare r1 = -d(2cu-bv) = -2*d*c * (u)   +    d*b * (v)
      //minus_d2c = -2*d*c;
      Integer minus_d2c = (Integer)(-2)*d*c;  

      std::cout << "   -2dc = " << minus_d2c <<"\n";  
      //db = +d*b
      Integer db = (d*b);
      std::cout << "   db = " << db <<"\n";

      std::cout << "   v = " << v << "\n";



      // std::cout << "   dbv = " << G_.times(v,db) << "\n";




      // auto xxx = G_.times(u,minus_d2c);
      // auto yyy = G_.times(v,db);
      // auto zzz = G_.plus(xxx,yyy);
      // std::cout << "   -2dcu = " << xxx <<"\n";
      // std::cout << "   dbv   = " << yyy << "\n";
      // std::cout << "   -2dcu+dbv = " << zzz << "\n";




      Coefficient minus_r1 = G_.plus( G_.times(u,minus_d2c), G_.times(v,db)  );
      //prepare r2 = -d(2av-bu) = -2*d*a * (v)   +    d*b * (u)
      Integer minus_d2a = (Integer)(-2)*d*a;
      Coefficient minus_r2 = G_.plus( G_.times(v,minus_d2a), G_.times(u,db)  );
    

      std::cout << "    r1=-d(2cu-bv) = " << minus_r1 << "     r2=-d(2av-bu) = " << minus_r2 << "\n";




      // auto aaa = G_.times(r1,mat_[idx][idx]);
      // auto bbb = G_.plus(u, aaa);
      // auto ccc = G_.times(r2,mat_[idx][idx+1]);
      // auto ddd = G_.plus( bbb, ccc );

      // std::cout << "=== " << ddd << "\n";

      std::cout << "row_" << i << " <- row_" << i << " + (" << minus_r1 << ")*row_" << idx << "\n";
      plus_equal_row(i,idx,minus_r1);


      std::cout << "col_" << i << " <- col_" << i << " + (" << minus_r1 << ")*col_" << idx << "\n";

      plus_equal_column(i,idx,minus_r1);


      std::cout << *this << "\n\n";




      std::cout << "row_" << i << " <- row_" << i << " + (" << minus_r2 << ")*row_" << idx+1 << "\n";
      plus_equal_row(i,idx+1,minus_r2);

      std::cout << "col_" << i << " <- col_" << i << " + (" << minus_r2 << ")*col_" << idx+1 << "\n";

      plus_equal_column(i,idx+1,minus_r2);

      std::cout << *this << "\n\n";

    }



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
      // os << std::setw(1) << std::left << mat[i][j] << " ";
      os << std::setw(14) << mat[i][j] ;// << " ";
    }
    os << "\n";
  }
  return os;
}




// /** \brief Write the an elementary matrix operation in os.*/
// template<class DenseMatrix >
// std::ostream & operator<<(std::ostream & os, 
//                           typename DenseMatrix::Elementary_matrix_operation & op)
// {
//   switch(op.type()) {
//   case DenseMatrix::Operation_type::plus_equal_column_t : 
//     os << "col_" << op.lhs() << " <- " << "col_" << op.lhs() << " + " << op.coefficient() << " * " << "col_" << op.rhs();
//     break;
//   case DenseMatrix::Operation_type::plus_equal_row_t :
//     os << "row_" << op.lhs() << " <- " << "row_" << op.lhs() << " + " << op.coefficient() << " * " << "row_" << op.rhs();
//     break;
//   case DenseMatrix::Operation_type::exchange_column_t :
//     os << "col_" << op.lhs() << " <-> " << "col_" << op.rhs();
//     break;
//   case DenseMatrix::Operation_type::exchange_row_t :
//     os << "row_" << op.lhs() << " <-> " << "row_" << op.rhs();
//     break;
//   case DenseMatrix::Operation_type::times_equal_column_t :
//     os << "col_" << op.lhs() << " <- " << op.coefficient() << " * " << "col_" << op.lhs();
//     break;
//   case DenseMatrix::Operation_type::times_equal_row_t :
//     os << "row_" << op.lhs() << " <- " << op.coefficient() << " * " << "row_" << op.lhs();
//     break;
//   default: os << "unknown operation type.\n";     
//   }



}  //namespace kumquat

#endif //KUMQUAT_DENSE_MATRIX_H_