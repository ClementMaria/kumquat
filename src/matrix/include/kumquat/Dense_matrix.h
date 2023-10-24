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
// /** \brief Set \f$\col_i \leftarrow z \times \col_i\f$. */
//   void times_equal_column(size_t i, Coefficient z) {
//     for(size_t k=0; k<n_; ++k) {
//       G_.times_equal(mat_[k][i],z);
//     }
//   }
// /** \brief Set row_i <- z * row_i. */
//   void times_equal_row(size_t i, Coefficient z) {
//     for(size_t k=0; k<m_; ++k) {
//       G_.times_equal(mat_[i][k],z);
//     }
//   }
/** \brief Set \f$\col_i \leftarrow z \times \col_i\f$. */
  void times_equal_column(size_t i, Integer z) {
    for(size_t k=0; k<n_; ++k) {
      G_.times_equal(mat_[k][i],z);
    }
  }
/** \brief Set row_i <- z * row_i. */
  void times_equal_row(size_t i, Integer z) {
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
  void plus_equal_column(size_t i, size_t j, Integer z) {
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
  void plus_equal_row(size_t i, size_t j, Integer z) {
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
  void diagonalize_gram_matrix_Qp_mod_Z_p_odd() {
    auto n = num_rows();
    auto m = num_columns();
    if(n != m) { 
      std::cerr << "Gram matrix is not square.\n"; return; 
    }
    auto p = G_.p();

    if((p % 2) == 0) { std::cerr << "No support for even p.\n"; return; }

    //after x iterations, the columns and rows 0 ... x-1 are block diagonal
    for(size_t num_iteration = 0; num_iteration < n; ++num_iteration) {
      std::cout << " ### iteration " << num_iteration << "\n";
      //check if the bottom right block B is uniformly 0, and if not compute the minimum r such that p^r B = 0;
      Integer min_order = (Integer)std::numeric_limits<signed long long>::max();
      int pivot_i = -1; int pivot_j = -1;//find element of minimal order
      for(size_t i = num_iteration; i<n; ++i) {//try to find a pivot in the diagonal
        for(size_t j = num_iteration; j<n; ++j) {
          if(!G_.trivial(mat_[i][j])) {//!= 0
            
            auto ord = G_.order(mat_[i][j]);

            std::cout << "      min_order = " << min_order << "    order(" << mat_[i][j] << ") == " << ord << "\n";

            if(ord < min_order) { min_order = ord; pivot_i = i; pivot_j = j; }
            else{//prioritize diagonal element of min order
              if((ord == min_order) && (i==j) ) { pivot_i = i; pivot_j = j; }
            }
          }
        }
      }
      if(pivot_i == -1) { return; }//the remaining matrix is trivial
      
      std::cout << "     found pivot (" << pivot_i << "," << pivot_j << ")\n";

      //now, the matrix is non trivial, element B[i][j] has minimal order, and if there is a diagonal coefficient of minimal order, then i==j
      if(pivot_i == pivot_j) {//put on top left corner
          exchange_row(num_iteration,pivot_i);
          exchange_column(num_iteration,pivot_i);

          std::cout << "     pivot in diagonal\n";
      }
      else {//pivot_i != pivot_j and B[i][j]==B[j][i] has strictly smaller order than B[i][i] and B[j][j]

        std::cout << "     pivot NOT in diagonal\n";
        std::cout << "       row_" << pivot_i << " += row_" << pivot_j << "\n"; 
        std::cout << "       col_" << pivot_i << " += col_" << pivot_j << "\n"; 

        plus_equal_row(pivot_i,pivot_j);  
        plus_equal_column(pivot_i,pivot_j);
        //now, if p odd, B[i][i] has order the minimal order, and is in the diagonal
        exchange_row(num_iteration,pivot_i);
        exchange_column(num_iteration,pivot_i);
      }

      std::cout << "       result where top left element " << num_iteration << "," << num_iteration << " is of minimal order \n";

      std::cout << *this << "\n\n\n";

      //now, top left element B[num_iteration][num_iteration] has minimal order
      // and is written as top_left == a p^{-r} with gcd(a,p)==1 and 0 < a < p
      Coefficient top_left = mat_[num_iteration][num_iteration];
      //now, enforce top left corner B[num_iteration][num_iteration]== eps/p^k, for 
      //eps = 1 or eps quadratic non-residue mod p. Compute also eps^{-1} mod p^r
      Integer eps, eps_inv;
      //if a quadratic residue mod p, i.e., a = x^2 mod p -> set eps = 1
      if(quadratic_residue(top_left.first,p)) {//compute a solution x
        auto x = solve_quadratic_residue(top_left.first,p);
        //compute s = x^{-1} mod p^k, which exists because gcd(a,p)=1 => gcd(x,p)=1
        auto s = inverse(x,top_left.second);
        //set row_num_iteration <- s*row_num_iteration and col_num_iteration <- s*col_num_iteration such that B[num_iteration][num_iteration] = 1/p^r mod Z
        times_equal_row(num_iteration,s);
        eps = 1; eps_inv = 1;
      }
      else {//else quadratic non-residue -> top_left is already a/p^k, set eps=a
        eps = top_left.first;
        eps_inv = inverse(eps,top_left.second);//in number_theory.h
      }
      //we do have B[num_iteration][num_iteration]== eps/p^k, for eps = 1 or 
      //eps quadratic non-residue mod p, and eps_inv = eps^-1 mod p^r
      for(size_t i=num_iteration+1; i<n; ++i) {
        //compute beta such that B[num_iteration][i]== beta * b / p^r with beta > 1
        // G_.p_normalize(mat_[num_iteration][i]);
        Integer beta = division(top_left.second, mat_[num_iteration][i].second);
                                                              //in number_theory.h
        //do row_i <- row_i - beta * eps^-1 * row_num_iteration and
        //   col_i <- col_i - beta * eps^-1 * col_num_iteration
        auto z = kumquat::times(kumquat::times( beta, eps_inv ), (Integer(-1)));
        plus_equal_row(i,num_iteration,z);
        plus_equal_column(i,num_iteration,z);
      }
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
      os << std::setw(7) << mat[i][j] ;// << " ";
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