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

#ifndef KUMQUAT_DENSE_MATRIX_H_ 
#define KUMQUAT_DENSE_MATRIX_H_

#include <iomanip>
#include <vector>
#include <map>
#include <kumquat/number_theory.h>
#include <kumquat/Q_U1.h>

namespace kumquat {

/** \class Dense_matrix Dense_matrix.h kumquat/Dense_matrix.h 
 * \brief A dense matrix type, suited for normalizing small matrices.
 * 
 * The template type represents the group of coefficients for the matrix entries, 
 * and must be a model of the concept AbelianGroup. 
 * 
 * All the values of a same row are represented contiguously in memory ; hence 
 * row operations are more efficient than column operations.
 */
template< class CoefficientStructure >
class Dense_matrix : public ScalarSetOperations {
public:
/** The algebraic structure containing the coefficients for the matrix entries.
 * 
 * Must be a model of AbelianGroup.*/
  typedef CoefficientStructure Coeff_struct;
/** \brief The type of coefficients for the matrix entries.*/
  typedef typename Coeff_struct::Element Coefficient;
/** \brief Signed integer type. Must be a model of SignedInteger.*/
  typedef typename Coeff_struct::Integer Integer;


/** \name Model of ScalarSetOperations, and additional constructor.
 * 
 * @{ \*

/** \brief Creates an n by m matrix with uninitialized coefficients. */
  Dense_matrix(size_t n, size_t m, CoefficientStructure G) 
  : n_(n), m_(m), G_(G), dim_ker_(-1), dim_im_(-1) {
    mat_.resize(n);
    for(size_t i=0; i<n; ++i) { mat_[i].resize(m); } 
    row_idx_.reserve(n);
    for(size_t i=0; i<n; ++i) { row_idx_.push_back(i); }
    col_idx_.reserve(m);
    for(size_t j=0; j<m; ++j) { col_idx_.push_back(j); }
  }
/** \brief Copy constructor. */
  Dense_matrix(const Dense_matrix& other) 
  : n_(other.n_), m_(other.m_), mat_(other.mat_), G_(other.G_), dim_ker_(other.dim_ker_), dim_im_(other.dim_im_), row_idx_(other.row_idx_), col_idx_(other.col_idx_) {}
/** \brief Move constructor.*/
  Dense_matrix(Dense_matrix&& other) noexcept {
    move_from(other);
  } 
/** Destructor. */
  ~Dense_matrix() { clear_matrix(); }
/** \brief Copy assignment. */  
  Dense_matrix& operator=(const Dense_matrix& other) {
    copy_from(other);
    return *this;
  }
/** \brief Move assignment relocates the whole matrix. */
  Dense_matrix& operator=(Dense_matrix&& other) noexcept {
    if(&other != this) {
      clear_matrix();
      move_from(other);
    }
    return *this;
  }
/** \brief Test for equality.*/
  inline bool operator==(const Dense_matrix& lhs, const Dense_matrix& rhs) {
    return (   (lhs.mat_ == rhs.mat_) 
            && (lhs.row_idx_ == rhs.row_idx_) 
            && (lhs.col_idx_ && rhs.col_idx_) );
  }
/** \brief Test for inequality.*/
  inline bool operator!=(const Dense_matrix& lhs, const Dense_matrix& rhs) {
    return !(lhs == rhs);
  }
private:
//copy other into this
  void copy_from(Dense_matrix& other) {
    if(this != &other) {
      n_ = other.n_;
      m_ = other.m_;
      mat_ = other.mat_;
      G_ = other.G_;
      dim_ker_ = other.dim_ker_;
      dim_im_ = other.dim_im_;
      row_idx_ = other.row_idx_;
      col_idx_ = other.col_idx_;
      // n_ = other.num_rows();
      // m_ = other.num_columns();
      // //copy the mat_
      // mat_.clear();   mat_.resize(n_);
      // for(size_t i=0; i < n_; ++i) {
      //   mat_[i].resize(m_);
      //   for(size_t j=0; j < m_; ++j) {
      //     mat_[i][j] = other[i][j];//CoefficientStructure::Element must be assignable
      //   }
      // }
      // //copy the coefficient structure
      // G_ = other.G_;
      // //copy remaining fields
      // dim_ker_ = other.dim_ker_; 
      // dim_im_ = other.dim_im_;
      // row_idx_.resize(n_); 
      // for(size_t i=0; i < n_; ++i) { row_idx_[i] = other.row_idx_[i]; }
      // col_idx_.resize(m_);
      // for(size_t j=0; j < m_; ++j) { col_idx_[j] = other.col_idx_[j]; }
    }
  }
//move from other to this.
  void move_from(Dense_matrix &other) {
    n_ = std::move(other.n_);
    m_ = std::move(other.m_);
    mat_ = std::move(other.mat_);
    G_ = std::move(other.G_);
    dim_ker_ = std::move(other.dim_ker_);
    dim_im_ = std::move(other.dim_im_);
    row_idx_ = std::move(other.row_idx_);
    col_idx_ = std::move(other.col_idx_);
  }
//clear the memory for the matrix
  void clear_matrix() {
    mat_.swap(std::vector< std::vector< Coefficient > >());
    row_idx_.swap(std::vector<size_t>());
    col_idx_.swap(std::vector<size_t>());
  }

/* @} */ // end ScalarSetOperations

public:

/** \name Accessors and simple calculations.
 * 
 * @{ */
/** \brief Access the element in row i and column j in the matrix. 
 * 
 * Undefined behavior if i or j is outside the range.
 * */
  Coefficient& operator()(size_t i, size_t j) {
    return mat_[row_idx_[i]][col_idx_[j]];
  }
/** \brief Return the total number of rows in the matrix.*/
  size_t num_rows() const { return n_; }
/** \brief Return the total number of columns in the matrix.*/
  size_t num_columns() const { return m_; }
/** \brief Return a reference to the algebraic structure of the entries of 
 * the matrix.*/
  Coeff_struct & coefficient_structure() { return G_; }
/** \brief Set the matrix to 0, while maintaining its size. **/
  void set_to_zero() {
    for(size_t i = 0; i < num_rows(); ++i) {
      for(size_t j = 0; j < num_columns(); ++j) {
        (*this)(i,j) = G_.additive_identity();
      }
    }
  }
/** \brief Access the vector encoding the row at index idx.*/
  std::vector< Coefficient > & operator[](size_t idx) {
    return mat_[row_idx_[idx]];
  };
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
/**\brief Return the trace of a square matrix.*/
  Coefficient trace() {
    if(num_rows() != num_columns()) {
      std::cout << "Error trace on non-square matrix.\n";
      return 0;
    }
    Coefficient tr = G_.additive_identity(0);
    for(size_t i=0; i<num_rows(); ++i) {
      G_.plus_equal(tr,(*this)(i,i));
    }
    return tr;
  }

/* @} */ // end accessors

/** \name Global modifications of the matrix.
 * @{ */

/** \brief Set the matrix to the identity. 
 * 
 * Works also if the matrix is rectangular, and turn the matrix into a block diagonal matrix, with a maximal top-left block identity, and all other coefficients to 0.
 * */
  void set_to_identity() {
    for(size_t i = 0; i < num_rows(); ++i) {
      for(size_t j = 0; j < num_columns(); ++j) {
        if(i==j) { (*this)(i,j) = G_.multiplicative_identity(); }
        else { 
          (*this)(i,j) = G_.additive_identity();
        }
      }
    }
  }
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
        (*this)(i,j) = f(i,j);
      }
    }
  }

/* @} */ // end global modifications of the matrix

/** \name Row and column operations.
 * 
 * @{ */
/** \brief Set \f$\col_i \leftarrow z \times \col_i\f$. 
 * 
 * If CoefficientStructure implements AbelianGroup, WeightType can be of type CoefficientStructure::Integer.
 * 
 * If CoefficientStructure also implements PseudoRing, WeightType can also be of type CoefficientStructure::Element (i.e., Coefficient).
 * */
  template<typename WeightType>
  void times_equal_column(size_t i, WeightType z) {
    for(size_t k=0; k<n_; ++k) {
      G_.times_equal( (*this)(k,i) , z );
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
      G_.times_equal( (*this)(i,k),z);
    }
  }
/** \brief Set col_i <- col_i + col_j. */
  void plus_equal_column(size_t i, size_t j) {
    for(size_t k=0; k<n_; ++k) {
      G_.plus_equal((*this)(k,i), (*this)(k,j));
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
      G_.plus_equal( (*this)(k,i), G_.times((*this)(k,j),z) );
    }    
  }
/** \brief Set row_i <- row_i + row_j. */
  void plus_equal_row(size_t i, size_t j) {
    for(size_t k=0; k<m_; ++k) {
      G_.plus_equal((*this)(i,k),(*this)(j,k));
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
      G_.plus_equal((*this)(i,k), G_.times((*this)(j,k), z) );
    }    
  }
/** \brief Exchange the columns of index i and j lazily.*/
  void exchange_column(size_t i, size_t j) {
    if(i==j) { return; }
    std::swap(col_idx_[i],col_idx_[j]);
    // for(size_t k=0; k<n_; ++k) {
    //   std::swap(mat_[k][i],mat_[k][j]);
    // }

  }
/** \brief Exchange the rows of index i and j lazily.*/
  void exchange_row(size_t i, size_t j) {
    if(i==j) { return; }
    std::swap(row_idx_[i],row_idx_[j]);
  }

/** \brief Labels for row and column operations.*/
  enum Operation_type {// ? means undefined
    plus_equal_column_t,//(i,j,z)  col_i <- col_i + z*col_j 
    plus_equal_row_t,//(i,j,z)     row_i <- row_i + z*row_j
    exchange_column_t,//(i,j,?)    col_i <-> col_j 
    exchange_row_t,//(i,j,?)       row_i <-> row_j 
    times_equal_column_t,//(i,?,z) col_i <- z*col_i 
    times_equal_row_t };//(i,?,z)  row_i <- z*row_i
/** \brief Type of row and column operations, to keep track of a reduction algorithm.*/
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

/* @} */ // end row column operations

/** \name Matrix reduction algorithms.
 * 
 * @{ */
/** \brief Reduce the matrix to column echelon form.
 * 
 * Keep track of the sequence of operation as a vector of 
 * tuple (i,j,z), meaning operation col_i <- col_i + z * col_j, 
 * pushed into the input col_ops.
 * 
 * The Dense_matrix is in column echelon form after the operation.
 * 
 * The coefficient structure must be a model of concept Field.
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
                          G_.times((*this)(curr_low,i),
                            G_.multiplicative_inverse(
                              (*this)(curr_low,it_conflict->second)
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
 * The coefficient structure must be a model of concept Field.
 * */
  void row_echelon_form(std::vector< 
                              std::tuple<size_t,size_t,Coefficient> > & row_ops) {
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
                          G_.times((*this)(i,curr_low),
                             G_.multiplicative_inverse(
                               (*this)(it_conflict->second,curr_low)
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
      if(!G_.trivial((*this)(j,i))) { return j; }
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
      if(!G_.trivial((*this)(i,j))) { return j; }
    }
    return -1;
  }

public:
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
  //call x=(*this)(r_idx,c_idx) and y=(*this)(k,c_idx)
        if(!G_.trivial((*this)(k,c_idx))) {//if y!=0
  //write u*x + v*y = gcd(x,y), with gcd > 0
          auto u_v_gcd = G_.extended_gcd((*this)(r_idx,c_idx), 
                                         (*this)(k,c_idx));
  //if x does not divide y, i.e., |gcd| < |x|
          if( std::get<2>(u_v_gcd) < G_.abs((*this)(r_idx,c_idx))) {
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
  //call x=(*this)(r_idx,c_idx) and y=(*this)(r_idx,k)
        if(!G_.trivial((*this)(r_idx,k))) {//if y!=0
  //write u*x + v*y = gcd(x,y), with gcd > 0
          auto u_v_gcd = G_.extended_gcd((*this)(r_idx,c_idx), 
                                         (*this)(r_idx,k)   );
  //if x does not divide y, i.e., |gcd| < |x|
          if( std::get<2>(u_v_gcd) < G_.abs((*this)(r_idx,c_idx)) ) {
  //perform col_{c_idx} <- u * col_{c_idx} + v * col_k
            times_equal_column(c_idx, std::get<0>(u_v_gcd));
            ops.emplace_back(c_idx,c_idx,std::get<0>(u_v_gcd),times_equal_column_t);
            plus_equal_column(c_idx, k, std::get<1>(u_v_gcd));
            ops.emplace_back(c_idx, k, std::get<1>(u_v_gcd),plus_equal_column_t);
          }//now x <- gcd(x,y) and new_x divides y
        }
      }
  //now (*this)(r_idx,c_idx) divides all entries in its column and row -> cancel the column then the row
      for(size_t k = r_idx+1; k<num_rows(); ++k) {
  //call x=(*this)(r_idx,c_idx) and y=(*this)(k,c_idx)
        if(!G_.trivial((*this)(k,c_idx))) {//if y!=0
  //perform row_k <- y/x * row_{r_idx}
          auto z = G_.times(G_.division((*this)(k,c_idx),(*this)(r_idx,c_idx)),-1);
          plus_equal_row(k, r_idx, z);
          ops.emplace_back(k, r_idx, z, plus_equal_row_t);
        }
      }
  //we can directly trivialize the row r_idx now
      for(size_t k = c_idx+1; k<num_columns(); ++k) {
        if(!G_.trivial((*this)(r_idx,k))) {//record the operation
          auto z = G_.times(G_.division((*this)(r_idx,k),(*this)(r_idx,c_idx)),-1);
          ops.emplace_back(k, c_idx, z, plus_equal_column_t);
          (*this)(r_idx,k) = G_.additive_identity();//set directly to 0
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
      if(!G_.trivial((*this)(i,idx))) { return i; }
    }
    return -1;
  }

public:
/** \brief Diagonalize the Gram matrix of a bilinear form.
 * 
 * The matrix mat_ must contain the values b(x_i,x_j) of a bilinear form from a 
 * finite abelian p-group G = <x_1, ... x_n> to the additive abelian group 
 * \f$\mathbb{Q}_{(p)}/\mathbb{Z}\f$, for a prime number \f$p \geq 2\f$. In 
 * particular, the matrix is square symmetric.
 * 
 * The result is:
 * for p odd, a diagonal matrix with diagonal elements of the form \f$a/p^r\f$ with 
 * \$f\gcd(a,p)=1\$f and \f$a\f$ is either 1 or a quadratic non-residue modulo p.
 * 
 * for p = 2, a block diagonal matrix with blocks of size 1 or 2. The 1-block are 
 * of the form 
 * \f$a/2^r\f$ with \$f\gcd(a,2)=1\$f and \f$a\f$ is either 1 or a quadratic 
 * non-residue modulo 2^r. The 2-blocks are of the form:
 * 
 * | 0      1/2^r |          | 1/2^(r-1)   1/2^r     |
 * |                   or    |                       |
 * | 1/2^r  0     |          | 1/2^r       1/2^(r-1) |
 * */
  void diagonalize_gram_matrix_Qp_mod_Z() {
    auto n = num_rows();
    auto m = num_columns();
    if(n != m) { 
      std::cerr << "Gram matrix is not square.\n"; return; 
    }
    auto p = G_.p();
    //distinguish odd and even prime p
    if((p % 2) == 0) { diagonalize_gram_matrix_Qp_mod_Z_p_even(); }
    else { diagonalize_gram_matrix_Qp_mod_Z_p_odd(); }
  }
/** \brief Compute the Gauss sum associated to a general bilinear map on a p-group.
 * */
  typename Q_U1<Integer>::Element gauss_sum_Qp_mod_Z() {
    Q_U1<Integer> q_u1;
    diagonalize_gram_matrix_Qp_mod_Z();
    auto gauss_sum = q_u1.additive_identity();//1 Q_U! is a multiplicative group

    auto n = num_rows();//square nxn matrix
    auto p = G_.p();
    if((p % 2) == 0) { 
      for(size_t i=0; i<n; ) {

        if( (i != n-1) && !(G_.trivial((*this)(i,i+1)))) {//2 x 2 block of the form:
        /*
         * | 0      1/2^r |          | 1/2^(r-1)   1/2^r     |
         * |              |    or    |                       |
         * | 1/2^r  0     |          | 1/2^r       1/2^(r-1) |
         */
          if(!G_.trivial((*this)(i,i))) {//type 2, Gauss sum is (-1)^r
            auto two_to_r = G_.denominator((*this)(i,i+1));//2^r
            Integer r = logp(two_to_r,(Integer)2);
            if(r % 2 == 1) {
              q_u1.plus_equal(gauss_sum, q_u1.element(-1,1,0,1));
            }
          }
          // else type 1 and Gauss sum is 1, do nothing
          i+=2;
        }
        else { //1 x 1 block of value a/2^r
          auto a = G_.numerator((*this)(i,i));
          auto res_a = (a % 8);
          auto res_a2 = res_a * res_a;// (a mod 8)^2
     
          if(res_a2 == 1) {//gauss_sum *= exp(i 2 pi a/8)
            q_u1.plus_equal( gauss_sum, q_u1.element(1,1,res_a,8) );
          }
          else {//res_a2 == 25, gauss_sum *= (-1)^{r*(res_a2-1)/8} exp(i 2 pi a/8)
                //                        *= (-1)^r exp(i 2 pi a/8)  
            auto two_to_r = G_.denominator((*this)(i,i));//2^r
            Integer r = logp(two_to_r,(Integer)2);//r
            if(r % 2 == 1) {
              q_u1.plus_equal(gauss_sum, q_u1.element(-1,1,res_a,8));
            }
            else {
              q_u1.plus_equal(gauss_sum, q_u1.element(1,1,res_a,8));
            }
          }
          ++i;
        }
      }
    }
    else {//p odd, the matrix is diagonal
      for(size_t i=0; i<n; ++i) {
        auto p = G_.p();
        auto a = G_.numerator((*this)(i,i));
        auto p_to_r = G_.denominator((*this)(i,i));//p^r

        auto legendre_2a_p = legendre_symbol((Integer)(2)*a, p);//in 0,1,-1
        if(legendre_2a_p == -1) {
          auto r = logp(p_to_r,p);
          if(r%2 == 0) { legendre_2a_p = 1; }
        }
        //now legendre_2a_p == (legendre(2a,p))^r
        if(p_to_r % 4 == 1) {
          q_u1.plus_equal(gauss_sum, q_u1.element(legendre_2a_p,1,0,1));
        }
        else {
          q_u1.plus_equal(gauss_sum, q_u1.element(legendre_2a_p,1,1,4));
        }
      }
    }
    return gauss_sum;
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
    for(size_t num_iteration = 0; num_iteration < n; ) {
      //check if the bottom right block B is uniformly 0, and if not compute the minimum r such that p^r B = 0;
      int pivot_i, pivot_j;
      find_pivot_Qp_mod_Z(pivot_i, pivot_j, num_iteration);

      if(pivot_i == -1) { return; }//the remaining matrix is trivial
      
      //now, the matrix is non trivial, element B[i][j] has minimal order, and if there is a diagonal coefficient of minimal order, then i==j
      if(pivot_i == pivot_j) {//put on top left corner
        exchange_row(num_iteration,pivot_i);
        exchange_column(num_iteration,pivot_i);
        //use the pivot to cancel row[num_iteration...m] and col[num_iteration...n] 
        cancel_row_column_Qp_mod_Z(num_iteration);
        ++num_iteration;
      }
      else {//pivot_i != pivot_j and B[i][j]==B[j][i] has strictly larger order than B[i][i] and B[j][j]
        //put the future block in top left corner
        exchange_row(num_iteration,pivot_i);
        exchange_column(num_iteration,pivot_i);
        exchange_row(num_iteration+1,pivot_j);
        exchange_column(num_iteration+1,pivot_j);

        cancel_2_2_block_Q2_mod_Z(num_iteration);
        num_iteration += 2;
      }
    }
  }
/* Find the element in the sub-matrix M[idx...n][idx...n] with maximal order in Qp_mod_Z.
If such element is in the diagonal, always return a diagonal element. If the matrix is uniformly 0, the pivot indices are returned as -1.
*/
  void find_pivot_Qp_mod_Z(int & pivot_i, int & pivot_j, int idx=0) {
    Integer max_order = 0;
    pivot_i = -1; pivot_j = -1;//find element of minimal order
    for(size_t i = idx; i < num_rows(); ++i) {//try to find a pivot in the diagonal
      for(size_t j = idx; j < num_columns(); ++j) {
        if(!G_.trivial((*this)(i,j))) {//!= 0
          auto ord = G_.order((*this)(i,j));
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
/* Cancel the row and column of given index idx assuming that the element M[idx][idx] has maximal order in Qp_mod_Z over all other elements in M[idx...n][idx...n]. Subroutine of diagonalization of Gram matrices.*/ 
  void cancel_row_column_Qp_mod_Z(size_t idx) {
        auto p = G_.p(); auto n = num_rows();
    //top left element M[idx][idx] has minimal order
    // and is written as top_left == a p^{-r} with gcd(a,p)==1 and 0 < a < p
    Coefficient top_left = (*this)(idx,idx);
    //now, enforce top left corner B[idx][idx]== eps/p^k, for 
    //eps = 1 or eps quadratic non-residue mod p. Compute also eps^{-1} mod p^r
    Integer eps, eps_inv;
    //if a quadratic residue mod p, i.e., a = x^2 mod p -> set eps = 1
    if(quadratic_residue(top_left.first,top_left.second,p)) {//compute a solution x
      auto x = solve_quadratic_residue(top_left.first,top_left.second,p);
      //compute s = x^{-1} mod p^k, which exists because gcd(a,p)=1 => gcd(x,p)=1
      auto s = inverse(x,top_left.second);
      //set row_idx <- s*row_idx and col_idx <- s*col_idx such that B[idx][idx] = 1/p^r mod Z
      times_equal_row(idx,s);
      times_equal_column(idx,s);
      eps = 1; eps_inv = 1;
    }
    else {//else quadratic non-residue -> top_left is already a/p^k, set eps=a
      eps = top_left.first;
      eps_inv = kumquat::inverse(eps, top_left.second);//in number_theory.h
    }

    top_left = (*this)(idx,idx);
  //we do have M[idx][idx]== eps/p^k, for eps = 1 or 
  //eps quadratic non-residue mod p, and eps_inv = eps^-1 mod p^r
    for(size_t i=idx+1; i<n; ++i) {
  //compute beta such that B[idx][i]== beta * b / p^r with beta > 1
  // G_.p_normalize((*this)(idx,i));
      Integer beta = kumquat::division(top_left.second, (*this)(idx,i).second);//in number_theory.h
      //do row_i <- row_i - beta * b * eps^-1 * row_idx and
      //   col_i <- col_i - beta * eps^-1 * col_idx
      auto z = kumquat::times(
                 kumquat::times( 
                   kumquat::times( beta, 
                                   (*this)(idx,i).first ) 
                   , eps_inv )
                 , (Integer(-1)));

      plus_equal_row(i,idx,z);
      plus_equal_column(i,idx,z);
    }
  } 
/* Subroutine of the diagonalization of the Gram matrix with coefficients in Q_2 / Z, where the block B[idx,idx+1][idx,idx+1] is of the form:
*  |2a/2^m   b/2^m |   where b is odd and a,c are arbitrary
*  |               |
*  |b/2^m    2c/2^m|
*
*  in row/col i>idx+1, we have the first elements equal to 
*  |u/2^m    v/2^m | for cancellation.
*
*  as a result, we get a canonical block of the form:
*  | 0      1/2^r |          | 1/2^(r-1)   1/2^r     |
*  |                   or    |                       |
*  | 1/2^r  0     |          | 1/2^r       1/2^(r-1) |
*/
  void cancel_2_2_block_Q2_mod_Z(size_t idx) {
    //extract integers a,b,c and 2^m
    Integer two_to_m = (*this)(idx,idx+1).second;
    Integer a = (*this)(idx,idx).first * (two_to_m / (2*(*this)(idx,idx).second)) ;
    Integer b = (*this)(idx,idx+1).first;
    Integer c = (*this)(idx+1,idx+1).first * 
                              (two_to_m / ((Integer)(2)*(*this)(idx+1,idx+1).second));
    //d is the inverse of 4ac-b^2 modulo 2^m
    Integer d_inv = ( (Integer)4 * a * c - b * b ) % two_to_m;
    if(d_inv < 0) { d_inv += two_to_m; }
    Integer d = kumquat::inverse(d_inv, two_to_m);

    //for all i >idx+1, cancel M[i][idx,idx+1] and M[idx,idx+1][i]
    for(size_t i=idx+2; i<num_rows(); ++i) {
      // find u s.t. (*this)(i,idx) == u / 2^m and v s.t. (*this)(i,idx+1) == v / 2^m
      Integer u = ((*this)(i,idx)).first * (two_to_m / (*this)(i,idx).second); 
      Integer v = ((*this)(i,idx+1)).first * (two_to_m / (*this)(i,idx+1).second);

      Integer minus_r1 = (-1) * d * ( 2 * c * u - b * v);
      Integer minus_r2 = (-1) * d * ( 2 * a * v - b * u);

      plus_equal_row(i,idx,minus_r1);
      plus_equal_column(i,idx,minus_r1);
      plus_equal_row(i,idx+1,minus_r2);
      plus_equal_column(i,idx+1,minus_r2);
    }
    /* now, put the 2x2 block
     *  |2a/2^m   b/2^m |   where b is odd and a,c are arbitrary
     *  |               |
     *  |b/2^m    2c/2^m|
     * into canonical form.*/
    //lemma 2.2 in https://arxiv.org/pdf/1405.7950.pdf
    //if ac even then equivalent to:
    /*  |0     1/2^m |   
     *  |            |
     *  |1/2^m      0|
     */
    //if ac odd, then equivalent to:
    /*  |1/2^m-1  1/2^m  |   
     *  |                |
     *  |1/2^m    1/2^m-1|
     */

    if((a*c) % 2 == 0) {
      (*this)(idx,idx) = G_.additive_identity();
      (*this)(idx+1,idx) = G_.element(1,two_to_m);
      (*this)(idx,idx+1) = G_.element(1,two_to_m);
      (*this)(idx+1,idx+1) = G_.additive_identity();
    }
    else {
      (*this)(idx,idx) = G_.element(2,two_to_m);
      (*this)(idx+1,idx) = G_.element(1,two_to_m);
      (*this)(idx,idx+1) = G_.element(1,two_to_m);
      (*this)(idx+1,idx+1) = G_.element(2,two_to_m);
    }
  }

/* @} */ //end matrix reductions

public:


public:  
/** \name Model of ScalarRingOperations, and additional matrix specific multiplication operations.
 * @{ */

/** \brief Matrix multiplication on the right. 
 * 
 * Naive cubic algorithm.
 **/
  Dense_matrix rtimes(const Dense_matrix& rhs) {
    assert( num_columns() == rhs.num_rows() );
    Dense_matrix prod_mat(num_rows(),rhs.num_columns(),G_); 
    for(size_t i = 0; i < num_rows(); ++i) {
      for(size_t j = 0; j < rhs.num_columns(); ++j) {
        prod_mat(i,j) = G_.additive_identity();
        for(size_t k = 0 ; k < num_columns() ; ++k) {
          G_.plus_equal( prod_mat(i,j), 
                         G_.times( (*this)(i,k) , rhs(k,j) ) ) ;  
        }
      }
    }
    return prod_mat;
  }
/** \brief Matrix multiplication on the right.
 * 
 * Set \mathtt{(*this) <- (*this) * rhs}, based on \f$\mathtt{rtimes}\f$.*/
  void rtimes_equal(const Dense_matrix& rhs) {
    *this = rtimes(rhs);
  }
/** \brief Matrix multiplication on the left. 
 * 
 * Naive cubic algorithm.
 * */
  Dense_matrix ltimes(const Dense_matrix& lhs) {
    assert( num_rows() == lhs.num_columns() );
    Dense_matrix prod_mat(lhs.num_rows(),num_columns(),G_); 
    for(size_t i = 0; i < lhs.num_rows(); ++i) {
      for(size_t j = 0; j < num_columns(); ++j) {
        prod_mat(i,j) = G_.additive_identity();
        for(size_t k = 0 ; k < lhs.num_columns() ; ++k) {
          G_.plus_equal( prod_mat(i,j), 
                         G_.times( lhs(i,k) , (*this)(k,j) ) ) ;  
        }
      }
    }
    return prod_mat; //move
  }
/** \brief Matrix multiplication on the left.
 * 
 * Set \f$\mathtt{(*this) <- lhs * (*this)}, based on \f$\mathtt{ltimes}.*/
  void ltimes_equal(const Dense_matrix& lhs) {
    *this = ltimes(rhs);
  }
/** \brief Matrix multiplication on the right by a scalar. 
 * 
 * The scalar may be an integer type, in which case coefficient-wise multiplication are the ones from the Z-module structure of group coefficients, or may be of the same type as the matrix coefficients, in which case the matrix coefficients must belong to a ring. 
 * */
  template<typename Scalar>
  void rtimes_equal(Scalar x) {
    for(size_t i = 0; i < num_rows(); ++i) {
      for(size_t j = 0; j < rhs.num_columns(); ++j) {
          G_.times_equal( (*this)(i,j), x ) ;  
      }
    }
  }

/** \brief Multiplication on the right by a scalar.
 * 
 * Based on \f$\mathtt{rtimes_equal}\f$. */
  template<typename Scalar>
  Dense_matrix& operator*=(Scalar x) {
    this->rtimes_equal(x);
    return *this;
  }
/** \brief Multiplication on the right by a matrix. 
 * 
 * Based on \f$\mathtt{rtimes_equal}\f$. */
  Dense_matrix& operator*=(const Dense_matrix& rhs) {
    this->rtimes_equal(rhs);
    return *this;
  }
/* } */ // end ScalarRingOperations 






  /** \brief Compute the tensor product on the right with the input matrix.
   * 
   * this <- this tensor rhs.
   **/
  Dense_matrix rtensor(const Dense_matrix &rhs) {
    Dense_matrix res( num_rows()*rhs.num_rows()
                    , num_columns()*rhs.num_columns(), G_);
    for(size_t i=0; i<num_rows(); ++i) {
      for(size_t j=0; j<num_columns(); ++j) {
        for(size_t k=0; k<rhs.num_rows(); ++k) {
          for(size_t l=0; l<rhs.num_columns(); ++l) {
            res(i*rhs.num_rows() + k, j*rhs.num_columns()+l) = 
                                  G_.times((*this)(i,j), rhs(k,l));
          }
        }
      }
    }
    return res;
  }
/** \brief Tensoring to the right.
 * 
 * Set *this <- *this \f$\otimes\f$ rhs, based on rtensor.*/
  void rtensor_equal(const Dense_matrix& rhs) {
    *this = rtensor(rhs);
  }

  /** \brief Compute the tensor product on the left with the input matrix.
   * 
   * this <- lhs tensor this.
   **/
  Dense_matrix ltensor(const Dense_matrix& lhs) {
    Dense_matrix res( num_rows()*lhs.num_rows()
                    , num_columns()*lhs.num_columns(), G_);
    for(size_t i=0; i<lhs.num_rows(); ++i) {
      for(size_t j=0; j<lhs.num_columns(); ++j) {
        for(size_t k=0; k<num_rows(); ++k) {
          for(size_t l=0; l<num_columns(); ++l) {
            res(i*num_rows() + k, j*num_columns()+l) = 
                                G_.times(lhs(i,j), (*this)(k,l));
          }
        }
      }
    } 
    return res;
  }

/** \brief Return the tensor product \f$\operatorname{id}_m \otimes M \otimes \operatorname{id}_n\f$,
   * where \f$\operatorname{id}_k\f$ is the k by k identity matrix.
  **/
  Dense_matrix tensor_id(int m, int n) {
    return (this->rtensor(identity_matrix(n))).ltensor(identity_matrix(m));
  }

/** \brief Return the n by n identity matrix.*/
  Dense_matrix identity_matrix(size_t n) {
    Dense_matrix id(n,n, G_);
    for(size_t i=0; i<n; ++i) { id(i,i) = G_.multiplicative_identity(); }
    return id;
  }

private:
  //number of rows of the matrix
  size_t n_;
  //number of columns of the matrix
  size_t m_;
  //encoding of the matrix as a bi-dimensional array. The element at row i and col j is in (*this)(i,j)
  std::vector< std::vector< Coefficient > > mat_;
  //the group to which the coefficients of the matrix belong.
  CoefficientStructure G_;
  //dim of the kernel, -1 if not initialized
  int dim_ker_;
  //dim of the image, -1 if not initialized
  int dim_im_;
  //encodes a permutation of row indices, such that the row of index i is in position mat_[row_idx_[i]] in the data structure
  std::vector<size_t> row_idx_;
  //encodes a permutation of col indices, such that the col of index i is in position mat_[...][col_idx_[i]] in the data structure
  std::vector<size_t> col_idx_;

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