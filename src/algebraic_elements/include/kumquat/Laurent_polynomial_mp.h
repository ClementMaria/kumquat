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

#ifndef KUMQUAT_LAURENT_POLYNOMIAL_MP_H_ 
#define KUMQUAT_LAURENT_POLYNOMIAL_MP_H_

#include <vector>
#include <numeric>
#include <kumquat/number_theory.h>
#include <boost/multiprecision/gmp.hpp>//boost multiprecision wrap over gmp
// #include <boost/integer/extended_euclidean.hpp>//boost extended gcd
// #include <boost/math/tools/polynomial.hpp>
#include <list>
#include <map>

namespace kumquat {

/** \brief A data type for Laurent polynomial with multi-precision integer coefficients.
 * 
 * Is model of concept ScalarRingOperations
 */
class Laurent_polynomial_mp {
public:
  
/** \brief An integer type for the coefficients.*/
  typedef boost::multiprecision::mpz_int Integer;//type of coefficients

/** \brief An integer type for the coefficients.*/
  typedef boost::multiprecision::mpz_int Coefficient;//type of coefficients

/** \brief Default initialization to 0.*/
  Laurent_polynomial_mp() : poly_() {}
/** \brief Initialization with an integer z, z/1.*/
  // template<typename IntegerType>
  Laurent_polynomial_mp(int z) {
    if(z!=0) { poly_.emplace_back(0,z); }
  } 
/** \brief Initialization with a range.*/
  // template<typename RangeType>
  // Laurent_polynomial_mp(RangeType &poly) 
  // {
  //   poly_ = std::list< std::pair<int, Coefficient> >(poly.begin(),poly.end());
  // }
/** \brief Initialize as a monomial.*/
  template<typename IntegerType>
  Laurent_polynomial_mp(int deg, IntegerType a) {
    if(a != 0) { poly_.emplace_back(deg,(Coefficient)(a)); }
  }
/** \name Model of ScalarSetOperations
 * 
 * @{ */
/** \brief Copy constructor.*/
  Laurent_polynomial_mp(const Laurent_polynomial_mp& other) {
    poly_ = other.poly_;
  }
/** \brief Move constructor.*/
  Laurent_polynomial_mp(Laurent_polynomial_mp&& other) noexcept {
    poly_ = std::move(other.poly_);
  }
/** brief Destructor.*/  
  ~Laurent_polynomial_mp() {}
/** \brief Copy assignment.*/
  Laurent_polynomial_mp& operator=(const Laurent_polynomial_mp& other) 
  {
    poly_ = other.poly_;
    return *this;
  }
/** \brief Move assignment.*/
  Laurent_polynomial_mp& operator=(Laurent_polynomial_mp&& other) noexcept
  {
    poly_ = other.poly_;
    return *this;
  }

  inline bool operator==(const Laurent_polynomial_mp& rhs) const { 
    return (poly_ == rhs.poly_);
   }
  inline bool operator!=(const Laurent_polynomial_mp& rhs) const { 
    return !((*this) == rhs); 
  }

  inline bool operator==(const int rhs) const { 
    Laurent_polynomial_mp rf_rhs(0,rhs);
    return (*this) == rf_rhs;
   }
  inline bool operator!=(const int rhs) const { 
    return !((*this) == rhs); 
  }
/* @} */

/** \name Model of ScalarGroupOperations
 * 
 * @{ */
  Laurent_polynomial_mp& operator+=(const Laurent_polynomial_mp& rhs) 
  { 
    auto itl = poly_.begin(); auto itr = rhs.poly_.begin();
    while(itl != poly_.end() && itr != rhs.poly_.end()) {
      if(itl->first == itr->first) {//same degree, sum coefficients
        itl->second += itr->second;
        if(itl->second == 0) {
          auto tmp = itl;
          ++itl; ++itr;
          poly_.erase(tmp);
        }
        else { ++itl; ++itr; }
      }
      else {
        if(itl->first < itr->first) { ++itl; }
        else { poly_.insert(itl, *itr); ++itr; }
      }
    }

    while(itr != rhs.poly_.end()) {
      poly_.push_back(*itr); ++itr;
    }
    return *this;
  }

  friend Laurent_polynomial_mp operator+(Laurent_polynomial_mp lhs,        
                     const Laurent_polynomial_mp& rhs) 
  {
    lhs += rhs;
    return lhs;
  }


  Laurent_polynomial_mp& operator-=(const Laurent_polynomial_mp& rhs) 
  { 
    auto itl = poly_.begin(); auto itr = rhs.poly_.begin();
    while(itl != poly_.end() && itr != rhs.poly_.end()) {
      if(itl->first == itr->first) {//same degree, sum coefficients
        itl->second -= itr->second;
        if(itl->second == 0) {
          auto tmp = itl;
          ++itl; ++itr;
          poly_.erase(tmp);
        }
        else { ++itl; ++itr; }
      }
      else {
        if(itl->first < itr->first) { ++itl; }
        else { poly_.emplace(itl, itr->first, (-1)*itr->second); ++itr; }
      }
    }

    while(itr != rhs.poly_.end()) {
      poly_.emplace_back(itr->first, (-1)*itr->second); ++itr;
    }

    return *this;
  }

  friend Laurent_polynomial_mp operator-(Laurent_polynomial_mp lhs, const Laurent_polynomial_mp& rhs) 
  {
    lhs -= rhs;
    return lhs;
  }

  /**
   * @brief      Multiplication by an integer (Z-module structure).
   *
   * @param[in]  rhs   The right hand side
   *
   * @return     The result of the multiplication assignment
   */
  Laurent_polynomial_mp& operator*=(Coefficient rhs)
  { 
    if(rhs == 0) { poly_.clear(); }
    else {
      for(auto it=poly_.begin(); it!= poly_.end(); ++it) {
        it->second *= rhs;
      }
    }
    return *this; 
  }
/* @} */

/** \name Model of ScalarRingOperations
 * 
 * @{ */
  Laurent_polynomial_mp& operator*=(const Laurent_polynomial_mp& rhs)
  {                 
    std::map<int,Coefficient> prod;
    for(auto itl = poly_.begin(); itl != poly_.end(); ++itl) {
      for(auto itr = rhs.poly_.begin(); itr != rhs.poly_.end(); ++itr) {

        Coefficient c = itl->second * itr->second;
        auto res_ins = prod.emplace(itl->first + itr->first, c);
        if(!res_ins.second) {//key already there
          res_ins.first->second += c;
          if(res_ins.first->second == 0) {
            prod.erase(res_ins.first);
          }
        }

      }
    }
    poly_.clear();
    for(auto pp : prod) { poly_.push_back(pp); }
    return *this; 
  }

  friend Laurent_polynomial_mp operator*(Laurent_polynomial_mp lhs, const Laurent_polynomial_mp& rhs) 
  {
    lhs *= rhs;
    return lhs;
  }
/* @} */

/** \name Model of ScalarFieldOperations
 * 
 * @{ */
  // Laurent_polynomial_mp& operator/=(const Laurent_polynomial_mp& rhs) 
  // {                     
  //   numerator_ *= rhs.denominator_;
  //   denominator_ *= rhs.numerator_;
  //   normalize();             
  //   return *this;
  // }

  // friend Laurent_polynomial_mp operator/(Laurent_polynomial_mp lhs,  const Laurent_polynomial_mp& rhs) 
  // {
  //   lhs /= rhs;
  //   return lhs;
  // }
/* @} */

  // Polynomial polynomial() const { return poly_; }

  std::string to_string() const {
    std::stringstream ss;
    for(auto pp : poly_) {
 
      if(pp.first == 0) { ss << "(" << pp.second << ")"; }//x^0
      else {//deg != 0
        if(pp.second == 1) { ss << "(x^" << pp.first << ")"; }
        else {
          if(pp.second == -1) { ss << "(-x^" << pp.first << ")"; }
          else { ss << "(" << pp.second << "x^" << pp.first << ")"; }
        }
      }
    }
    if(poly_.empty()) { ss << "(0)"; }
    return ss.str();
  }


  int sq_norm() {
    Integer nrm = 0;
    for(auto pp : poly_) {
      nrm += pp.first * pp.first;//non zero deg squared//std::abs(pp.first);
      // if(pp.second < 0) { nrm -= pp.second; }
      // else { nrm += pp.second; }//plus coefficient
    }

    if((int)nrm < 0) { std::cout << "sq_norm < 0?\n";}

    return (int)nrm;
  }

  int norm() {
    Integer nrm = 0;
    for(auto pp : poly_) {
      nrm += std::abs(pp.first);// * pp.first;//non zero deg squared//std::abs(pp.first);
      // if(pp.second < 0) { nrm -= pp.second; }
      // else { nrm += pp.second; }//plus coefficient
    }

    if(nrm < 0) { std::cout << "sq_norm < 0?\n";}

    return (int)nrm;
  }

  int spread() {
    
    if(poly_.empty()) { std::cout << "empty poly_\n"; return 0; }
    int sprd = (poly_.rbegin()->first - poly_.begin()->first);
    if(sprd <= 0 ) { std::cout << "spread <= 0\n"; }
    return sprd;//highest degree - smallest degree
  }

private:
/**
 * The list of non-zero monomials, represented by pairs (degree, coefficient). The list is always sorted by increasing degree.
 */
  std::list< std::pair<int, Coefficient> > poly_;

};

  

std::ostream& operator<<(std::ostream& os, const Laurent_polynomial_mp& x)
{
  os << x.to_string();
  return os;
}


}  //namespace kumquat

#endif //KUMQUAT_LAURENT_POLYNOMIAL_MP_INTEGRAL_MP_H_