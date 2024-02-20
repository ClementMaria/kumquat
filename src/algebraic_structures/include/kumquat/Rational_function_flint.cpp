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

#ifndef RATIONAL_FUNCTION_FLINT_
#define RATIONAL_FUNCTION_FLINT_

#include "fmpz_poly_qxx.h"
#include <Eigen/SparseCore>

//Wrapper around flint::fmpz_poly_qxx to be used as matrix coefficients in eigen
//maintains scalar elements of \Z(x)
class Rational_function_flint {
public:
  typedef slong                Basic_scalar;
  typedef flint::fmpz_poly_qxx Fraction;

private:
  Fraction val_;
public:
  //default initialize at 0
  Rational_function_flint() { val_ = 0; }
  //copy constructor
  Rational_function_flint(
                          const Rational_function_flint& other) 
  { val_ = other(); }

  //monom of (d,a) -> a*X^d)
  Rational_function_flint(int d, slong a) 
  { 
    Fraction monom("2  0 1");//X
    // val_ = Fraction(a);
    if(d < 0) { 
      d = -1 * d;
      monom /= (monom*monom);//X^-1
    } //a X^{-|d|}
    pow(monom, (unsigned int)d);
    val_ = monom;
    val_ *= a;
  }  
 
 //  //initialize to a integer type, that fits in slong
 //  template<typename IntegerFractionType>
 //  Rational_function_flint(IntegerFractionType a) {
 //   val_.set_zero(); 
 //   val_ *= a; 
 // }
  //explicit list of monomials in numerator and denominator.
  //return the fraction lnum / lden. A pair (d, a) means a*X^d
  Rational_function_flint( 
                 std::vector< std::pair<unsigned int, Basic_scalar> > & lnum
               , std::vector< std::pair<unsigned int, Basic_scalar> > & lden )
  {
    Fraction z, num, den;  z = 0; num = 0; den = 0;
    Fraction x("2  0 1");//monom X
    //num
    if(lnum.empty()) { num = 1; }
    else {
      for(auto monom : lnum) {//(a,d) -> a*X^d
        Fraction ax; 
        ax   = monom.second * pow(x, monom.first);
        num += ax;
      }
    }
    //den
    if(lden.empty()) { den = 1; }
    else {
      for(auto monom : lden) {//(a,d) -> a*x^d
        Fraction ax; ax = monom.second * pow(x, monom.first);
        den += ax;
      }
    }
    val_ = num/den;
  }

bool is_zero() {return val_.is_zero();}

  // //construct a polynomial from interpolation
  // Rational_function_flint( 
  //                                std::vector< std::pair<...,...> interpol > )
  // {
  //   //...
  // }


//assign the value val to the member flint::fmpz_poly_qxx.
  void assign(Fraction val) { val_ = val; }

//evaluate numerator and denominator at number z (int, float, complex, etc)
  // template< NumericType >
  // std::pair< NumericType, NumericType > 
//exact symbolic evaluation at exp(q i pi / N)
  // void evaluate(unsigned int q, unsigned int N) {
    // auto num = val_.num(); auto den = val_.den();
    // auto x = num.get_coeff(0).evaluate();

    // int xx = num(0);

    // std::cout << x << "\n";
    // std::cout << 

    // for(int i=0 ; i<10; ++i) {
    //   std::cout << num.get_coeff(i) << "  ";
    // }
    // std::cout << std::endl;

    // flint::fmpzxx xx = num.get_coeff(0);

    // int a0 = num.coeff(2);
    // std::cout << "a0 = " << a0 << std::endl;

    //num.degree() den.degree()
    //(k, a) means a* exp(q i pi / N)^k, s.t. 0 \geq k*q / N < 2
    // std::map< unsigned int, int > eval_num;
    // for(unsigned int k=0; k<= num.degree(); ++k) {
    //   slong a = num.get_coeff(k);
    //   if(a != 0) {
    //     int exp_pow = (k*q);
    //     while( exp_pow > 2*N ) { exp_pow -= 2*N; }
    //     auto res_ins = eval_num.insert(exp_pow, a);
    //     if(!(res_ins.second)) { res_ins.first->second += a; }
    //   }
    // }
    // std::map< unsigned int, int > eval_den;
    // for(unsigned int k=0; k<= den.degree(); ++k) {
    //   slong a = den.get_coeff(k);
    //   if(a != 0) {
    //     int exp_pow = (k*q);
    //     while( exp_pow > 2*N ) { exp_pow -= 2*N; }
    //     auto res_ins = eval_den.insert(exp_pow, a);
    //     if(!(res_ins.second)) { res_ins.first->second += a; }
    //   }
    // }
    // std::cout << "Num: ";
    // for(auto pp : eval_num) { 
    //   if(pp.second!=0) {
    //     std::cout << pp.second << ".eipi^(" 
    //                                   << pp.first << "*" << q << "/" << N << ") + ";
    //   } 
    // }
    // std::cout << std::endl;
    // std::cout << "Den: ";
    // for(auto pp : eval_den) { 
    //   if(pp.second!=0) {
    //     std::cout << pp.second << ".eipi^(" 
    //                                   << pp.first << "*" << q << "/" << N << ") + ";
    //   } 
    // }

  // }

//returns a copy of the member flint::fmpz_poly_qxx.
  Fraction operator()() const { return val_; }
//copy assignment
  Rational_function_flint& operator=(
                                      const Rational_function_flint& other)
  {
    if (this != &other) { val_ = other(); }
    return *this;
  }
//copy assignment with a scalar of type slong
  Rational_function_flint& operator=(const slong& other)
  {
    val_ = other;
    return *this;
  }

  Rational_function_flint& operator+=(
                                        const Rational_function_flint& rhs) 
  {
    val_ += rhs();
    return *this; // return the result by reference
  }
  Rational_function_flint& operator-=(
                                        const Rational_function_flint& rhs) 
  {
    val_ -= rhs();
    return *this; // return the result by reference
  }
  Rational_function_flint& operator*=(
                                        const Rational_function_flint& rhs) 
  {
    val_ *= rhs();
    return *this; // return the result by reference
  }
  Rational_function_flint& operator/=(
                                        const Rational_function_flint& rhs) 
  {
    if(rhs().is_zero()) { std::cout << "DIVISION BY ZERO.\n"; return *this; }
    val_ /= rhs();
    return *this; // return the result by reference
  }
};



  inline 
  Rational_function_flint operator+(
                                      Rational_function_flint  lhs,
                                const Rational_function_flint& rhs)  
  {
    // Rational_function_flint res(*this);
    lhs += rhs;
    return lhs;
  }
  inline
  Rational_function_flint operator-(
                                      Rational_function_flint  lhs,
                                const Rational_function_flint& rhs) 
  {
    // Rational_function_flint res(*this);
    lhs -= rhs;
    return lhs;
  }
  inline
  Rational_function_flint operator*(
                                      Rational_function_flint  lhs,
                                const Rational_function_flint& rhs)
  {
    // Rational_function_flint res(1);//*this);
    lhs *= rhs;
    // res *= other;
    // Rational_function_flint res;
    // res.assign( (*this)() * other() );
    return lhs;
  }
  inline
  Rational_function_flint operator/(
                                      Rational_function_flint  lhs,
                                const Rational_function_flint& rhs)
  {
    // Rational_function_flint res(*this);
    if(rhs().is_zero()) { std::cout << "DIVISION NY ZERO.\n"; return lhs;}

    lhs /= rhs;
    return lhs;
  }




std::ostream& operator<<( std::ostream& out
                        , const Rational_function_flint& f)
// { return out << f().pretty("t"); }
{ return out << f(); }
bool operator==( const Rational_function_flint& lhs
               , const Rational_function_flint& rhs)
{ return lhs() == rhs(); }
// bool operator==( const Rational_function_flint& lhs
//                , const slong& rhs)
// { 
//   Rational_function_flint rhs_qfrac(rhs);
//   return lhs == rhs_qfrac; 
// }
bool operator!=( const Rational_function_flint& lhs
               , const Rational_function_flint& rhs)
{ return !(lhs == rhs); }
// bool operator!=( const Rational_function_flint& lhs
//                , const slong& rhs)
// { return !(lhs == rhs); }

/** Returns a Rational_function_flint whose value is z^d. **/
Rational_function_flint pow( Rational_function_flint z
                                    , unsigned int d) 
{
  Rational_function_flint::Fraction qfrac;
  qfrac = z();
  qfrac = pow(qfrac,d);
  Rational_function_flint z_pow_d(qfrac);
  // z_pow_d.assign(pow(z(),d));
  // z.assign(pow(z(),d));
  return z_pow_d;
}


namespace Eigen {
template<> struct NumTraits<Rational_function_flint> 
                   : GenericNumTraits<Rational_function_flint>
{
    typedef Rational_function_flint Real;
    typedef Rational_function_flint NonInteger;
    typedef Rational_function_flint Nested;
    typedef slong Literal;

    static inline Real epsilon() { 
      // flint::fmpz_poly_qxx zero; zero = 0;
      // return zero; 
      return 0;
      // return Rational_function_flint(0); 
    }
    static inline Real dummy_precision() { 
      return 0;
      // Rational_function_flint zero(0); 
      // return Rational_function_flint(0); 
   }
    static inline int digits10() { return 0; }
 
    enum {
      IsInteger = 1,
      IsComplex = 0,
      IsSigned = 1,
      RequireInitialization = 1,
      ReadCost = 50,
      AddCost = 100,
      MulCost = 400
    };
};
 
}
// namespace Eigen {
//   template<> struct NumTraits<mpq_class> : GenericNumTraits<mpq_class>
//   {
//     typedef mpq_class Real;
//     typedef mpq_class NonInteger;
//     typedef mpq_class Nested;
 
//     static inline Real epsilon() { return 0; }
//     static inline Real dummy_precision() { return 0; }
//     static inline int digits10() { return 0; }
 
//     enum {
//       IsInteger = 0,
//       IsSigned = 1,
//       IsComplex = 0,
//       RequireInitialization = 1,
//       ReadCost = 6,
//       AddCost = 150,
//       MulCost = 100
//     };
//   };
// }



// namespace flint {//extra functions on flint::fmpz_poly_qxx, required by Eigen
 
// // inline const adouble& conj(const adouble& x)  { return x; }
// // inline const adouble& real(const adouble& x)  { return x; }
// // inline adouble imag(const adouble&)    { return 0.; }
// // inline adouble abs(const adouble&  x)  { return fabs(x); }
// // inline adouble abs2(const adouble& x)  { return x*x; }
 
// }
 
#endif // RATIONAL_FUNCTION_EXACT_INTEGERS_
