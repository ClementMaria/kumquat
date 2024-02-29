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

#include <iostream>
#include <kumquat/Dense_matrix.h>
#include <kumquat/Z_mod_nZ.h>
// #include <tbb/tbb.h>
// #include <omp.h>
#include <tbb/parallel_for.h>

#include <algorithm>
#include <vector>
// #include <execution>

using namespace kumquat;

int main() {
//tbb
  tbb::parallel_for(size_t(0), size_t(10000), [&](size_t i){
    std::cout << i << " ";
  });

    // tbb::parallel_for(0, size, 
    //   [&](int idx){
    //     new (&complex_[idx]) Hasse_simp(cpx, cpx.simplex(idx));
    //   });
    // for (int idx=0; idx < size; ++idx)
    //   if (complex_[idx].boundary_.empty())
    //     vertices_.push_back(idx);

//  std::execution::seq,
// std::execution::par,
// std::execution::par_unseq, and
// std::execution::unseq


// {
//   std::vector<int> v(1000000,0);
//   int i=0;
//   for(auto it=v.begin(); it!=v.end(); ++it) {
//     *it = i++;
//   }
//   // increment elements in-place
//   std::for_each(std::execution::seq, v.begin(), v.end(), [](int &n) { 
//     auto cp = n;
//     for(int i=0; i<cp; ++i) { n += i; }
//   });
// } 


// {
//   std::vector<int> v(1000000,0);
//   int i=0;
//   for(auto it=v.begin(); it!=v.end(); ++it) {
//     *it = i++;
//   }
//   // increment elements in-place
//   std::for_each(v.begin(), v.end(), [](int &n) { 
//     auto cp = n;
//     for(int i=0; i<cp; ++i) { n += i; }
//   });
// } 

 
    // std::cout << "after:\t";
    // // std::for_each(v.cbegin(), v.cend(), print);
    // std::cout << '\n';
 
    // struct Sum
    // {
    //     void operator()(int n) { sum += n; }
    //     int sum {0};
    // };
 
    // // invoke Sum::operator() for each element
    // Sum s = std::for_each(v.cbegin(), v.cend(), Sum());    
    // std::cout << "sum:\t" << s.sum << '\n';


  // #pragma omp parallel for
  // for(int i=0; i<100; ++i) {
  //   std::cout << i << "\n";
  // }

 return 0;
}