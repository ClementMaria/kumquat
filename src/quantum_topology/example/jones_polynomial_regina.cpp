#include <iostream>
#include <fstream> 
#include <cmath>
#include <iomanip>
// #include <cppunit/extensions/HelperMacros.h>
// #include "census/census.h"
// #include "link/examplelink.h"
#include "link/link.h"
// #include "maths/laurent.h"
// #include "maths/laurent2.h"
// #include "packet/container.h"
// #include "surface/normalsurface.h"
// #include "surface/normalsurfaces.h"
// #include "triangulation/dim3.h"

// using regina::Crossing;
// using regina::ExampleLink;
using regina::Link;
// using regina::Triangulation;
// using regina::StrandRef;
// using regina::Edge;


/* Convert the Jones polynomial from Regina to Bar-Nathan convnetion in "On Khovanov categorification of the Jones polynomial", which just consists of replacing the varaible x from the Regina output by (-q).
*/
template<class LaurentPoly>
void convert_jones_regina_to_bar_nathan(LaurentPoly &Lp) {
  typename LaurentPoly::Coefficient c;
  for(int i = Lp.minExp(); i <= Lp.maxExp(); ++i) {
    if( (i%2) == 1) { Lp.set(i,Lp[i]*(-1)); }
  }  
}



int main() {

  std::string gausscode = "-<2 +<7 ->8 +<9 ->10 +<11 ->12 +<13 ->14 +<15 -<16 +>17 -<18 +>19 -<20 +>21 -<22 +>23 -<24 +>27 -<28 +>29 -<30 +>31 +>1 +>2 +<26 -<1 +<3 ->4 +<5 ->6 ->7 +<8 ->9 +<10 ->11 +<12 ->13 +<14 ->15 -<25 -<27 +>28 -<29 +>30 -<31 ->3 +<4 ->5 +<6 +>16 -<17 +>18 -<19 +>20 -<21 +>22 -<23 +>24 +>25 ->26";
  //compute the jones polynomial
  Link curr_link(gausscode);
  // std::ostringstream s2; 
  auto jones_poly = curr_link.jones();//regina::ALG_TREEWIDTH);
  // std::cout << jones_poly << "\n";
  convert_jones_regina_to_bar_nathan(jones_poly);

  std::cout << jones_poly << "\n";

  return 0;
}
