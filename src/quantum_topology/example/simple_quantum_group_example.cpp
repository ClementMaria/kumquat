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
#include <kumquat/Quantum_group.h>

using namespace kumquat;

int main() {
	int num_strands = 4;
	int max_twist = 10;
	std::cout << "Initialization of Uq(sl2(C) with num_strands = " << num_strands << " and max_twist = " << max_twist << "\n";
	Quantum_group_Uqsl2_gen_q<2> uq(num_strands,max_twist);
	
	
	std::vector< std::pair<int,int> > br;
	br.emplace_back(2,3); //braid s_3^{-1} s_1^2 s_2^3
	br.emplace_back(1,2);
	br.emplace_back(-3,1);
	Braid b(num_strands,br);

	std::cout << uq.quantum_invariant(b).to_string() << "\n";

	std::vector< std::map<int,int> > pbr(3, std::map<int,int>() );
	pbr[0][2] = 3;//layer 0: s_2^3
	pbr[1][1] = 2; pbr[1][-3] = 2; //layer 1: s_1^2 s_3^{-2}
	pbr[2][-2] = 1;//layer 0: s_2^{-1}

	Plat_braid pb(num_strands, pbr);

	std::cout << uq.quantum_invariant(pb).to_string() << "\n";

	return 0;
}