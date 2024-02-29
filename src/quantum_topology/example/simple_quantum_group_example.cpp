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
	
	std::cout << "Quantum invariant of the unknot: " << uq.quantum_invariant_unknot().to_string() << "\n";

	uq.display();



	Plat_braid pbu(num_strands);
	pbu.add_twist(0,1);
	pbu.add_twist(1,1);
	pbu.add_twist(2,1);

	std::cout << "X\n";

	auto q_pbu = uq.quantum_invariant(pbu);
	q_pbu /= uq.quantum_invariant_unknot();

	std::cout << "Quantum invariant of unknot = " << q_pbu << "\n";





	// std::vector< std::pair<int,int> > br;
	// br.emplace_back(2,3); //braid s_3^{-1} s_1^2 s_2^3
	// br.emplace_back(1,2);
	// br.emplace_back(0,-1);
	// Braid b(num_strands,br);

	// std::cout << "Braid b = " << b << "\n";

	// auto q_b = uq.quantum_invariant(b);
	// q_b /= uq.quantum_invariant_unknot();
	// std::cout << q_b.to_string() << "\n";
	
	// std::vector< std::map<int,int> > pbr(3, std::map<int,int>() );
	// pbr[0][1] = 3;//layer 0: s_1^3
	// pbr[1][0] = 2; pbr[1][2] = -2; //layer 1: s_0^2 s_2^{-2}
	// pbr[2][1] = -1;//layer 0: s_1^{-1}

	// Plat_braid pb(num_strands, pbr);

	// std::cout << "Braid pb = " << pb << "\n";

	// std::cout << uq.quantum_invariant(pb).to_string() << "\n";


	// std::cout << "X\n";

	// Plat_braid pbu(num_strands);
	// pbu.add_twist(0,1);
	// pbu.add_twist(1,1);
	// pbu.add_twist(2,1);

	// std::cout << "X\n";

	// auto q_pbu = uq.quantum_invariant(pbu);
	// q_pbu /= uq.quantum_invariant_unknot();

	// std::cout << "Quantum invariant of unknot = " << q_pbu << "\n";
	return 0;
}