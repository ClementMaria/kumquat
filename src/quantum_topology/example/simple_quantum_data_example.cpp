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
#include <kumquat/Quantum_data.h>
#include <kumquat/Rational_function_integral_mp.h>

using namespace kumquat;

int main() {

	Quantum_data<Rational_function_integral_mp> qd;
	std::cout << "Quantum half X^2-X^{-2} = " << qd.quantum_half() << "\n\n";
	std::cout << "Quantum monomials X^k, and 1 if k<0: \n";
	for(int k=-1; k<=10; ++k) {
		std::cout << "   k = " << k << "  ->  " << qd.quantum_monomial(k) << "\n";
	}
	std::cout << "\n";
	std::cout << "Quantum integers [n] = (X^{2n}-X^{-2n}) / (X^2 - X^{-2}): \n";
	for(int n=-1; n<=10; ++n) {
		std::cout << "   n = " << n << "  ->  " << qd.quantum_integer(n) << "\n";
	}
	std::cout << "\n";

	std::cout << "Quantum factorials [n]!: \n";
	for(int n=-1; n<=10; ++n) {
		std::cout << "   n = " << n << "  ->  " << qd.quantum_factorial(n) << "\n";
	}
	std::cout << "\n";

	return 0;
}