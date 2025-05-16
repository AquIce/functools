#include "src/function.hpp"

#include <iostream>
#include <memory>

int main(int argc, char** argv) {
	
	auto func = std::make_shared<functools::PolynomialFunction>(
		2,
		functools::Upcast<functools::ConstantFunction>(
			functools::CoeffsToConstFunctions(
				std::vector<Type>({1, 0, -1})
			)
		)
	);

	auto func2 = std::make_shared<functools::PolynomialFunction>(
		1,
		functools::Upcast<functools::ConstantFunction>(
			functools::CoeffsToConstFunctions(
				std::vector<Type>({1, 1})
			)
		)
	);

	std::cout << func->Repr() << "\n";
	std::cout << func2->Repr() << "\n";

	functools::DivisionResult res = func / func2;

	std::cout << "Quotient: " << "\n";
	std::cout << res.Quotient->Repr() << "\n";

	std::cout << "Remainder: " << "\n";
	std::cout << res.Remainder->Repr() << "\n";

	return 0;
}
