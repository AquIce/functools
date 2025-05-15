#include "src/function.hpp"

#include <iostream>
#include <memory>

int main(int argc, char** argv) {
	
	auto func = std::make_shared<functools::PolynomialFunction>(
		4,
		functools::Upcast<functools::ConstantFunction>(
			functools::CoeffsToConstFunctions(
				std::vector<Type>({5, 4, 3, 2, 1})
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

	func = func2 * func;

	std::cout << func->Evaluate(10) << "\n";
	std::cout << func->Repr() << "\n";

	/*std::unique_ptr<functools::Function> funcPrime = func->GetDerivative();

	std::cout << funcPrime->Evaluate(10) << "\n";
	std::cout << funcPrime->Repr() << "\n";

	std::unique_ptr<functools::Function> Func = func->GetPrimitive();

	std::cout << Func->Evaluate(10) << "\n";
	std::cout << Func->Repr() << "\n";
*/
	return 0;
}
