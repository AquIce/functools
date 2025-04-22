#include "src/function.hpp"

#include <iostream>
#include <memory>

int main(int argc, char** argv) {
	
	auto func = std::make_unique<functools::PolynomialFunction<4>>(
		std::array<Type, 5>({1, 2, 3, 4, 5})
	);

	std::cout << func->Evaluate(10) << "\n";
	std::cout << func->Repr() << "\n";

	std::unique_ptr<functools::Function> funcPrime = func->GetDerivative();

	std::cout << funcPrime->Evaluate(10) << "\n";
	std::cout << funcPrime->Repr() << "\n";

	std::unique_ptr<functools::Function> Func = func->GetPrimitive();

	std::cout << Func->Evaluate(10) << "\n";
	std::cout << Func->Repr() << "\n";

	return 0;
}
