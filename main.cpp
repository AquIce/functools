#include "src/function.hpp"

#include <iostream>
#include <memory>

int main(int argc, char** argv) {
	auto func = std::make_shared<functools::TrigonometryFunction>(
		functools::TrigonometryFunctionType::COT,
		std::make_shared<functools::PolynomialFunction>(
			3,
			functools::Upcast<functools::ConstantFunction>(
				functools::CoeffsToConstFunctions(
					std::vector<Type>({3, 2, 1, 0})
				)
			)
		)
	);

	auto func2 = std::make_shared<functools::PolynomialFunction>(
		5,
		functools::Upcast<functools::ConstantFunction>(
			functools::CoeffsToConstFunctions(
				std::vector<Type>({1, 2, 3, 4, 5, 6})
			)
		)
	);

	std::shared_ptr<functools::PolynomialFunction> func3 = func * func2;

	auto fprime = func3->GetDerivative();

	LOG("f");
	LOG(func3->Repr());
	LOG("f'");
	LOG(fprime->Repr());

	return 0;
}
