#include "src/function.hpp"

#include <iostream>
#include <memory>

int main(int argc, char** argv) {
	auto func = std::make_shared<functools::TrigonometryFunction>(
		functools::TrigonometryFunctionType::COS,
		std::make_shared<functools::PolynomialFunction>(
			1,
			functools::Upcast<functools::ConstantFunction>(
				functools::CoeffsToConstFunctions(
					std::vector<Type>({1, 0})
				)
			)
		)
	);

	auto func2 = std::make_shared<functools::PolynomialFunction>(
		1,
		functools::Upcast<functools::ConstantFunction>(
			functools::CoeffsToConstFunctions(
				std::vector<Type>({1, 2})
			)
		)
	);

	std::shared_ptr<functools::PolynomialFunction> func3 = func * func2;

	LOG(func3->Repr());
	LOG(func3->GetDerivative()->Repr());

	return 0;
}
