#include "src/function.hpp"

#include <iostream>
#include <memory>

int main(int argc, char** argv) {
	
	auto func = std::make_shared<functools::PolynomialFunction>(
		1,
		std::vector<std::shared_ptr<functools::Function>>({
			std::dynamic_pointer_cast<functools::Function>(
				std::make_shared<functools::PolynomialFunction>(
					1,
					functools::Upcast<functools::ConstantFunction>(
						functools::CoeffsToConstFunctions(
							std::vector<Type>({1, 2})
						)
					)
				)
			),
			std::dynamic_pointer_cast<functools::Function>(
				std::make_shared<functools::ConstantFunction>(3)
			)
		})
	);

	auto func2 = std::make_shared<functools::PolynomialFunction>(
		1,
		functools::Upcast<functools::ConstantFunction>(
			functools::CoeffsToConstFunctions(
				std::vector<Type>({1, 2})
			)
		)
	);

	std::cout << func->Repr() << "\n";
	std::cout << func->GetDerivative()->Repr() << "\n";
	std::cout << func->GetPrimitive()->Repr() << "\n";

	std::cout << func2->Repr() << "\n";
	std::cout << func2->GetDerivative()->Repr() << "\n";
	std::cout << func2->GetPrimitive()->Repr() << "\n";

	return 0;
}
