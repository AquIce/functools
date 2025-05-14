#include "src/function.hpp"

#include <iostream>
#include <memory>

int main(int argc, char** argv) {
	
	auto func = std::make_shared<functools::PolynomialFunction<4>>(
		functools::Upcast<functools::ConstantFunction, 5>(
			functools::CoeffsToConstFunctions<5>(
				std::array<Type, 5>({5, 4, 3, 2, 1})
			)
		)
	);

	auto func2 = std::make_shared<functools::PolynomialFunction<3>>(
		functools::Upcast<functools::ConstantFunction, 4>(
			functools::CoeffsToConstFunctions<4>(
				std::array<Type, 4>({2, 1, 0, -1})
			)
		)
	);

	func = func2 + func;

	// auto func = std::make_shared<functools::PolynomialFunction<3>>(
	// 	std::array<std::shared_ptr<functools::Function>, 4>({
	// 		std::dynamic_pointer_cast<functools::Function>(
	// 			std::make_shared<functools::PolynomialFunction<1>>(
	// 				functools::Upcast<functools::ConstantFunction, 2>(
	// 					functools::CoeffsToConstFunctions<2>(
	// 						std::array<Type, 2>({1, -1})
	// 					)
	// 				)
	// 			)
	// 		),
	// 		std::dynamic_pointer_cast<functools::Function>(
	// 			std::make_shared<functools::ConstantFunction>(5)
	// 		),
	// 		std::dynamic_pointer_cast<functools::Function>(
	// 			std::make_shared<functools::PolynomialFunction<3>>(
	// 				functools::Upcast<functools::ConstantFunction, 4>(
	// 					functools::CoeffsToConstFunctions<4>(
	// 						std::array<Type, 4>({3, -2, 1, 0})
	// 					)
	// 				)
	// 			)
	// 		),
	// 		std::dynamic_pointer_cast<functools::Function>(
	// 			std::make_shared<functools::ConstantFunction>(7)
	// 		)
	// 	})
	// );

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
