#include "src/function.hpp"

#include <iostream>
#include <memory>

int main(int argc, char** argv) {

	// cot(3x^3 + 2x^2 + x)
	auto func = std::make_shared<functools::TrigonometryFunction>(
		functools::TrigonometryFunctionType::COT,
		std::make_shared<functools::ComplexFunction>(

			std::make_shared<functools::ComplexFunction>(
				std::make_shared<functools::ConstantFunction>(3),
				functools::FunctionOperator::TIMES,
				std::make_shared<functools::ComplexFunction>(
					std::make_shared<functools::IdentityFunction>(),
					functools::FunctionOperator::POWER,
					std::make_shared<functools::ConstantFunction>(3)
				)
			),

			functools::FunctionOperator::PLUS,
			
			std::make_shared<functools::ComplexFunction>(
				std::make_shared<functools::ComplexFunction>(
					std::make_shared<functools::ConstantFunction>(2),
					functools::FunctionOperator::TIMES,
					std::make_shared<functools::ComplexFunction>(
						std::make_shared<functools::IdentityFunction>(),
						functools::FunctionOperator::POWER,
						std::make_shared<functools::ConstantFunction>(2)
					)
				),

				functools::FunctionOperator::PLUS,

				std::make_shared<functools::IdentityFunction>()
			)
		)
	);

	// x^5 + 2x^4 + 3x^3 + 4x^2 + 5x + 6
	auto func2 = std::make_shared<functools::ComplexFunction>(

		std::make_shared<functools::ComplexFunction>(
			std::make_shared<functools::IdentityFunction>(),
			functools::FunctionOperator::POWER,
			std::make_shared<functools::ConstantFunction>(5)
		),

		functools::FunctionOperator::PLUS,
		
		std::make_shared<functools::ComplexFunction>(
			std::make_shared<functools::ComplexFunction>(
				std::make_shared<functools::ConstantFunction>(2),
				functools::FunctionOperator::TIMES,
				std::make_shared<functools::ComplexFunction>(
					std::make_shared<functools::IdentityFunction>(),
					functools::FunctionOperator::POWER,
					std::make_shared<functools::ConstantFunction>(4)
				)
			),

			functools::FunctionOperator::PLUS,

			std::make_shared<functools::ComplexFunction>(
				std::make_shared<functools::ComplexFunction>(
					std::make_shared<functools::ConstantFunction>(3),
					functools::FunctionOperator::TIMES,
					std::make_shared<functools::ComplexFunction>(
						std::make_shared<functools::IdentityFunction>(),
						functools::FunctionOperator::POWER,
						std::make_shared<functools::ConstantFunction>(3)
					)
				),

				functools::FunctionOperator::PLUS,

				std::make_shared<functools::ComplexFunction>(
					std::make_shared<functools::ComplexFunction>(
						std::make_shared<functools::ConstantFunction>(2),
						functools::FunctionOperator::TIMES,
						std::make_shared<functools::ComplexFunction>(
							std::make_shared<functools::IdentityFunction>(),
							functools::FunctionOperator::POWER,
							std::make_shared<functools::ConstantFunction>(4)
						)
					),

					functools::FunctionOperator::PLUS,

					std::make_shared<functools::ComplexFunction>(
						std::make_shared<functools::ComplexFunction>(
							std::make_shared<functools::ConstantFunction>(5),
							functools::FunctionOperator::TIMES,
							std::make_shared<functools::IdentityFunction>()
						),

						functools::FunctionOperator::PLUS,

						std::make_shared<functools::ConstantFunction>(6)
					)
				)
			)
		)
	);

	std::shared_ptr<functools::ComplexFunction> func3 = func * func2;
	auto func4 = std::make_shared<functools::LogarithmFunction>(
		std::make_shared<functools::ConstantFunction>(std::exp(1)),
		func3
	);

	auto fprime = func4->GetDerivative();

	LOG("f");
	LOG(func4->Repr());
	LOG("f'");
	LOG(fprime->Repr());

	return 0;
}
