#pragma once

#include "binds.hpp"
#include "debug.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>

#define DegreeType uint8_t
#define Type double

template <DegreeType lhs, DegreeType rhs>
struct MaxDegree {
public:
	static const DegreeType value = (lhs > rhs) ? lhs : rhs;
};

namespace functools {

	enum class FunctionType {
		NONE = 0,
		CONSTANT,
		POLYNOMIAL,
		TRIGONOMETRY
	};

	// ---
	// Prototypes
	// ---

	class Function {
	public:

		virtual Type Evaluate(Type x);

		// virtual std::shared_ptr<Function> operator+(Type num);
		// virtual std::shared_ptr<Function> operator+(std::shared_ptr<Function> other);
		// virtual std::shared_ptr<Function> operator-(Type num);
		// virtual std::shared_ptr<Function> operator-(std::shared_ptr<Function> other);
		// virtual std::shared_ptr<Function> operator*(Type num);
		// virtual std::shared_ptr<Function> operator*(std::shared_ptr<Function> other);
		// virtual std::shared_ptr<Function> operator/(Type num);
		// virtual std::shared_ptr<Function> operator/(std::shared_ptr<Function> other);

		virtual std::shared_ptr<Function> GetDerivative() const = 0; 
		virtual std::shared_ptr<Function> GetPrimitive() const = 0;

		virtual std::string Repr() const = 0;

		virtual FunctionType GetType() const;

		virtual DegreeType GetPolynomialDegree() const = 0;
		virtual std::shared_ptr<Function> GetPolynomialCoefficientAt(DegreeType index) const = 0;
	};

	class ConstantFunction : public Function {
	public:
		ConstantFunction(Type value);

		Type Evaluate(Type x) override;

		// std::shared_ptr<Function> operator+(Type num) override;
		// std::shared_ptr<Function> operator+(std::shared_ptr<Function> other) override;
		// std::shared_ptr<Function> operator-(Type num) override;
		// std::shared_ptr<Function> operator-(std::shared_ptr<Function> other) override;
		// std::shared_ptr<Function> operator*(Type num) override;
		// std::shared_ptr<Function> operator*(std::shared_ptr<Function> other) override;
		// std::shared_ptr<Function> operator/(Type num) override;
		// std::shared_ptr<Function> operator/(std::shared_ptr<Function> other) override;

		std::shared_ptr<Function> GetDerivative() const override;

		std::shared_ptr<Function> GetPrimitive() const override;

		std::string Repr() const override;

		FunctionType GetType() const override;	
		
		Type GetValue();

		DegreeType GetPolynomialDegree() const override;
		std::shared_ptr<Function> GetPolynomialCoefficientAt(DegreeType index) const override;
	
	private:
		Type value;
	};

	template <DegreeType Degree = 0>
	class PolynomialFunction : public Function {
	
	public:

		PolynomialFunction(
			std::array<std::shared_ptr<Function>, Degree + 1> coefficients
		);

		Type Evaluate(Type x) override;

		// std::shared_ptr<Function> operator+(Type num) override;
		// std::shared_ptr<Function> operator+(std::shared_ptr<Function> other) override;
		// std::shared_ptr<Function> operator-(Type num) override;
		// std::shared_ptr<Function> operator-(std::shared_ptr<Function> other) override;
		// std::shared_ptr<Function> operator*(Type num) override;
		// std::shared_ptr<Function> operator*(std::shared_ptr<Function> other) override;
		// std::shared_ptr<Function> operator/(Type num) override;
		// std::shared_ptr<Function> operator/(std::shared_ptr<Function> other) override;

		std::shared_ptr<Function> GetDerivative() const override;

		std::shared_ptr<Function> GetPrimitive() const override;

		std::string Repr() const override;

		FunctionType GetType() const override;

		std::array<std::shared_ptr<Function>, Degree + 1> GetCoefficients();

		DegreeType GetPolynomialDegree() const override;
		std::shared_ptr<Function> GetPolynomialCoefficientAt(DegreeType index) const override;

	private:
		std::array<std::shared_ptr<Function>, Degree + 1> coefficients;
	};

	template <typename T, size_t S>
	std::array<std::shared_ptr<Function>, S> Upcast(
		std::array<std::shared_ptr<T>, S> funcs
	) {
		std::array<std::shared_ptr<Function>, S> res;
		
		std::transform(
			funcs.begin(), funcs.end(),
			res.begin(),
			[](std::shared_ptr<T> ptr) -> std::shared_ptr<Function> {
				return std::dynamic_pointer_cast<Function>(ptr);
			}
		);

		return res;
	}

	// ---
	// General Function
	// ---

	Type Function::Evaluate(Type x) { throw std::runtime_error("Invalid function"); }

	// std::shared_ptr<Function> Function::operator+(Type num) { throw std::runtime_error("Invalid function"); }
	// std::shared_ptr<Function> Function::operator+(std::shared_ptr<Function> other) { throw std::runtime_error("Invalid function"); }
	// std::shared_ptr<Function> Function::operator-(Type num) { throw std::runtime_error("Invalid function"); }
	// std::shared_ptr<Function> Function::operator-(std::shared_ptr<Function> other) { throw std::runtime_error("Invalid function"); }
	// std::shared_ptr<Function> Function::operator*(Type num) { throw std::runtime_error("Invalid function"); }
	// std::shared_ptr<Function> Function::operator*(std::shared_ptr<Function> other) { throw std::runtime_error("Invalid function"); }
	// std::shared_ptr<Function> Function::operator/(Type num) { throw std::runtime_error("Invalid function"); }
	// std::shared_ptr<Function> Function::operator/(std::shared_ptr<Function> other) { throw std::runtime_error("Invalid function"); }

	FunctionType Function::GetType() const {
		return FunctionType::NONE;
	}

	// ---
	// Helper Functions
	// ---

	template <size_t S>
	std::array<std::shared_ptr<ConstantFunction>, S> CoeffsToConstFunctions(
		std::array<Type, S> values
	) {
		std::array<std::shared_ptr<ConstantFunction>, S> constFunctions;
		for(size_t i = 0; i < values.size(); i++) {
			constFunctions.at(i) = std::make_shared<ConstantFunction>(values.at(i));
		}
		return constFunctions;
	}

	// ---
	// Constant Functions
	// ---

	ConstantFunction::ConstantFunction(Type value) :
		Function(),
		value(value)
	{}

	Type ConstantFunction::Evaluate(Type x) {
		return value;
	}

	// std::shared_ptr<Function> ConstantFunction::operator+(Type num) {
	// 	return std::make_shared<ConstantFunction>(value + num);
	// }
	// std::shared_ptr<Function> ConstantFunction::operator+(std::shared_ptr<Function> other) {
	// 	return *other.get() + value;
	// }
	// std::shared_ptr<Function> ConstantFunction::operator-(Type num) {
	// 	return std::make_shared<ConstantFunction>(value - num);
	// }
	// std::shared_ptr<Function> ConstantFunction::operator-(std::shared_ptr<Function> other) {
	// 	return *other.get() - value;
	// }
	// std::shared_ptr<Function> ConstantFunction::operator*(Type num) {
	// 	return std::make_shared<ConstantFunction>(value * num);
	// }
	// std::shared_ptr<Function> ConstantFunction::operator*(std::shared_ptr<Function> other) {
	// 	return *other.get() * value;
	// }
	// std::shared_ptr<Function> ConstantFunction::operator/(Type num) {
	// 	return std::make_shared<ConstantFunction>(value / num);
	// }
	// std::shared_ptr<Function> ConstantFunction::operator/(std::shared_ptr<Function> other) {
	// 	return *other.get() / value;
	// }

	std::shared_ptr<Function> ConstantFunction::GetDerivative() const {
		return std::make_shared<ConstantFunction>(0);
	}

	std::shared_ptr<Function> ConstantFunction::GetPrimitive() const {
		return std::make_shared<PolynomialFunction<1>>(
			functools::Upcast<functools::ConstantFunction, 2>(
				CoeffsToConstFunctions<2>(
					std::array<Type, 2>({ value, 0 })
				)
			)
		);
	}

	std::string ConstantFunction::Repr() const {
		return std::to_string(value);
	}

	FunctionType ConstantFunction::GetType() const {
		return FunctionType::CONSTANT;
	}

	Type ConstantFunction::GetValue() {
		return value;
	}

	DegreeType ConstantFunction::GetPolynomialDegree() const {
		throw std::runtime_error("Trying to get ConstantFunction's polynomial degree");
	}

	std::shared_ptr<Function> ConstantFunction::GetPolynomialCoefficientAt(DegreeType index) const {
		throw std::runtime_error("Trying to get one of the polynomial coefficients of a ConstantFunction");
	}

	// --- 
	// Polynomial Functions
	// ---

	template <DegreeType Degree>
	PolynomialFunction<Degree>::PolynomialFunction(
		std::array<std::shared_ptr<Function>, Degree + 1> coefficients
	) :
		Function(),
		coefficients(coefficients)
	{}

	template <DegreeType Degree>
	Type PolynomialFunction<Degree>::Evaluate(Type x) {
		Type result = 0;
		DegreeType currentDegree = Degree;

		for(const auto& coefficient : coefficients) {
			result += coefficient->Evaluate(x) * std::pow(x, currentDegree--);
		}
		return result;
	}

	// template <DegreeType Degree>
	// std::shared_ptr<Function> PolynomialFunction<Degree>::operator+(Type num) {
	// 	std::array<std::shared_ptr<Function>, Degree + 1> coeffs;
	// 	std::copy(coefficients.begin(), coefficients.end(), coeffs.begin());

	// 	coeffs.at(Degree) = *(coefficients.at(Degree).get()) + num;
	// 	return std::make_shared<PolynomialFunction<Degree>>(coeffs);
	// }

	// template <DegreeType Degree>
	// std::shared_ptr<Function> PolynomialFunction<Degree>::operator+(std::shared_ptr<Function> other) {
	// 	switch(other->GetType()) {

	// 		case FunctionType::NONE: {
	// 			throw std::runtime_error("Invalid function");
	// 		}

	// 		case FunctionType::CONSTANT: {
	// 			return *this + std::dynamic_pointer_cast<ConstantFunction>(other)->GetValue();
	// 		}

	// 		case FunctionType::POLYNOMIAL: {

	// 			DegreeType otherDegree = other->GetPolynomialDegree();

	// 			if(Degree == otherDegree) {

	// 				std::array<std::shared_ptr<Function>, Degree + 1> coeffs;

	// 				for(size_t i = 0; i < Degree + 1; i++) {
	// 					coeffs.at(i) = *(coefficients.at(i).get()) + other->GetPolynomialCoefficientAt(i);
	// 				}

	// 				return std::make_shared<PolynomialFunction<Degree>>(coeffs);
	// 			}

	// 			if(Degree > otherDegree) {
	// 				std::array<std::shared_ptr<Function>, Degree + 1> coeffs;
	// 				std::copy(coefficients.begin(), coefficients.end(), coeffs.begin());

	// 				for(DegreeType i = 0; i <= otherDegree; i++) {
	// 					coeffs.at(Degree - otherDegree + i) =
	// 						*(coeffs.at(Degree - otherDegree + i).get())
	// 						+ other->GetPolynomialCoefficientAt(i);
	// 				}

	// 				return std::make_shared<PolynomialFunction<Degree>>(
	// 					coeffs
	// 				);
	// 			}

	// 			throw std::runtime_error("To add two polynomial functions, you need to put the higher degree one first.");

	// 			// // TODO : Fix here (non-const)
	// 			// std::array<std::shared_ptr<Function>, otherDegree + 1> coeffs;
	// 			// for(DegreeType i = 0; i <= otherDegree; i++) {
	// 			// 	coeffs.at(i) = other->GetPolynomialCoefficientAt(i);
	// 			// }

	// 			// for(DegreeType i = 0; i <= Degree; i++) {
	// 			// 	coeffs.at(otherDegree - Degree + i) =
	// 			// 		*(coeffs.at(otherDegree - Degree + i).get())
	// 			// 		+ coefficients.at(i);
	// 			// }

	// 			// return std::make_shared<PolynomialFunction<Degree>>(
	// 			// 	coeffs
	// 			// );
	// 		}

	// 		case FunctionType::TRIGONOMETRY: {
	// 			NOIMP;
	// 			break;
	// 		}

	// 		default: {
	// 			NOIMP;
	// 			break;
	// 		}
	// 	}
	// }
	// template <DegreeType Degree>
	// std::shared_ptr<Function> PolynomialFunction<Degree>::operator-(Type num) {
	// 	NOIMP;
	// }
	// template <DegreeType Degree>
	// std::shared_ptr<Function> PolynomialFunction<Degree>::operator-(std::shared_ptr<Function> other) {
	// 	NOIMP;
	// }
	// template <DegreeType Degree>
	// std::shared_ptr<Function> PolynomialFunction<Degree>::operator*(Type num) {
	// 	NOIMP;
	// }
	// template <DegreeType Degree>
	// std::shared_ptr<Function> PolynomialFunction<Degree>::operator*(std::shared_ptr<Function> other) {
	// 	NOIMP;
	// }
	// template <DegreeType Degree>
	// std::shared_ptr<Function> PolynomialFunction<Degree>::operator/(Type num) {
	// 	NOIMP;	
	// }
	// template <DegreeType Degree>
	// std::shared_ptr<Function> PolynomialFunction<Degree>::operator/(std::shared_ptr<Function> other) {
	// 	NOIMP;
	// }

	// Add inner derivative
	template <DegreeType Degree>
	std::shared_ptr<Function> PolynomialFunction<Degree>::GetDerivative() const {
		NOIMP;

		/*std::array<std::shared_ptr<Function>, Degree> derivedCoefficients;

		for(DegreeType i = 0; i < Degree; i++) {
			derivedCoefficients.at(i) = coefficients.at(i) * (Degree - i);
		}

		return std::make_shared<PolynomialFunction<Degree - 1>>(derivedCoefficients);
*/		}

	template <DegreeType Degree>
	std::shared_ptr<Function> PolynomialFunction<Degree>::GetPrimitive() const {
		NOIMP;

		// std::array<Type, Degree + 2> primitiveCoefficients;

		// for(DegreeType i = 0; i <= Degree; i++) {
		// 	primitiveCoefficients.at(i) = coefficients.at(i) / (Degree - i + 1);
		// }

		// return std::make_shared<PolynomialFunction<Degree + 1>>(primitiveCoefficients);
	}

	template <DegreeType Degree>
	std::string PolynomialFunction<Degree>::Repr() const {

		std::string repr = "";
		DegreeType currentDegree = Degree;
		for(const auto& coefficient : coefficients) {

			std::string xForm = "x^" + std::to_string(currentDegree);
			
			if(currentDegree == 0) {
				xForm = "";
			} else if(currentDegree == 1) {
				xForm = "x";
			}

			repr +=
				std::string(currentDegree != Degree ? "+ " : "")
				+ "(" + coefficient->Repr() + ") " + xForm
				+ std::string(currentDegree != 0 ? " " : ""); 

			currentDegree--;
		}

		return repr;
	}

	template <DegreeType Degree>
	FunctionType PolynomialFunction<Degree>::GetType() const {
		return FunctionType::POLYNOMIAL;
	}

	template <DegreeType Degree>
	std::array<std::shared_ptr<Function>, Degree + 1> PolynomialFunction<Degree>::GetCoefficients() {
		return coefficients;
	}

	template <DegreeType Degree>
	DegreeType PolynomialFunction<Degree>::GetPolynomialDegree() const {
		return Degree;
	}

	template <DegreeType Degree>
	std::shared_ptr<Function> PolynomialFunction<Degree>::GetPolynomialCoefficientAt(DegreeType index) const {
		LOG(std::string("Getting polynomial coeff at ") + std::to_string(index))
		return coefficients.at(index);
	}

	// For Degree = 0

	template <>
	Type PolynomialFunction<0>::Evaluate(Type x) {
		return coefficients.at(0)->Evaluate(x);
	}

	template <>
	std::shared_ptr<Function> PolynomialFunction<0>::GetDerivative() const {
		NOIMP;
		// return std::make_shared<PolynomialFunction<0>>(
		// 	std::array<Type, 1>({
		// 		static_cast<Type>(0)
		// 	})
		// );
	}

	template <>
	std::string PolynomialFunction<0>::Repr() const {
		return coefficients.at(0)->Repr();
	}

	// For Degree = std::numeric_limits<DegreeType>::max()

	template <>
	std::shared_ptr<Function> PolynomialFunction<std::numeric_limits<DegreeType>::max()>::GetPrimitive() const {
		throw std::runtime_error("Exceeded maximum allowable degree for polynomial functions.");
	}

	// --- 
	// Trigonometry Functions
	// ---

	// enum class TrigonometryFunctionType {
	// 	SIN,
	// 	COS,
	// 	TAN,
	// 	SINH,
	// 	COSH,
	// 	TANH
	// };
	
	// class TrigonometryFunction : public Function {
	
	// public:

	// 	TrigonometryFunction(
	// 		Type coefficient,
	// 		const TrigonometryFunctionType type,
	// 		std::shared_ptr<Function> inner
	// 	) :
	// 		coefficient(coefficient),
	// 		type(type),
	// 		inner(inner)
	// 	{}

	// 	Type Evaluate(Type x) override {
	// 		switch(type) {
	// 			case TrigonometryFunctionType::SIN:
	// 				return static_cast<Type>(std::sin(x));
	// 			case TrigonometryFunctionType::COS:
	// 				return static_cast<Type>(std::cos(x));
	// 			case TrigonometryFunctionType::TAN:
	// 				return static_cast<Type>(std::tan(x));
	// 			case TrigonometryFunctionType::SINH:
	// 				return static_cast<Type>(std::sinh(x));
	// 			case TrigonometryFunctionType::COSH:
	// 				return static_cast<Type>(std::cosh(x));
	// 			case TrigonometryFunctionType::TANH:
	// 				return static_cast<Type>(std::tanh(x));
	// 			default:
	// 				throw std::runtime_error("Invalid trigonometry operation");
	// 		}
	// 	}

	// 	// TODO : Add internal derivative
	// 	std::shared_ptr<Function> GetDerivative() const override {
	// 		switch(type) {
	// 			case TrigonometryFunctionType::SIN:
	// 				return std::make_shared<TrigonometryFunction>(
	// 					1,
	// 					TrigonometryFunctionType::COS
	// 				);
	// 			case TrigonometryFunctionType::COS:
	// 				return std::make_shared<TrigonometryFunction>(
	// 					- coefficient / std::abs(coefficient),
	// 					TrigonometryFunctionType::SIN
	// 				);
	// 			case TrigonometryFunctionType::TAN:
	// 				NOIMP;
	// 			case TrigonometryFunctionType::SINH:
	// 				NOIMP;
	// 			case TrigonometryFunctionType::COSH:
	// 				NOIMP;
	// 			case TrigonometryFunctionType::TANH:
	// 				NOIMP;
	// 			default:
	// 				NOIMP;
	// 		}
	// 	}

	// 	std::shared_ptr<Function> GetPrimitive() const override;

	// 	std::string Repr() const override;

	// 	FunctionType GetType() const override {
	// 		return FunctionType::TRIGONOMETRY;
	// 	}

	// private:
	// 	Type coefficient;
	// 	TrigonometryFunctionType type;
	// 	std::shared_ptr<Function> inner;
	// };

}

// ---
// Operators
// ---

std::shared_ptr<functools::ConstantFunction> operator+(
	std::shared_ptr<functools::ConstantFunction> lhs,
	Type rhs
) {
	return std::make_shared<functools::ConstantFunction>(lhs->GetValue() + rhs);
}
std::shared_ptr<functools::ConstantFunction> operator+(
	std::shared_ptr<functools::ConstantFunction> lhs,
	std::shared_ptr<functools::ConstantFunction> rhs
) {
	return std::make_shared<functools::ConstantFunction>(lhs->GetValue() + rhs->GetValue());
}
template <DegreeType Degree>
std::shared_ptr<functools::PolynomialFunction<Degree>> operator+(
	std::shared_ptr<functools::ConstantFunction> lhs,
	std::shared_ptr<functools::PolynomialFunction<Degree>> rhs
) {
	
}
template <DegreeType Degree>
std::shared_ptr<functools::PolynomialFunction<Degree>> operator+(
	std::shared_ptr<functools::PolynomialFunction<Degree>> lhs,
	Type rhs
) {
	NOIMP;
}
template <DegreeType Degree>
std::shared_ptr<functools::PolynomialFunction<Degree>> operator+(
	std::shared_ptr<functools::PolynomialFunction<Degree>> lhs,
	std::shared_ptr<functools::ConstantFunction> rhs
) {
	NOIMP;
}
template <DegreeType lhsDegree, DegreeType rhsDegree>
std::shared_ptr<functools::PolynomialFunction<MaxDegree<lhsDegree, rhsDegree>::value>> operator+(
	std::shared_ptr<functools::PolynomialFunction<lhsDegree>> lhs,
	std::shared_ptr<functools::PolynomialFunction<rhsDegree>> rhs
) {
	NOIMP;
}