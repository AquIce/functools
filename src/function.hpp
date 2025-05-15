#pragma once

#include "debug.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#define DegreeType uint8_t
#define Type double

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

		virtual ~Function() = default;

		virtual Type Evaluate(Type x);

		virtual std::shared_ptr<Function> GetDerivative() const = 0; 
		virtual std::shared_ptr<Function> GetPrimitive() const = 0;

		virtual std::string Repr() const = 0;

		virtual FunctionType GetType() const;
	};

	class ConstantFunction : public Function {
	public:
		ConstantFunction(Type value);

		Type Evaluate(Type x) override;

		std::shared_ptr<Function> GetDerivative() const override;

		std::shared_ptr<Function> GetPrimitive() const override;

		std::string Repr() const override;

		FunctionType GetType() const override;	
		
		Type GetValue();
	
	private:
		Type m_value;
	};

	class PolynomialFunction : public Function {
	
	public:

		PolynomialFunction(
			DegreeType degree,
			std::vector<std::shared_ptr<Function>> coefficients
		);

		Type Evaluate(Type x) override;

		std::shared_ptr<Function> GetDerivative() const override;

		std::shared_ptr<Function> GetPrimitive() const override;

		std::string Repr() const override;

		FunctionType GetType() const override;

		DegreeType GetDegree();

		std::vector<std::shared_ptr<Function>> GetCoefficients();

	private:
		const DegreeType m_degree;
		std::vector<std::shared_ptr<Function>> m_coefficients;
	};
}

// ---
// Operators Prototypes
// ---

std::shared_ptr<functools::Function> operator+(
	std::shared_ptr<functools::Function> lhs,
	std::shared_ptr<functools::Function> rhs
);
std::shared_ptr<functools::Function> operator+(
	std::shared_ptr<functools::Function> lhs,
	Type rhs
);

std::shared_ptr<functools::Function> operator*(
	std::shared_ptr<functools::Function> lhs,
	std::shared_ptr<functools::Function> rhs
);
std::shared_ptr<functools::Function> operator*(
	std::shared_ptr<functools::Function> lhs,
	Type rhs
);

#include "binds.hpp"

namespace functools {

	template <typename T>
	std::vector<std::shared_ptr<Function>> Upcast(
		std::vector<std::shared_ptr<T>> funcs
	) {
		std::vector<std::shared_ptr<Function>> res;

		for(auto func : funcs) {
			res.push_back(std::dynamic_pointer_cast<Function>(func));
		}

		return res;
	}

	// ---
	// General Function
	// ---

	Type Function::Evaluate(Type x) { throw std::runtime_error("Invalid function"); }

	FunctionType Function::GetType() const {
		return FunctionType::NONE;
	}

	// ---
	// Helper Functions
	// ---

	std::vector<std::shared_ptr<ConstantFunction>> CoeffsToConstFunctions(
		std::vector<Type> values
	) {
		std::vector<std::shared_ptr<ConstantFunction>> constFunctions;
		for(size_t i = 0; i < values.size(); i++) {
			constFunctions.push_back(std::make_shared<ConstantFunction>(values.at(i)));
		}
		return constFunctions;
	}

	// ---
	// Constant Functions
	// ---

	ConstantFunction::ConstantFunction(Type value) :
		Function(),
		m_value(value)
	{}

	Type ConstantFunction::Evaluate(Type x) {
		return m_value;
	}

	std::shared_ptr<Function> ConstantFunction::GetDerivative() const {
		return std::make_shared<ConstantFunction>(0);
	}

	std::shared_ptr<Function> ConstantFunction::GetPrimitive() const {
		return std::make_shared<PolynomialFunction>(
			1,
			functools::Upcast<functools::ConstantFunction>(
				CoeffsToConstFunctions(
					std::vector<Type>({ m_value, 0 })
				)
			)
		);
	}

	std::string ConstantFunction::Repr() const {
		return std::to_string(m_value);
	}

	FunctionType ConstantFunction::GetType() const {
		return FunctionType::CONSTANT;
	}

	Type ConstantFunction::GetValue() {
		return m_value;
	}

	// --- 
	// Polynomial Functions
	// ---

	PolynomialFunction::PolynomialFunction(
		DegreeType degree,
		std::vector<std::shared_ptr<Function>> coefficients
	) :
		Function(),
		m_degree(degree),
		m_coefficients(coefficients)
	{
		if(m_degree != m_coefficients.size() - 1) {
			throw std::runtime_error("Invalid polynomial degree.");
		}
	}

	Type PolynomialFunction::Evaluate(Type x) {
		Type result = 0;
		DegreeType currentDegree = m_degree;

		for(const auto& coefficient : m_coefficients) {
			result += coefficient->Evaluate(x) * std::pow(x, currentDegree--);
		}
		return result;
	}

	// Add inner derivative
	std::shared_ptr<Function> PolynomialFunction::GetDerivative() const {
		NOIMP;

		/*std::array<std::shared_ptr<Function>, Degree> derivedCoefficients;

		for(DegreeType i = 0; i < Degree; i++) {
			derivedCoefficients.at(i) = coefficients.at(i) * (Degree - i);
		}

		return std::make_shared<PolynomialFunction<Degree - 1>>(derivedCoefficients);
*/		}

	std::shared_ptr<Function> PolynomialFunction::GetPrimitive() const {
		NOIMP;

		// std::array<Type, Degree + 2> primitiveCoefficients;

		// for(DegreeType i = 0; i <= Degree; i++) {
		// 	primitiveCoefficients.at(i) = coefficients.at(i) / (Degree - i + 1);
		// }

		// return std::make_shared<PolynomialFunction<Degree + 1>>(primitiveCoefficients);
	}

	std::string PolynomialFunction::Repr() const {

		std::string repr = "";
		DegreeType currentDegree = m_degree;
		for(const auto& coefficient : m_coefficients) {

			std::string xForm = "x^" + std::to_string(currentDegree);
			
			if(currentDegree == 0) {
				xForm = "";
			} else if(currentDegree == 1) {
				xForm = "x";
			}

			repr +=
				std::string(currentDegree != m_degree ? "+ " : "")
				+ "(" + coefficient->Repr() + ") " + xForm
				+ std::string(currentDegree != 0 ? " " : ""); 

			currentDegree--;
		}

		return repr;
	}

	FunctionType PolynomialFunction::GetType() const {
		return FunctionType::POLYNOMIAL;
	}

	DegreeType PolynomialFunction::GetDegree() {
		return m_degree;
	}

	std::vector<std::shared_ptr<Function>> PolynomialFunction::GetCoefficients() {
		return m_coefficients;
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
// Operators Prototypes
// ---

// +

std::shared_ptr<functools::ConstantFunction> operator+(
	std::shared_ptr<functools::ConstantFunction> lhs,
	Type rhs
);
std::shared_ptr<functools::ConstantFunction> operator+(
	std::shared_ptr<functools::ConstantFunction> lhs,
	std::shared_ptr<functools::ConstantFunction> rhs
);
std::shared_ptr<functools::PolynomialFunction> operator+(
	std::shared_ptr<functools::ConstantFunction> lhs,
	std::shared_ptr<functools::PolynomialFunction> rhs
);
std::shared_ptr<functools::PolynomialFunction> operator+(
	std::shared_ptr<functools::PolynomialFunction> lhs,
	Type rhs
);
std::shared_ptr<functools::PolynomialFunction> operator+(
	std::shared_ptr<functools::PolynomialFunction> lhs,
	std::shared_ptr<functools::ConstantFunction> rhs
);
std::shared_ptr<functools::PolynomialFunction> operator+(
	std::shared_ptr<functools::PolynomialFunction> lhs,
	std::shared_ptr<functools::PolynomialFunction> rhs
);

// -

std::shared_ptr<functools::ConstantFunction> operator-(
	std::shared_ptr<functools::ConstantFunction> lhs,
	Type rhs
);
std::shared_ptr<functools::ConstantFunction> operator-(
	std::shared_ptr<functools::ConstantFunction> lhs,
	std::shared_ptr<functools::ConstantFunction> rhs
);
std::shared_ptr<functools::PolynomialFunction> operator-(
	std::shared_ptr<functools::ConstantFunction> lhs,
	std::shared_ptr<functools::PolynomialFunction> rhs
);
std::shared_ptr<functools::PolynomialFunction> operator-(
	std::shared_ptr<functools::PolynomialFunction> lhs,
	Type rhs
);
std::shared_ptr<functools::PolynomialFunction> operator-(
	std::shared_ptr<functools::PolynomialFunction> lhs,
	std::shared_ptr<functools::ConstantFunction> rhs
);
std::shared_ptr<functools::PolynomialFunction> operator-(
	std::shared_ptr<functools::PolynomialFunction> lhs,
	std::shared_ptr<functools::PolynomialFunction> rhs
);

// *

std::shared_ptr<functools::ConstantFunction> operator*(
	std::shared_ptr<functools::ConstantFunction> lhs,
	Type rhs
);
std::shared_ptr<functools::ConstantFunction> operator*(
	std::shared_ptr<functools::ConstantFunction> lhs,
	std::shared_ptr<functools::ConstantFunction> rhs
);
std::shared_ptr<functools::PolynomialFunction> operator*(
	std::shared_ptr<functools::ConstantFunction> lhs,
	std::shared_ptr<functools::PolynomialFunction> rhs
);
std::shared_ptr<functools::PolynomialFunction> operator*(
	std::shared_ptr<functools::PolynomialFunction> lhs,
	Type rhs
);
std::shared_ptr<functools::PolynomialFunction> operator*(
	std::shared_ptr<functools::PolynomialFunction> lhs,
	std::shared_ptr<functools::ConstantFunction> rhs
);
std::shared_ptr<functools::PolynomialFunction> operator*(
	std::shared_ptr<functools::PolynomialFunction> lhs,
	std::shared_ptr<functools::PolynomialFunction> rhs
);

// /

std::shared_ptr<functools::ConstantFunction> operator/(
	std::shared_ptr<functools::ConstantFunction> lhs,
	Type rhs
);
std::shared_ptr<functools::ConstantFunction> operator/(
	std::shared_ptr<functools::ConstantFunction> lhs,
	std::shared_ptr<functools::ConstantFunction> rhs
);
// Makes no sense
std::shared_ptr<functools::PolynomialFunction> operator/(
	std::shared_ptr<functools::ConstantFunction> lhs,
	std::shared_ptr<functools::PolynomialFunction> rhs
);
std::shared_ptr<functools::PolynomialFunction> operator/(
	std::shared_ptr<functools::PolynomialFunction> lhs,
	Type rhs
);
std::shared_ptr<functools::PolynomialFunction> operator/(
	std::shared_ptr<functools::PolynomialFunction> lhs,
	std::shared_ptr<functools::ConstantFunction> rhs
);
// Change return type (P(X) = D(X)Q(X) + R(X))
std::shared_ptr<functools::PolynomialFunction> operator/(
	std::shared_ptr<functools::PolynomialFunction> lhs,
	std::shared_ptr<functools::PolynomialFunction> rhs
);

// ---
// Operators
// ---

// +

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
std::shared_ptr<functools::PolynomialFunction> operator+(
	std::shared_ptr<functools::ConstantFunction> lhs,
	std::shared_ptr<functools::PolynomialFunction> rhs
) {
	DegreeType degree = rhs->GetDegree();
	std::vector<std::shared_ptr<functools::Function>> res = rhs->GetCoefficients();
	res.at(degree) = res.at(degree) + lhs->GetValue();
	return std::make_shared<functools::PolynomialFunction>(
		degree,
		res
	);
}
std::shared_ptr<functools::PolynomialFunction> operator+(
	std::shared_ptr<functools::PolynomialFunction> lhs,
	Type rhs
) {
	DegreeType degree = lhs->GetDegree();
	std::vector<std::shared_ptr<functools::Function>> res = lhs->GetCoefficients();
	res.at(degree) = res.at(degree) + rhs;
	return std::make_shared<functools::PolynomialFunction>(
		degree,
		res
	);
}
std::shared_ptr<functools::PolynomialFunction> operator+(
	std::shared_ptr<functools::PolynomialFunction> lhs,
	std::shared_ptr<functools::ConstantFunction> rhs
) {
	return rhs + lhs;
}
std::shared_ptr<functools::PolynomialFunction> operator+(
	std::shared_ptr<functools::PolynomialFunction> lhs,
	std::shared_ptr<functools::PolynomialFunction> rhs
) {

	DegreeType lhsDegree = lhs->GetDegree();
	DegreeType rhsDegree = rhs->GetDegree();

	if(lhsDegree == rhsDegree) {

		std::vector<std::shared_ptr<functools::Function>> res = lhs->GetCoefficients() + rhs->GetCoefficients();

		return std::make_shared<functools::PolynomialFunction>(
			lhsDegree,
			res
		);
	}

	std::vector<std::shared_ptr<functools::Function>> lhsCoeffs = lhs->GetCoefficients();
	std::vector<std::shared_ptr<functools::Function>> rhsCoeffs = rhs->GetCoefficients();

	if(lhsDegree > rhsDegree) {
		std::vector<std::shared_ptr<functools::Function>> res;
		for(auto func : lhsCoeffs) {
			res.push_back(func);
		}

		for(DegreeType i = 0; i <= rhsDegree; i++) {
			res.at(lhsDegree - rhsDegree + i) = res.at(lhsDegree - rhsDegree + i) + rhsCoeffs.at(i);
		}

		return std::make_shared<functools::PolynomialFunction>(
			lhsDegree,
			res
		);	
	}

	std::vector<std::shared_ptr<functools::Function>> res;
	for(auto func : rhsCoeffs) {
		res.push_back(func);
	}

	for(DegreeType i = 0; i <= lhsDegree; i++) {
		res.at(rhsDegree - lhsDegree + i) = res.at(rhsDegree - lhsDegree + i) + lhsCoeffs.at(i);
	}

	return std::make_shared<functools::PolynomialFunction>(
		rhsDegree,
		res
	);
}

// -

std::shared_ptr<functools::ConstantFunction> operator-(
	std::shared_ptr<functools::ConstantFunction> lhs,
	Type rhs
) {
	return lhs + (-rhs);
}
std::shared_ptr<functools::ConstantFunction> operator-(
	std::shared_ptr<functools::ConstantFunction> lhs,
	std::shared_ptr<functools::ConstantFunction> rhs
) {
	return lhs + rhs * (-1);
}
std::shared_ptr<functools::PolynomialFunction> operator-(
	std::shared_ptr<functools::ConstantFunction> lhs,
	std::shared_ptr<functools::PolynomialFunction> rhs
) {
	return lhs + rhs * (-1);
}
std::shared_ptr<functools::PolynomialFunction> operator-(
	std::shared_ptr<functools::PolynomialFunction> lhs,
	Type rhs
) {
	return lhs + (-rhs);
}
std::shared_ptr<functools::PolynomialFunction> operator-(
	std::shared_ptr<functools::PolynomialFunction> lhs,
	std::shared_ptr<functools::ConstantFunction> rhs
) {
	return lhs + rhs * (-1);
}
std::shared_ptr<functools::PolynomialFunction> operator-(
	std::shared_ptr<functools::PolynomialFunction> lhs,
	std::shared_ptr<functools::PolynomialFunction> rhs
) {
	return lhs + rhs * (-1);
}

// *

std::shared_ptr<functools::ConstantFunction> operator*(
	std::shared_ptr<functools::ConstantFunction> lhs,
	Type rhs
){
	return std::make_shared<functools::ConstantFunction>(
		lhs->GetValue() * rhs
	);	
}
std::shared_ptr<functools::ConstantFunction> operator*(
	std::shared_ptr<functools::ConstantFunction> lhs,
	std::shared_ptr<functools::ConstantFunction> rhs
){
	return std::make_shared<functools::ConstantFunction>(
		lhs->GetValue() * rhs->GetValue()
	);
}
std::shared_ptr<functools::PolynomialFunction> operator*(
	std::shared_ptr<functools::ConstantFunction> lhs,
	std::shared_ptr<functools::PolynomialFunction> rhs
){
	std::vector<std::shared_ptr<functools::Function>> res;
	for(auto func : rhs->GetCoefficients()) {
		res.push_back(func * lhs);
	}
	return std::make_shared<functools::PolynomialFunction>(
		rhs->GetDegree(),
		res
	);
}
std::shared_ptr<functools::PolynomialFunction> operator*(
	std::shared_ptr<functools::PolynomialFunction> lhs,
	Type rhs
){
	std::vector<std::shared_ptr<functools::Function>> res;
	for(auto func : lhs->GetCoefficients()) {
		res.push_back(func * rhs);
	}
	return std::make_shared<functools::PolynomialFunction>(
		lhs->GetDegree(),
		res
	);
}
std::shared_ptr<functools::PolynomialFunction> operator*(
	std::shared_ptr<functools::PolynomialFunction> lhs,
	std::shared_ptr<functools::ConstantFunction> rhs
){
	return rhs * lhs;
}
std::shared_ptr<functools::PolynomialFunction> operator*(
	std::shared_ptr<functools::PolynomialFunction> lhs,
	std::shared_ptr<functools::PolynomialFunction> rhs
){
	std::vector<std::shared_ptr<functools::Function>> res;
	for(uint8_t i = 0; i < lhs->GetDegree() + rhs->GetDegree() + 1; i++) {
		res.push_back(
			std::make_shared<functools::ConstantFunction>(0)
		);
	}
	std::vector<std::shared_ptr<functools::Function>> lhsCoeffs = lhs->GetCoefficients();
	std::vector<std::shared_ptr<functools::Function>> rhsCoeffs = rhs->GetCoefficients();

	for(DegreeType i = 0; i < lhsCoeffs.size(); i++) {
		for(DegreeType j = 0; j < rhsCoeffs.size(); j++) {
			res.at(i + j) = res.at(i + j) + lhsCoeffs.at(i) * rhsCoeffs.at(j);
		}
	}

	return std::make_shared<functools::PolynomialFunction>(
		lhs->GetDegree() + rhs->GetDegree(),
		res
	);
}

// /

std::shared_ptr<functools::ConstantFunction> operator/(
	std::shared_ptr<functools::ConstantFunction> lhs,
	Type rhs
){
	NOIMP;
}
std::shared_ptr<functools::ConstantFunction> operator/(
	std::shared_ptr<functools::ConstantFunction> lhs,
	std::shared_ptr<functools::ConstantFunction> rhs
){
	NOIMP;
}
// Makes no sense
std::shared_ptr<functools::PolynomialFunction> operator/(
	std::shared_ptr<functools::ConstantFunction> lhs,
	std::shared_ptr<functools::PolynomialFunction> rhs
){
	NOIMP;
}
std::shared_ptr<functools::PolynomialFunction> operator/(
	std::shared_ptr<functools::PolynomialFunction> lhs,
	Type rhs
){
	NOIMP;
}
std::shared_ptr<functools::PolynomialFunction> operator/(
	std::shared_ptr<functools::PolynomialFunction> lhs,
	std::shared_ptr<functools::ConstantFunction> rhs
){
	NOIMP;
}
// Change return type (P(X) = D(X)Q(X) + R(X))
std::shared_ptr<functools::PolynomialFunction> operator/(
	std::shared_ptr<functools::PolynomialFunction> lhs,
	std::shared_ptr<functools::PolynomialFunction> rhs
){
	NOIMP;
}

// ---
// General
// ---

std::shared_ptr<functools::Function> operator+(
	std::shared_ptr<functools::Function> lhs,
	std::shared_ptr<functools::Function> rhs
) {
	if(auto lhsCast = std::dynamic_pointer_cast<functools::ConstantFunction>(lhs)) {
		if(auto rhsCast = std::dynamic_pointer_cast<functools::ConstantFunction>(rhs)) {
			return lhsCast + rhsCast;
		}
		if(auto rhsCast = std::dynamic_pointer_cast<functools::PolynomialFunction>(rhs)) {
			return lhsCast + rhsCast;
		}
	}
	if(auto lhsCast = std::dynamic_pointer_cast<functools::PolynomialFunction>(lhs)) {
		if(auto rhsCast = std::dynamic_pointer_cast<functools::ConstantFunction>(rhs)) {
			return lhsCast + rhsCast;
		}
		if(auto rhsCast = std::dynamic_pointer_cast<functools::PolynomialFunction>(rhs)) {
			return lhsCast + rhsCast;
		}
	}
	
	throw std::runtime_error("Invalid function");
}

std::shared_ptr<functools::Function> operator+(
	std::shared_ptr<functools::Function> lhs,
	Type rhs
) {
	if(auto lhsCast = std::dynamic_pointer_cast<functools::ConstantFunction>(lhs)) {
		return lhsCast + rhs;
	}
	if(auto lhsCast = std::dynamic_pointer_cast<functools::PolynomialFunction>(lhs)) {
		return lhsCast + rhs;
	}

	throw std::runtime_error("Invalid function");
}

std::shared_ptr<functools::Function> operator*(
	std::shared_ptr<functools::Function> lhs,
	std::shared_ptr<functools::Function> rhs
) {
	if(auto lhsCast = std::dynamic_pointer_cast<functools::ConstantFunction>(lhs)) {
		if(auto rhsCast = std::dynamic_pointer_cast<functools::ConstantFunction>(rhs)) {
			return lhsCast * rhsCast;
		}
		if(auto rhsCast = std::dynamic_pointer_cast<functools::PolynomialFunction>(rhs)) {
			return lhsCast * rhsCast;
		}
	}
	if(auto lhsCast = std::dynamic_pointer_cast<functools::PolynomialFunction>(lhs)) {
		if(auto rhsCast = std::dynamic_pointer_cast<functools::ConstantFunction>(rhs)) {
			return lhsCast * rhsCast;
		}
		if(auto rhsCast = std::dynamic_pointer_cast<functools::PolynomialFunction>(rhs)) {
			return lhsCast * rhsCast;
		}
	}
	
	throw std::runtime_error("Invalid function");
}

std::shared_ptr<functools::Function> operator*(
	std::shared_ptr<functools::Function> lhs,
	Type rhs
) {
	if(auto lhsCast = std::dynamic_pointer_cast<functools::ConstantFunction>(lhs)) {
		return lhsCast * rhs;
	}
	if(auto lhsCast = std::dynamic_pointer_cast<functools::PolynomialFunction>(lhs)) {
		return lhsCast * rhs;
	}

	throw std::runtime_error("Invalid function");
}
