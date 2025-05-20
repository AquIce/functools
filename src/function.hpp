/*
 * TODO :
 * - Add operators to ComplexFunction (from -)
 * - Replace PowerN with ComplexFunction
 * - Switch general operators to have macro definition (reduces drastically the number of lines)
 * - Add comparison checks to functions (for ComplexFunction::isZero())
*/

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
		TRIGONOMETRY,
		COMPLEX
	};

	enum class FunctionOperator {
		PLUS,
		MINUS,
		TIMES,
		DIVIDED,
		POWER
	};

	enum class TrigonometryFunctionType {
		SIN = 0,
		CSC,
		COS,
		SEC,
		TAN,
		COT,
		SINH,
		CSCH,
		COSH,
		SECH,
		TANH,
		COTH
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
		virtual bool isZero() const = 0;
	};

	class ConstantFunction : public Function {
	public:
		ConstantFunction(Type value);

		Type Evaluate(Type x) override;

		std::shared_ptr<Function> GetDerivative() const override;

		std::shared_ptr<Function> GetPrimitive() const override;

		std::string Repr() const override;

		FunctionType GetType() const override;

		bool isZero() const override;

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

		PolynomialFunction();

		Type Evaluate(Type x) override;

		std::shared_ptr<Function> GetDerivative() const override;

		std::shared_ptr<Function> GetPrimitive() const override;

		std::string Repr() const override;

		FunctionType GetType() const override;

		bool isZero() const override;

		DegreeType GetDegree();

		std::vector<std::shared_ptr<Function>> GetCoefficients();

	private:
		const DegreeType m_degree;
		std::vector<std::shared_ptr<Function>> m_coefficients;
	};

	class TrigonometryFunction : public Function {
	
	public:

		TrigonometryFunction(
			const TrigonometryFunctionType type,
			std::shared_ptr<Function> inner
		);

		Type Evaluate(Type x) override;

		std::shared_ptr<Function> GetDerivative() const override;

		std::shared_ptr<Function> GetPrimitive() const override;

		std::string Repr() const override;

		FunctionType GetType() const override;

		bool isZero() const override;

	private:
		TrigonometryFunctionType m_type;
		std::shared_ptr<Function> m_inner;
	};

	class ComplexFunction : public Function {

	public:

		ComplexFunction(
			std::shared_ptr<Function> lhs,
			const FunctionOperator op,
			std::shared_ptr<Function> rhs
		);

		Type Evaluate(Type x) override;

		std::shared_ptr<Function> GetDerivative() const override;

		std::shared_ptr<Function> GetPrimitive() const override;

		std::string Repr() const override;

		FunctionType GetType() const override;

		bool isZero() const override;

	private:
		std::shared_ptr<Function> m_lhs;		
		std::shared_ptr<Function> m_rhs;
		FunctionOperator op;		
	};

	struct DivisionResult {
		std::shared_ptr<functools::Function> Quotient;
		std::shared_ptr<functools::Function> Remainder;
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

std::shared_ptr<functools::Function> operator-(
	std::shared_ptr<functools::Function> lhs,
	std::shared_ptr<functools::Function> rhs
);
std::shared_ptr<functools::Function> operator-(
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

struct functools::DivisionResult operator/(
	std::shared_ptr<functools::Function> lhs,
	std::shared_ptr<functools::Function> rhs
);
struct functools::DivisionResult operator/(
	std::shared_ptr<functools::Function> lhs,
	Type rhs
);

std::shared_ptr<functools::Function> PowerN(
	std::shared_ptr<functools::Function> func,
	Type exponent
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

	std::shared_ptr<PolynomialFunction> Simplify(
		std::shared_ptr<PolynomialFunction> func
	) {
		DegreeType i = 0;
		std::vector<std::shared_ptr<Function>> coeffs = func->GetCoefficients();
		std::vector<std::shared_ptr<Function>> res;
		while(i < coeffs.size() && coeffs.at(i)->isZero()) {
			i++;
		}
		for(; i < coeffs.size(); i++) {
			res.push_back(coeffs.at(i));
		}
		if(res.size() == 0) {
			res.push_back(
				std::make_shared<ConstantFunction>(0)
			);
		}
		return std::make_shared<PolynomialFunction>(
			res.size() - 1,
			res
		);
	}

	std::shared_ptr<PolynomialFunction> XPowerN(Type n) {
		auto res = std::vector<std::shared_ptr<Function>>(
			n + 1,
			std::dynamic_pointer_cast<Function>(
				std::make_shared<ConstantFunction>(0)
			)
		);
		res.at(0) = std::dynamic_pointer_cast<Function>(
			std::make_shared<ConstantFunction>(1)
		);
		return std::make_shared<PolynomialFunction>(n, res);
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

	bool ConstantFunction::isZero() const {
		return m_value == 0;
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

	PolynomialFunction::PolynomialFunction() :
		Function(),
		m_degree(0),
		m_coefficients(
			functools::Upcast<functools::ConstantFunction>(
				functools::CoeffsToConstFunctions(
					std::vector<Type>({0})
				)
			)
		)
	{}

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
		auto derivative = std::make_shared<PolynomialFunction>();

		for(DegreeType i = 0; i < m_degree; i++) {
			if(auto coefficientCast = std::dynamic_pointer_cast<ConstantFunction>(m_coefficients.at(i))) {
				derivative = std::dynamic_pointer_cast<PolynomialFunction>(
					derivative + (m_coefficients.at(i) * (m_degree - i) * XPowerN(m_degree - i - 1))
				);
				continue;
			}

			derivative = std::dynamic_pointer_cast<PolynomialFunction>(
				derivative + (
					(XPowerN(m_degree - i) * m_coefficients.at(i)->GetDerivative()) +
					(XPowerN(m_degree - i - 1) * m_coefficients.at(i))
				)
			);
		}

		return derivative;
	}

	std::shared_ptr<Function> PolynomialFunction::GetPrimitive() const {
		auto primitive = std::make_shared<PolynomialFunction>();

		for(DegreeType i = 0; i < m_degree + 1; i++) {
			if(auto coefficientCast = std::dynamic_pointer_cast<ConstantFunction>(m_coefficients.at(i))) {
				primitive = std::dynamic_pointer_cast<PolynomialFunction>(
					primitive + (
						(m_coefficients.at(i) / (m_degree - i + 1)).Quotient * XPowerN(m_degree - i + 1)
					)
				);
				continue;
			}

			if(auto coefficientCast = std::dynamic_pointer_cast<PolynomialFunction>(m_coefficients.at(i))) {
				primitive = std::dynamic_pointer_cast<PolynomialFunction>(
					primitive + (
						m_coefficients.at(i) * XPowerN(m_degree - i)
					)->GetPrimitive()
				);
				continue;
			} 

			// TODO Use integration by parts
			NOIMP;
		}

		return primitive;
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

	bool PolynomialFunction::isZero() const {
		for(const auto& coefficient : m_coefficients) {
			if(!coefficient->isZero()) {
				return false;
			}
		}
		return true;
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
	
	TrigonometryFunction::TrigonometryFunction(
		const TrigonometryFunctionType type,
		std::shared_ptr<Function> inner
	) :
		m_type(type),
		m_inner(inner)
	{}

	Type TrigonometryFunction::Evaluate(Type x) {
		switch(m_type) {
			case TrigonometryFunctionType::SIN:
				return static_cast<Type>(std::sin(x));
			case TrigonometryFunctionType::CSC:
				return static_cast<Type>(1 / std::sin(x));
			case TrigonometryFunctionType::COS:
				return static_cast<Type>(std::cos(x));
			case TrigonometryFunctionType::SEC:
				return static_cast<Type>(1 / std::cos(x));
			case TrigonometryFunctionType::TAN:
				return static_cast<Type>(std::tan(x));
			case TrigonometryFunctionType::COT:
				return static_cast<Type>(1 / std::tan(x));
			case TrigonometryFunctionType::SINH:
				return static_cast<Type>(std::sinh(x));
			case TrigonometryFunctionType::CSCH:
				return static_cast<Type>(1 / std::sinh(x));
			case TrigonometryFunctionType::COSH:
				return static_cast<Type>(std::cosh(x));
			case TrigonometryFunctionType::SECH:
				return static_cast<Type>(1 / std::cosh(x));
			case TrigonometryFunctionType::TANH:
				return static_cast<Type>(std::tanh(x));
			case TrigonometryFunctionType::COTH:
				return static_cast<Type>(1 / std::tanh(x));
			default:
				throw std::runtime_error("Invalid trigonometry operation");
		}
	}

	std::shared_ptr<Function> TrigonometryFunction::GetDerivative() const {
		switch(m_type) {
			case TrigonometryFunctionType::SIN: {
				return (
					std::make_shared<functools::TrigonometryFunction>(
						functools::TrigonometryFunctionType::COS,
						m_inner
					) * m_inner->GetDerivative()
				);
			}
			case TrigonometryFunctionType::CSC: {
				return (
					std::make_shared<functools::TrigonometryFunction>(
						functools::TrigonometryFunctionType::CSC,
						m_inner
					) * std::make_shared<functools::TrigonometryFunction>(
						functools::TrigonometryFunctionType::COT,
						m_inner
					) * m_inner->GetDerivative() * (-1)
				);
			}
			case TrigonometryFunctionType::COS: {
				return (
					std::make_shared<functools::TrigonometryFunction>(
						functools::TrigonometryFunctionType::SIN,
						m_inner
					) * m_inner->GetDerivative() * (-1)
				);
			}
			case TrigonometryFunctionType::SEC: {
				return (
					std::make_shared<functools::TrigonometryFunction>(
						functools::TrigonometryFunctionType::SEC,
						m_inner
					) * std::make_shared<functools::TrigonometryFunction>(
						functools::TrigonometryFunctionType::TAN,
						m_inner
					) * m_inner->GetDerivative()
				);
			}
			case TrigonometryFunctionType::TAN: {
				return (
					PowerN(
						std::make_shared<functools::TrigonometryFunction>(
							functools::TrigonometryFunctionType::SEC,
							m_inner
						),
						2
					) * m_inner->GetDerivative()
				);
			}
			case TrigonometryFunctionType::COT: {
				return (
					PowerN(
						std::make_shared<functools::TrigonometryFunction>(
							functools::TrigonometryFunctionType::CSC,
							m_inner
						),
						2
					) * m_inner->GetDerivative() * (-1)
				);
			}
			case TrigonometryFunctionType::SINH:
				return (
					std::make_shared<functools::TrigonometryFunction>(
						functools::TrigonometryFunctionType::COSH,
						m_inner
					) * m_inner->GetDerivative()
				);
			case TrigonometryFunctionType::CSCH:
				return (
					std::make_shared<functools::TrigonometryFunction>(
						functools::TrigonometryFunctionType::CSCH,
						m_inner
					) *
					std::make_shared<functools::TrigonometryFunction>(
						functools::TrigonometryFunctionType::TANH,
						m_inner
					) * m_inner->GetDerivative() * (-1)
				);
			case TrigonometryFunctionType::COSH:
				return (
					std::make_shared<functools::TrigonometryFunction>(
						functools::TrigonometryFunctionType::SINH,
						m_inner
					) * m_inner->GetDerivative()
				);
			case TrigonometryFunctionType::SECH:
				return (
					std::make_shared<functools::TrigonometryFunction>(
						functools::TrigonometryFunctionType::SECH,
						m_inner
					) *
					std::make_shared<functools::TrigonometryFunction>(
						functools::TrigonometryFunctionType::TANH,
						m_inner
					) * m_inner->GetDerivative() * (-1)
				);
			case TrigonometryFunctionType::TANH:
				return (
					PowerN(
						std::make_shared<functools::TrigonometryFunction>(
							functools::TrigonometryFunctionType::SECH,
							m_inner
						),
						2
					) * m_inner->GetDerivative()
				);
			case TrigonometryFunctionType::COTH:
				return (
					PowerN(
						std::make_shared<functools::TrigonometryFunction>(
							functools::TrigonometryFunctionType::COSH,
							m_inner
						),
						2
					) * m_inner->GetDerivative() * (-1)
				);
			default:
				throw std::runtime_error("Invalid function : " + this->Repr());
		}
	}

	std::shared_ptr<Function> TrigonometryFunction::GetPrimitive() const {
		NOIMP;
	}

	std::string TrigonometryFunction::Repr() const {
		std::string name =
			m_type == TrigonometryFunctionType::SIN ?
				"sin":
			m_type == TrigonometryFunctionType::COS ?
				"cos":
			m_type == TrigonometryFunctionType::TAN ?
				"tan":
			m_type == TrigonometryFunctionType::SINH ?
				"sinh":
			m_type == TrigonometryFunctionType::COSH ?
				"cosh":
				"tanh";

		return name + "(" + m_inner->Repr() + ")";
	}

	FunctionType TrigonometryFunction::GetType() const {
		return FunctionType::TRIGONOMETRY;
	}

	bool TrigonometryFunction::isZero() const {
		if(auto innerCast = std::dynamic_pointer_cast<ConstantFunction>(m_inner)) {
			switch(m_type) {
				case TrigonometryFunctionType::SIN: {
					return std::fmod(innerCast->GetValue(), 2 * M_PI) == 0;	
				}
				case TrigonometryFunctionType::COS: {
					return std::fmod(innerCast->GetValue(), 2 * M_PI) == M_PI;	
				}
				case TrigonometryFunctionType::TAN: {
					return std::fmod(innerCast->GetValue(), 2 * M_PI) == 0;	
				}
				case TrigonometryFunctionType::SINH: {
					return innerCast->GetValue() == 0;	
				}
				case TrigonometryFunctionType::COSH: {
					return false;	
				}
				case TrigonometryFunctionType::TANH: {
					return innerCast->GetValue() == 0;	
				}
				default:
					throw std::runtime_error("Invalid function : " + this->Repr());
			}
		}
		return false;
	}

	// --- 
	// Complex Functions
	// ---
	
	ComplexFunction::ComplexFunction(
		std::shared_ptr<Function> lhs,
		const FunctionOperator op,
		std::shared_ptr<Function> rhs
	) :
		m_lhs(lhs),
		m_op(op),
		m_rhs(rhs)
	{}

	Type ComplexFunction::Evaluate(Type x) {
		switch(m_op) {
			case FunctionOperator::PLUS:
				return m_lhs->Evaluate(x) + m_rhs->Evaluate(x);
			case FunctionOperator::MINUS:
				return m_lhs->Evaluate(x) + m_rhs->Evaluate(x);
			case FunctionOperator::TIMES:
				return m_lhs->Evaluate(x) + m_rhs->Evaluate(x);
			case FunctionOperator::DIVIDED:
				return m_lhs->Evaluate(x) + m_rhs->Evaluate(x);
			case FunctionOperator::POWER:
				return std::pow(m_lhs->Evaluate(x), m_rhs->Evaluate(x));
			default:
				throw std::runtime_error("Invalid operator");
		}
	}

	std::shared_ptr<Function> ComplexFunction::GetDerivative() const {
		switch(m_op) {
			case FunctionOperator::PLUS:
				return m_lhs->GetDerivative() + m_rhs->GetDerivative();
			case FunctionOperator::MINUS:
				return m_lhs->GetDerivative() - m_rhs->GetDerivative();
			case FunctionOperator::TIMES:
				return m_lhs->GetDerivative() * m_rhs + m_rhs->GetDerivative() * m_lhs;
			case FunctionOperator::DIVIDED:
				return m_lhs->GetDerivative() * m_rhs - m_rhs->GetDerivative() * m_lhs;
			case FunctionOperator::POWER:
				return m_rhs * std::make_shared<ComplexFunction>(
					m_lhs, m_op, m_rhs - 1
				) * m_rhs->GetDerivative();
			default:
				throw std::runtime_error("Invalid operator");
		}
	}

	std::shared_ptr<Function> ComplexFunction::GetPrimitive() const {
		switch(m_op) {
			case FunctionOperator::PLUS:
				return m_lhs->GetPrimitive() + m_rhs->GetPrimitive();
			case FunctionOperator::MINUS:
				return m_lhs->GetPrimitive() - m_rhs->GetPrimitive();
			case FunctionOperator::TIMES:
				NOIMP;
			case FunctionOperator::DIVIDED:
				NOIMP;
			case FunctionOperator::POWER:
				NOIMP;
			default:
				throw std::runtime_error("Invalid operator");
		}

	}

	std::string ComplexFunction::Repr() const {
		switch(m_op) {
			case FunctionOperator::PLUS:
				return m_lhs->Repr() + " + " + m_rhs->Repr();
			case FunctionOperator::MINUS:
				return m_lhs->Repr() + " - " + m_rhs->Repr();
			case FunctionOperator::TIMES:
				return m_lhs->Repr() + " * " + m_rhs->Repr();
			case FunctionOperator::DIVIDED:
				return m_lhs->Repr() + " / " + m_rhs->Repr();
			case FunctionOperator::POWER:
				return m_lhs->Repr() + " ^ " + m_rhs->Repr();
			default:
				throw std::runtime_error("Invalid operator");
		}
	}

	FunctionType ComplexFunction::GetType() const {
		return FunctionType::COMPLEX;
	}

	bool ComplexFunction::isZero() const {
		return m_lhs->isZero() && m_rhs->isZero();
	}
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
std::shared_ptr<functools::ComplexFunction> operator+(
	std::shared_ptr<functools::ConstantFunction> lhs,
	std::shared_ptr<functools::TrigonometryFunction> rhs
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
std::shared_ptr<functools::PolynomialFunction> operator+(
	std::shared_ptr<functools::PolynomialFunction> lhs,
	std::shared_ptr<functools::TrigonometryFunction> rhs
);

std::shared_ptr<functools::ComplexFunction> operator+(
	std::shared_ptr<functools::TrigonometryFunction> lhs,
	Type rhs
);
std::shared_ptr<functools::ComplexFunction> operator+(
	std::shared_ptr<functools::TrigonometryFunction> lhs,
	std::shared_ptr<functools::ConstantFunction> rhs
);
std::shared_ptr<functools::PolynomialFunction> operator+(
	std::shared_ptr<functools::TrigonometryFunction> lhs,
	std::shared_ptr<functools::PolynomialFunction> rhs
);
std::shared_ptr<functools::ComplexFunction> operator+(
	std::shared_ptr<functools::TrigonometryFunction> lhs,
	std::shared_ptr<functools::TrigonometryFunction> rhs
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
std::shared_ptr<functools::ComplexFunction> operator-(
	std::shared_ptr<functools::ConstantFunction> lhs,
	std::shared_ptr<functools::TrigonometryFunction> rhs
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
std::shared_ptr<functools::PolynomialFunction> operator-(
	std::shared_ptr<functools::PolynomialFunction> lhs,
	std::shared_ptr<functools::TrigonometryFunction> rhs
);

std::shared_ptr<functools::ComplexFunction> operator-(
	std::shared_ptr<functools::TrigonometryFunction> lhs,
	Type rhs
);
std::shared_ptr<functools::ComplexFunction> operator-(
	std::shared_ptr<functools::TrigonometryFunction> lhs,
	std::shared_ptr<functools::ConstantFunction> rhs
);
std::shared_ptr<functools::PolynomialFunction> operator-(
	std::shared_ptr<functools::TrigonometryFunction> lhs,
	std::shared_ptr<functools::PolynomialFunction> rhs
);
std::shared_ptr<functools::ComplexFunction> operator-(
	std::shared_ptr<functools::TrigonometryFunction> lhs,
	std::shared_ptr<functools::TrigonometryFunction> rhs
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
std::shared_ptr<functools::ComplexFunction> operator*(
	std::shared_ptr<functools::ConstantFunction> lhs,
	std::shared_ptr<functools::TrigonometryFunction> rhs
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
std::shared_ptr<functools::PolynomialFunction> operator*(
	std::shared_ptr<functools::PolynomialFunction> lhs,
	std::shared_ptr<functools::TrigonometryFunction> rhs
);

std::shared_ptr<functools::ComplexFunction> operator*(
	std::shared_ptr<functools::TrigonometryFunction> lhs,
	Type rhs
);
std::shared_ptr<functools::ComplexFunction> operator*(
	std::shared_ptr<functools::TrigonometryFunction> lhs,
	std::shared_ptr<functools::ConstantFunction> rhs
);
std::shared_ptr<functools::PolynomialFunction> operator*(
	std::shared_ptr<functools::TrigonometryFunction> lhs,
	std::shared_ptr<functools::PolynomialFunction> rhs
);
std::shared_ptr<functools::ComplexFunction> operator*(
	std::shared_ptr<functools::TrigonometryFunction> lhs,
	std::shared_ptr<functools::TrigonometryFunction> rhs
);

// /

struct functools::DivisionResult operator/(
	std::shared_ptr<functools::ConstantFunction> lhs,
	Type rhs
);
struct functools::DivisionResult operator/(
	std::shared_ptr<functools::ConstantFunction> lhs,
	std::shared_ptr<functools::ConstantFunction> rhs
);
struct functools::DivisionResult operator/(
	std::shared_ptr<functools::ConstantFunction> lhs,
	std::shared_ptr<functools::PolynomialFunction> rhs
);
struct functools::DivisionResult operator/(
	std::shared_ptr<functools::ConstantFunction> lhs,
	std::shared_ptr<functools::TrigonometryFunction> rhs
);

struct functools::DivisionResult operator/(
	std::shared_ptr<functools::PolynomialFunction> lhs,
	Type rhs
);
struct functools::DivisionResult operator/(
	std::shared_ptr<functools::PolynomialFunction> lhs,
	std::shared_ptr<functools::ConstantFunction> rhs
);
struct functools::DivisionResult operator/(
	std::shared_ptr<functools::PolynomialFunction> lhs,
	std::shared_ptr<functools::PolynomialFunction> rhs
);
struct functools::DivisionResult operator/(
	std::shared_ptr<functools::PolynomialFunction> lhs,
	std::shared_ptr<functools::TrigonometryFunction> rhs
);

struct functools::DivisionResult operator/(
	std::shared_ptr<functools::TrigonometryFunction> lhs,
	Type rhs
);
struct functools::DivisionResult operator/(
	std::shared_ptr<functools::TrigonometryFunction> lhs,
	std::shared_ptr<functools::ConstantFunction> rhs
);
struct functools::DivisionResult operator/(
	std::shared_ptr<functools::TrigonometryFunction> lhs,
	std::shared_ptr<functools::PolynomialFunction> rhs
);
struct functools::DivisionResult operator/(
	std::shared_ptr<functools::TrigonometryFunction> lhs,
	std::shared_ptr<functools::TrigonometryFunction> rhs
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
std::shared_ptr<functools::ComplexFunction> operator+(
	std::shared_ptr<functools::ConstantFunction> lhs,
	std::shared_ptr<functools::TrigonometryFunction> rhs
) {
	return std::make_shared<functools::ComplexFunction>(
		lhs,
		functools::FunctionOperator::PLUS,
		rhs
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
std::shared_ptr<functools::ComplexFunction> operator+(
	std::shared_ptr<functools::PolynomialFunction> lhs,
	std::shared_ptr<functools::TrigonometryFunction> rhs
) {
	return std::make_shared<functools::ComplexFunction>(
		lhs,
		functools::FunctionOperator::PLUS,
		rhs
	);
}

std::shared_ptr<functools::ComplexFunction> operator+(
	std::shared_ptr<functools::TrigonometryFunction> lhs,
	Type rhs
) {
	return std::make_shared<functools::ComplexFunction>(
		lhs,
		functools::FunctionOperator::PLUS,
		rhs
	);
}
std::shared_ptr<functools::TrigonometryFunction> operator+(
	std::shared_ptr<functools::TrigonometryFunction> lhs,
	std::shared_ptr<functools::ConstantFunction> rhs
) {
	return std::make_shared<functools::ComplexFunction>(
		lhs,
		functools::FunctionOperator::PLUS,
		rhs
	);
}
std::shared_ptr<functools::PolynomialFunction> operator+(
	std::shared_ptr<functools::TrigonometryFunction> lhs,
	std::shared_ptr<functools::PolynomialFunction> rhs
) {
	std::vector<std::shared_ptr<functools::Function>> res = rhs->GetCoefficients();

	res.back() = std::make_shared<functools::ComplexFunction>(
		lhs,
		functools::FunctionOperator::PLUS,
		res.back()
	);

	return std::make_shared<PolynomialFunction>(
		res.size() - 1,
		res
	);
}
std::shared_ptr<functools::ComplexFunction> operator+(
	std::shared_ptr<functools::TrigonometryFunction> lhs,
	std::shared_ptr<functools::TrigonometryFunction> rhs
) {
	return std::make_shared<functools::ComplexFunction>(
		lhs,
		functools::FunctionOperator::PLUS,
		rhs
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
	return lhs - rhs->GetValue();
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
	return lhs - rhs->GetValue();
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
	for(DegreeType i = 0; i < lhs->GetDegree() + rhs->GetDegree() + 1; i++) {
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

struct functools::DivisionResult operator/(
	std::shared_ptr<functools::ConstantFunction> lhs,
	Type rhs
){
	return functools::DivisionResult{
		std::make_shared<functools::ConstantFunction>(
			lhs->GetValue() / rhs
		),
		std::make_shared<functools::ConstantFunction>(0)
	};
}
struct functools::DivisionResult operator/(
	std::shared_ptr<functools::ConstantFunction> lhs,
	std::shared_ptr<functools::ConstantFunction> rhs
){
	return lhs / rhs->GetValue();
}
// Makes no sense (until 1/x-like functions are implemented)
struct functools::DivisionResult operator/(
	std::shared_ptr<functools::ConstantFunction> lhs,
	std::shared_ptr<functools::PolynomialFunction> rhs
){
	NOIMP;
}
struct functools::DivisionResult operator/(
	std::shared_ptr<functools::PolynomialFunction> lhs,
	Type rhs
){
	std::vector<std::shared_ptr<functools::Function>> res;
	for(auto func : lhs->GetCoefficients()) {
		res.push_back((func / rhs).Quotient);
	}
	return functools::DivisionResult{
		std::make_shared<functools::PolynomialFunction>(
			lhs->GetDegree(),
			res
		),
		std::make_shared<functools::ConstantFunction>(0)
	};
}
struct functools::DivisionResult operator/(
	std::shared_ptr<functools::PolynomialFunction> lhs,
	std::shared_ptr<functools::ConstantFunction> rhs
){
	return lhs / rhs->GetValue();
}

struct functools::DivisionResult operator/(
	std::shared_ptr<functools::PolynomialFunction> lhs,
	std::shared_ptr<functools::PolynomialFunction> rhs
){
	DegreeType lhsDegree = lhs->GetDegree();
	DegreeType rhsDegree = rhs->GetDegree();

	// Makes no sense (until 1/x-like functions are implemented)
	if(rhsDegree > lhsDegree) {
		NOIMP;
	}

	auto quotient = std::make_shared<functools::PolynomialFunction>();
	auto remainder = std::make_shared<functools::PolynomialFunction>(
		lhsDegree,
		lhs->GetCoefficients()
	);

	std::vector<std::shared_ptr<functools::Function>> rhsCoeffs = rhs->GetCoefficients();

	for(DegreeType i = 0; i < lhsDegree - rhsDegree + 1; i++) {

		functools::DivisionResult intermediaryResult = remainder->GetCoefficients().at(i) / rhsCoeffs.at(0);

		if(!intermediaryResult.Remainder->isZero()) {
			NOIMP;
		}

		auto intermediaryQuotientCoeffs = std::vector<std::shared_ptr<functools::Function>>(
			lhsDegree - rhsDegree - i + 1,
			std::dynamic_pointer_cast<functools::Function>(
				std::make_shared<functools::ConstantFunction>(0)
			)
		);

		intermediaryQuotientCoeffs.at(0) = intermediaryResult.Quotient;
		auto intermediaryQuotient = std::make_shared<functools::PolynomialFunction>(
			lhsDegree - rhsDegree - i,
			intermediaryQuotientCoeffs
		);

		quotient = std::dynamic_pointer_cast<functools::PolynomialFunction>(
			std::dynamic_pointer_cast<functools::Function>(quotient) + intermediaryQuotient
		);
		remainder = std::dynamic_pointer_cast<functools::PolynomialFunction>(
			std::dynamic_pointer_cast<functools::Function>(remainder) - intermediaryQuotient * rhs
		);
	}

	return functools::DivisionResult{
		functools::Simplify(quotient), functools::Simplify(remainder)
	};
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
		if(auto rhsCast = std::dynamic_pointer_cast<functools::TrigonometryFunction>(rhs)) {
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
		if(auto rhsCast = std::dynamic_pointer_cast<functools::TrigonometryFunction>(rhs)) {
			return lhsCast + rhsCast;
		}
	}
	if(auto lhsCast = std::dynamic_pointer_cast<functools::TrigonometryFunction>(lhs)) {
		if(auto rhsCast = std::dynamic_pointer_cast<functools::ConstantFunction>(rhs)) {
			return lhsCast + rhsCast;
		}
		if(auto rhsCast = std::dynamic_pointer_cast<functools::PolynomialFunction>(rhs)) {
			return lhsCast + rhsCast;
		}
		if(auto rhsCast = std::dynamic_pointer_cast<functools::TrigonometryFunction>(rhs)) {
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
	if(auto lhsCast = std::dynamic_pointer_cast<functools::TrigonometryFunction>(lhs)) {
		return lhsCast + rhs;
	}

	throw std::runtime_error("Invalid function");
}

std::shared_ptr<functools::Function> operator-(
	std::shared_ptr<functools::Function> lhs,
	std::shared_ptr<functools::Function> rhs
) {
	if(auto lhsCast = std::dynamic_pointer_cast<functools::ConstantFunction>(lhs)) {
		if(auto rhsCast = std::dynamic_pointer_cast<functools::ConstantFunction>(rhs)) {
			return lhsCast - rhsCast;
		}
		if(auto rhsCast = std::dynamic_pointer_cast<functools::PolynomialFunction>(rhs)) {
			return lhsCast - rhsCast;
		}
		if(auto lhsCast = std::dynamic_pointer_cast<functools::TrigonometryFunction>(lhs)) {
			return lhsCast - rhsCast;
		}
	}
	if(auto lhsCast = std::dynamic_pointer_cast<functools::PolynomialFunction>(lhs)) {
		if(auto rhsCast = std::dynamic_pointer_cast<functools::ConstantFunction>(rhs)) {
			return lhsCast - rhsCast;
		}
		if(auto rhsCast = std::dynamic_pointer_cast<functools::PolynomialFunction>(rhs)) {
			return lhsCast - rhsCast;
		}
		if(auto lhsCast = std::dynamic_pointer_cast<functools::TrigonometryFunction>(lhs)) {
			return lhsCast - rhsCast;
		}
	}
	if(auto lhsCast = std::dynamic_pointer_cast<functools::TrigonometryFunction>(lhs)) {
		if(auto rhsCast = std::dynamic_pointer_cast<functools::ConstantFunction>(rhs)) {
			return lhsCast - rhsCast;
		}
		if(auto rhsCast = std::dynamic_pointer_cast<functools::PolynomialFunction>(rhs)) {
			return lhsCast - rhsCast;
		}
		if(auto lhsCast = std::dynamic_pointer_cast<functools::TrigonometryFunction>(lhs)) {
			return lhsCast - rhsCast;
		}
	}
	
	throw std::runtime_error("Invalid function");
}

std::shared_ptr<functools::Function> operator-(
	std::shared_ptr<functools::Function> lhs,
	Type rhs
) {
	if(auto lhsCast = std::dynamic_pointer_cast<functools::ConstantFunction>(lhs)) {
		return lhsCast - rhs;
	}
	if(auto lhsCast = std::dynamic_pointer_cast<functools::PolynomialFunction>(lhs)) {
		return lhsCast - rhs;
	}
	if(auto lhsCast = std::dynamic_pointer_cast<functools::TrigonometryFunction>(lhs)) {
		return lhsCast - rhs;
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
		if(auto rhsCast = std::dynamic_pointer_cast<functools::TrigonometryFunction>(rhs)) {
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
		if(auto rhsCast = std::dynamic_pointer_cast<functools::TrigonometryFunction>(rhs)) {
			return lhsCast * rhsCast;
		}
	}
	if(auto lhsCast = std::dynamic_pointer_cast<functools::TrigonometryFunction>(lhs)) {
		if(auto rhsCast = std::dynamic_pointer_cast<functools::ConstantFunction>(rhs)) {
			return lhsCast * rhsCast;
		}
		if(auto rhsCast = std::dynamic_pointer_cast<functools::PolynomialFunction>(rhs)) {
			return lhsCast * rhsCast;
		}
		if(auto rhsCast = std::dynamic_pointer_cast<functools::TrigonometryFunction>(rhs)) {
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
	if(auto lhsCast = std::dynamic_pointer_cast<functools::TrigonometryFunction>(lhs)) {
		return lhsCast * rhs;
	}

	throw std::runtime_error("Invalid function");
}

struct functools::DivisionResult operator/(
	std::shared_ptr<functools::Function> lhs,
	std::shared_ptr<functools::Function> rhs
) {
	if(auto lhsCast = std::dynamic_pointer_cast<functools::ConstantFunction>(lhs)) {
		if(auto rhsCast = std::dynamic_pointer_cast<functools::ConstantFunction>(rhs)) {
			return lhsCast / rhsCast;
		}
		if(auto rhsCast = std::dynamic_pointer_cast<functools::PolynomialFunction>(rhs)) {
			return lhsCast / rhsCast;
		}
		if(auto rhsCast = std::dynamic_pointer_cast<functools::TrigonometryFunction>(rhs)) {
			return lhsCast / rhsCast;
		}
	}
	if(auto lhsCast = std::dynamic_pointer_cast<functools::PolynomialFunction>(lhs)) {
		if(auto rhsCast = std::dynamic_pointer_cast<functools::ConstantFunction>(rhs)) {
			return lhsCast / rhsCast;
		}
		if(auto rhsCast = std::dynamic_pointer_cast<functools::PolynomialFunction>(rhs)) {
			return lhsCast / rhsCast;
		}
		if(auto rhsCast = std::dynamic_pointer_cast<functools::TrigonometryFunction>(rhs)) {
			return lhsCast / rhsCast;
		}
	}
	if(auto lhsCast = std::dynamic_pointer_cast<functools::TrigonometryFunction>(lhs)) {
		if(auto rhsCast = std::dynamic_pointer_cast<functools::ConstantFunction>(rhs)) {
			return lhsCast / rhsCast;
		}
		if(auto rhsCast = std::dynamic_pointer_cast<functools::PolynomialFunction>(rhs)) {
			return lhsCast / rhsCast;
		}
		if(auto rhsCast = std::dynamic_pointer_cast<functools::TrigonometryFunction>(rhs)) {
			return lhsCast / rhsCast;
		}
	}
	
	throw std::runtime_error("Invalid function");
}

struct functools::DivisionResult operator/(
	std::shared_ptr<functools::Function> lhs,
	Type rhs
) {
	if(auto lhsCast = std::dynamic_pointer_cast<functools::ConstantFunction>(lhs)) {
		return lhsCast / rhs;
	}
	if(auto lhsCast = std::dynamic_pointer_cast<functools::PolynomialFunction>(lhs)) {
		return lhsCast / rhs;
	}
	if(auto lhsCast = std::dynamic_pointer_cast<functools::TrigonometryFunction>(lhs)) {
		return lhsCast / rhs;
	}

	throw std::runtime_error("Invalid function");
}

std::shared_ptr<functools::Function> PowerN(
	std::shared_ptr<functools::Function> func,
	Type exponent
) {
	std::shared_ptr<functools::Function> fn = func;
	for(Type i = 0; i < exponent; i++) {
		fn = fn * func;
	}
	return fn;
}