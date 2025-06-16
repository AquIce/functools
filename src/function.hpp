/*
 * TODO :
 * - Add auto simplification to (...) * (...)^2
 * - Switch general operators to have macro definition (reduces drastically the number of lines)
 * - Add comparison checks to functions (for ComplexFunction::isZero())
 * - Move operators logic to ComplexFunction (all operators only create a ComplexFunction)
 * - Add primitives to trigonometric functions
 * - Add primitives that need by part integration
 * - Add log_{n} functions
 * - Add exponential functions
 * - Add sgn function
 * - Add abs function
 * - Add static evaluation to any function containing only ConstantFunctions
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
		IDENTITY,
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
		virtual bool needsParentheses() const = 0;
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
		bool needsParentheses() const override;

		Type GetValue();

	private:
		Type m_value;
	};

	class IdentityFunction : public Function {
	
	public:

		IdentityFunction();

		Type Evaluate(Type x) override;

		std::shared_ptr<Function> GetDerivative() const override;
		std::shared_ptr<Function> GetPrimitive() const override;

		std::string Repr() const override;

		FunctionType GetType() const override;
		bool isZero() const override;
		bool needsParentheses() const override;
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
		bool needsParentheses() const override;

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
		bool needsParentheses() const override;

		std::shared_ptr<Function> Operate(
			std::shared_ptr<Function> other,
			FunctionOperator op
		);

		std::shared_ptr<Function> Simplify();

	private:
		std::shared_ptr<Function> m_lhs;		
		FunctionOperator m_op;		
		std::shared_ptr<Function> m_rhs;
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

std::shared_ptr<functools::Function> operator/(
	std::shared_ptr<functools::Function> lhs,
	std::shared_ptr<functools::Function> rhs
);
std::shared_ptr<functools::Function> operator/(
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
		return std::make_shared<ComplexFunction>(
			std::make_shared<ConstantFunction>(m_value),
			FunctionOperator::TIMES,
			std::make_shared<IdentityFunction>()
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

	bool ConstantFunction::needsParentheses() const {
		return false;
	}

	Type ConstantFunction::GetValue() {
		return m_value;
	}

	// --- 
	// Identity Functions
	// ---

	IdentityFunction::IdentityFunction() {}

	Type IdentityFunction::Evaluate(Type x) {
		return x;
	}

	std::shared_ptr<Function> IdentityFunction::GetDerivative() const {
		return std::make_shared<ConstantFunction>(1);
	}

	std::shared_ptr<Function> IdentityFunction::GetPrimitive() const {
		return std::make_shared<ComplexFunction>(
			std::make_shared<ConstantFunction>(1.0 / 2),
			FunctionOperator::TIMES,
			std::make_shared<ComplexFunction>(
				std::make_shared<IdentityFunction>(),
				FunctionOperator::POWER,
				std::make_shared<ConstantFunction>(2)
			)
		);
	}

	std::string IdentityFunction::Repr() const {
		return "x";
	}

	FunctionType IdentityFunction::GetType() const {
		return FunctionType::IDENTITY;
	}

	bool IdentityFunction::isZero() const {
		return false;
	}

	bool IdentityFunction::needsParentheses() const {
		return false;
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
					std::make_shared<TrigonometryFunction>(
						TrigonometryFunctionType::COS,
						m_inner
					) * m_inner->GetDerivative()
				);
			}
			case TrigonometryFunctionType::CSC: {
				return (
					std::make_shared<TrigonometryFunction>(
						TrigonometryFunctionType::CSC,
						m_inner
					) * std::make_shared<TrigonometryFunction>(
						TrigonometryFunctionType::COT,
						m_inner
					) * m_inner->GetDerivative() * (-1)
				);
			}
			case TrigonometryFunctionType::COS: {
				return (
					std::make_shared<TrigonometryFunction>(
						TrigonometryFunctionType::SIN,
						m_inner
					) * m_inner->GetDerivative() * (-1)
				);
			}
			case TrigonometryFunctionType::SEC: {
				return (
					std::make_shared<TrigonometryFunction>(
						TrigonometryFunctionType::SEC,
						m_inner
					) * std::make_shared<TrigonometryFunction>(
						TrigonometryFunctionType::TAN,
						m_inner
					) * m_inner->GetDerivative()
				);
			}
			case TrigonometryFunctionType::TAN: {
				return (
					std::make_shared<ComplexFunction>(
						std::make_shared<TrigonometryFunction>(
							TrigonometryFunctionType::SEC,
							m_inner
						),
						FunctionOperator::POWER,
						std::make_shared<ConstantFunction>(2)
					) * m_inner->GetDerivative()
				);
			}
			case TrigonometryFunctionType::COT: {
				return (
					std::make_shared<ComplexFunction>(
						std::make_shared<TrigonometryFunction>(
							TrigonometryFunctionType::CSC,
							m_inner
						),
						FunctionOperator::POWER,
						std::make_shared<ConstantFunction>(2)
					) * m_inner->GetDerivative() * (-1)
				);
			}
			case TrigonometryFunctionType::SINH:
				return (
					std::make_shared<TrigonometryFunction>(
						TrigonometryFunctionType::COSH,
						m_inner
					) * m_inner->GetDerivative()
				);
			case TrigonometryFunctionType::CSCH:
				return (
					std::make_shared<TrigonometryFunction>(
						TrigonometryFunctionType::CSCH,
						m_inner
					) *
					std::make_shared<TrigonometryFunction>(
						TrigonometryFunctionType::TANH,
						m_inner
					) * m_inner->GetDerivative() * (-1)
				);
			case TrigonometryFunctionType::COSH:
				return (
					std::make_shared<TrigonometryFunction>(
						TrigonometryFunctionType::SINH,
						m_inner
					) * m_inner->GetDerivative()
				);
			case TrigonometryFunctionType::SECH:
				return (
					std::make_shared<TrigonometryFunction>(
						TrigonometryFunctionType::SECH,
						m_inner
					) *
					std::make_shared<TrigonometryFunction>(
						TrigonometryFunctionType::TANH,
						m_inner
					) * m_inner->GetDerivative() * (-1)
				);
			case TrigonometryFunctionType::TANH:
				return (
					std::make_shared<ComplexFunction>(
						std::make_shared<TrigonometryFunction>(
							TrigonometryFunctionType::SECH,
							m_inner
						),
						FunctionOperator::POWER,
						std::make_shared<ConstantFunction>(2)
					) * m_inner->GetDerivative()
				);
			case TrigonometryFunctionType::COTH:
				return (
					std::make_shared<ComplexFunction>(
						std::make_shared<TrigonometryFunction>(
							TrigonometryFunctionType::COSH,
							m_inner
						),
						FunctionOperator::POWER,
						std::make_shared<ConstantFunction>(2)
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

	bool TrigonometryFunction::needsParentheses() const {
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
	{
		if(m_op == FunctionOperator::DIVIDED && m_rhs->isZero()) {
			throw std::runtime_error("Division by zero"); 
		}
		if(auto lhsCast = std::dynamic_pointer_cast<ComplexFunction>(m_lhs)) {
			if(auto lhsSimplified = lhsCast->Simplify()) {
				m_lhs = lhsSimplified;
			}
		}
		if(auto rhsCast = std::dynamic_pointer_cast<ComplexFunction>(m_rhs)) {
			if(auto rhsSimplified = rhsCast->Simplify()) {
				m_rhs = rhsSimplified;
			}
		}
	}

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
			case FunctionOperator::TIMES: {
				auto a = m_lhs->GetDerivative() * m_rhs;
				auto b = m_rhs->GetDerivative() * m_lhs;

				return m_lhs->GetDerivative() * m_rhs + m_rhs->GetDerivative() * m_lhs;
			}
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
		std::string lhsRepr = m_lhs->needsParentheses() ?
			(std::string("(") + m_lhs->Repr() + ")") :
			m_lhs->Repr();
		std::string rhsRepr = m_rhs->needsParentheses() ?
			(std::string("(") + m_rhs->Repr() + ")") :
			m_rhs->Repr();

		switch(m_op) {
			case FunctionOperator::PLUS:
				return lhsRepr + " + " + rhsRepr;
			case FunctionOperator::MINUS:
				return lhsRepr + " - " + rhsRepr;
			case FunctionOperator::TIMES:
				return lhsRepr + " * " + rhsRepr;
			case FunctionOperator::DIVIDED:
				return lhsRepr + " / " + rhsRepr;
			case FunctionOperator::POWER:
				return lhsRepr + " ^ " + rhsRepr;
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

	bool ComplexFunction::needsParentheses() const {
		return false;
	}

	std::shared_ptr<Function> ComplexFunction::Operate(
		std::shared_ptr<Function> other,
		FunctionOperator op
	) {
		// See TODO
		return nullptr;
	}

	std::shared_ptr<Function> ComplexFunction::Simplify() {
		switch(m_op) {
			case FunctionOperator::PLUS: {
				if(m_lhs->isZero()) {
					return m_rhs;
				}
				if(m_rhs->isZero()) {
					return m_lhs;
				}
				if(m_lhs->Repr() == m_rhs->Repr()) {
					return std::make_shared<ComplexFunction>(
						m_lhs,
						FunctionOperator::TIMES,
						std::make_shared<ConstantFunction>(2)
					);
				}
				break;
			}
			case FunctionOperator::MINUS: {
				if(m_lhs->isZero()) {
					return m_rhs * (-1);
				}
				if(m_rhs->isZero()) {
					return m_lhs;
				}
				if(m_lhs->Repr() == m_rhs->Repr()) {
					return std::make_shared<ConstantFunction>(0);
				}
				break;
			}
			case FunctionOperator::TIMES: {
				if(m_lhs->isZero() || m_rhs->isZero()) {
					return std::make_shared<ConstantFunction>(0);
				}
				if(auto lhsCast = std::dynamic_pointer_cast<ConstantFunction>(m_lhs)) {
					if(lhsCast->GetValue() == 1) {
						return m_rhs;
					}
				}
				if(auto rhsCast = std::dynamic_pointer_cast<ConstantFunction>(m_lhs)) {
					if(rhsCast->GetValue() == 1) {
						return m_lhs;
					}
				}
				if(m_lhs->Repr() == m_rhs->Repr()) {
					return std::make_shared<ComplexFunction>(
						m_lhs,
						FunctionOperator::POWER,
						std::make_shared<ConstantFunction>(2)
					);
				}
				break;
			}
			case FunctionOperator::DIVIDED: {
				if(m_lhs->isZero()) {
					return std::make_shared<ConstantFunction>(0);
				}
				if(auto rhsCast = std::dynamic_pointer_cast<ConstantFunction>(m_lhs)) {
					if(rhsCast->GetValue() == 1) {
						return m_lhs;
					}
				}
				if(m_lhs->Repr() == m_rhs->Repr()) {
					return std::make_shared<ConstantFunction>(1);
				}
				break;
			}
			case FunctionOperator::POWER: {
				if(m_lhs->isZero()) {
					return std::make_shared<ConstantFunction>(0);
				}
				if(m_rhs->isZero()) {
					return std::make_shared<ConstantFunction>(1);
				}
				if(auto lhsCast = std::dynamic_pointer_cast<ConstantFunction>(m_lhs)) {
					if(lhsCast->GetValue() == 1) {
						return std::make_shared<ConstantFunction>(1);
					}
				}
				if(auto rhsCast = std::dynamic_pointer_cast<ConstantFunction>(m_lhs)) {
					if(rhsCast->GetValue() == 1) {
						return m_lhs;
					}
				}
				break;
			}
		}

		return nullptr;
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
std::shared_ptr<functools::ComplexFunction> operator+(
	std::shared_ptr<functools::ConstantFunction> lhs,
	std::shared_ptr<functools::IdentityFunction> rhs
);
std::shared_ptr<functools::ComplexFunction> operator+(
	std::shared_ptr<functools::ConstantFunction> lhs,
	std::shared_ptr<functools::TrigonometryFunction> rhs
);
std::shared_ptr<functools::ComplexFunction> operator+(
	std::shared_ptr<functools::ConstantFunction> lhs,
	std::shared_ptr<functools::ComplexFunction> rhs
);

std::shared_ptr<functools::ComplexFunction> operator+(
	std::shared_ptr<functools::IdentityFunction> lhs,
	Type rhs
);
std::shared_ptr<functools::ComplexFunction> operator+(
	std::shared_ptr<functools::IdentityFunction> lhs,
	std::shared_ptr<functools::ConstantFunction> rhs
);
std::shared_ptr<functools::ComplexFunction> operator+(
	std::shared_ptr<functools::IdentityFunction> lhs,
	std::shared_ptr<functools::IdentityFunction> rhs
);
std::shared_ptr<functools::ComplexFunction> operator+(
	std::shared_ptr<functools::IdentityFunction> lhs,
	std::shared_ptr<functools::TrigonometryFunction> rhs
);
std::shared_ptr<functools::ComplexFunction> operator+(
	std::shared_ptr<functools::IdentityFunction> lhs,
	std::shared_ptr<functools::ComplexFunction> rhs
);

std::shared_ptr<functools::ComplexFunction> operator+(
	std::shared_ptr<functools::TrigonometryFunction> lhs,
	Type rhs
);
std::shared_ptr<functools::ComplexFunction> operator+(
	std::shared_ptr<functools::TrigonometryFunction> lhs,
	std::shared_ptr<functools::ConstantFunction> rhs
);
std::shared_ptr<functools::ComplexFunction> operator+(
	std::shared_ptr<functools::TrigonometryFunction> lhs,
	std::shared_ptr<functools::IdentityFunction> rhs
);
std::shared_ptr<functools::ComplexFunction> operator+(
	std::shared_ptr<functools::TrigonometryFunction> lhs,
	std::shared_ptr<functools::TrigonometryFunction> rhs
);
std::shared_ptr<functools::ComplexFunction> operator+(
	std::shared_ptr<functools::TrigonometryFunction> lhs,
	std::shared_ptr<functools::ComplexFunction> rhs
);

std::shared_ptr<functools::ComplexFunction> operator+(
	std::shared_ptr<functools::ComplexFunction> lhs,
	Type rhs
);
std::shared_ptr<functools::ComplexFunction> operator+(
	std::shared_ptr<functools::ComplexFunction> lhs,
	std::shared_ptr<functools::ConstantFunction> rhs
);
std::shared_ptr<functools::ComplexFunction> operator+(
	std::shared_ptr<functools::ComplexFunction> lhs,
	std::shared_ptr<functools::IdentityFunction> rhs
);
std::shared_ptr<functools::ComplexFunction> operator+(
	std::shared_ptr<functools::ComplexFunction> lhs,
	std::shared_ptr<functools::TrigonometryFunction> rhs
);
std::shared_ptr<functools::ComplexFunction> operator+(
	std::shared_ptr<functools::ComplexFunction> lhs,
	std::shared_ptr<functools::ComplexFunction> rhs
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
std::shared_ptr<functools::ComplexFunction> operator-(
	std::shared_ptr<functools::ConstantFunction> lhs,
	std::shared_ptr<functools::IdentityFunction> rhs
);
std::shared_ptr<functools::ComplexFunction> operator-(
	std::shared_ptr<functools::ConstantFunction> lhs,
	std::shared_ptr<functools::TrigonometryFunction> rhs
);
std::shared_ptr<functools::ComplexFunction> operator-(
	std::shared_ptr<functools::ConstantFunction> lhs,
	std::shared_ptr<functools::ComplexFunction> rhs
);

std::shared_ptr<functools::ComplexFunction> operator-(
	std::shared_ptr<functools::IdentityFunction> lhs,
	Type rhs
);
std::shared_ptr<functools::ComplexFunction> operator-(
	std::shared_ptr<functools::IdentityFunction> lhs,
	std::shared_ptr<functools::ConstantFunction> rhs
);
std::shared_ptr<functools::ConstantFunction> operator-(
	std::shared_ptr<functools::IdentityFunction> lhs,
	std::shared_ptr<functools::IdentityFunction> rhs
);
std::shared_ptr<functools::ComplexFunction> operator-(
	std::shared_ptr<functools::IdentityFunction> lhs,
	std::shared_ptr<functools::TrigonometryFunction> rhs
);
std::shared_ptr<functools::ComplexFunction> operator-(
	std::shared_ptr<functools::IdentityFunction> lhs,
	std::shared_ptr<functools::ComplexFunction> rhs
);

std::shared_ptr<functools::ComplexFunction> operator-(
	std::shared_ptr<functools::TrigonometryFunction> lhs,
	Type rhs
);
std::shared_ptr<functools::ComplexFunction> operator-(
	std::shared_ptr<functools::TrigonometryFunction> lhs,
	std::shared_ptr<functools::ConstantFunction> rhs
);
std::shared_ptr<functools::ComplexFunction> operator-(
	std::shared_ptr<functools::TrigonometryFunction> lhs,
	std::shared_ptr<functools::IdentityFunction> rhs
);
std::shared_ptr<functools::ComplexFunction> operator-(
	std::shared_ptr<functools::TrigonometryFunction> lhs,
	std::shared_ptr<functools::TrigonometryFunction> rhs
);
std::shared_ptr<functools::ComplexFunction> operator-(
	std::shared_ptr<functools::TrigonometryFunction> lhs,
	std::shared_ptr<functools::ComplexFunction> rhs
);

std::shared_ptr<functools::ComplexFunction> operator-(
	std::shared_ptr<functools::ComplexFunction> lhs,
	Type rhs
);
std::shared_ptr<functools::ComplexFunction> operator-(
	std::shared_ptr<functools::ComplexFunction> lhs,
	std::shared_ptr<functools::ConstantFunction> rhs
);
std::shared_ptr<functools::ComplexFunction> operator-(
	std::shared_ptr<functools::ComplexFunction> lhs,
	std::shared_ptr<functools::IdentityFunction> rhs
);
std::shared_ptr<functools::ComplexFunction> operator-(
	std::shared_ptr<functools::ComplexFunction> lhs,
	std::shared_ptr<functools::TrigonometryFunction> rhs
);
std::shared_ptr<functools::ComplexFunction> operator-(
	std::shared_ptr<functools::ComplexFunction> lhs,
	std::shared_ptr<functools::ComplexFunction> rhs
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
std::shared_ptr<functools::ComplexFunction> operator*(
	std::shared_ptr<functools::ConstantFunction> lhs,
	std::shared_ptr<functools::IdentityFunction> rhs
);
std::shared_ptr<functools::ComplexFunction> operator*(
	std::shared_ptr<functools::ConstantFunction> lhs,
	std::shared_ptr<functools::TrigonometryFunction> rhs
);
std::shared_ptr<functools::ComplexFunction> operator*(
	std::shared_ptr<functools::ConstantFunction> lhs,
	std::shared_ptr<functools::ComplexFunction> rhs
);

std::shared_ptr<functools::ComplexFunction> operator*(
	std::shared_ptr<functools::IdentityFunction> lhs,
	Type rhs
);
std::shared_ptr<functools::ComplexFunction> operator*(
	std::shared_ptr<functools::IdentityFunction> lhs,
	std::shared_ptr<functools::ConstantFunction> rhs
);
std::shared_ptr<functools::ComplexFunction> operator*(
	std::shared_ptr<functools::IdentityFunction> lhs,
	std::shared_ptr<functools::IdentityFunction> rhs
);
std::shared_ptr<functools::ComplexFunction> operator*(
	std::shared_ptr<functools::IdentityFunction> lhs,
	std::shared_ptr<functools::TrigonometryFunction> rhs
);
std::shared_ptr<functools::ComplexFunction> operator*(
	std::shared_ptr<functools::IdentityFunction> lhs,
	std::shared_ptr<functools::ComplexFunction> rhs
);

std::shared_ptr<functools::ComplexFunction> operator*(
	std::shared_ptr<functools::TrigonometryFunction> lhs,
	Type rhs
);
std::shared_ptr<functools::ComplexFunction> operator*(
	std::shared_ptr<functools::TrigonometryFunction> lhs,
	std::shared_ptr<functools::ConstantFunction> rhs
);
std::shared_ptr<functools::ComplexFunction> operator*(
	std::shared_ptr<functools::TrigonometryFunction> lhs,
	std::shared_ptr<functools::IdentityFunction> rhs
);
std::shared_ptr<functools::ComplexFunction> operator*(
	std::shared_ptr<functools::TrigonometryFunction> lhs,
	std::shared_ptr<functools::TrigonometryFunction> rhs
);
std::shared_ptr<functools::ComplexFunction> operator*(
	std::shared_ptr<functools::TrigonometryFunction> lhs,
	std::shared_ptr<functools::ComplexFunction> rhs
);

std::shared_ptr<functools::ComplexFunction> operator*(
	std::shared_ptr<functools::ComplexFunction> lhs,
	Type rhs
);
std::shared_ptr<functools::ComplexFunction> operator*(
	std::shared_ptr<functools::ComplexFunction> lhs,
	std::shared_ptr<functools::ConstantFunction> rhs
);
std::shared_ptr<functools::ComplexFunction> operator*(
	std::shared_ptr<functools::ComplexFunction> lhs,
	std::shared_ptr<functools::IdentityFunction> rhs
);
std::shared_ptr<functools::ComplexFunction> operator*(
	std::shared_ptr<functools::ComplexFunction> lhs,
	std::shared_ptr<functools::TrigonometryFunction> rhs
);
std::shared_ptr<functools::ComplexFunction> operator*(
	std::shared_ptr<functools::ComplexFunction> lhs,
	std::shared_ptr<functools::ComplexFunction> rhs
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
std::shared_ptr<functools::ComplexFunction> operator/(
	std::shared_ptr<functools::ConstantFunction> lhs,
	std::shared_ptr<functools::IdentityFunction> rhs
);
std::shared_ptr<functools::ComplexFunction> operator/(
	std::shared_ptr<functools::ConstantFunction> lhs,
	std::shared_ptr<functools::TrigonometryFunction> rhs
);
std::shared_ptr<functools::ComplexFunction> operator/(
	std::shared_ptr<functools::ConstantFunction> lhs,
	std::shared_ptr<functools::ComplexFunction> rhs
);

std::shared_ptr<functools::ComplexFunction> operator/(
	std::shared_ptr<functools::IdentityFunction> lhs,
	Type rhs
);
std::shared_ptr<functools::ComplexFunction> operator/(
	std::shared_ptr<functools::IdentityFunction> lhs,
	std::shared_ptr<functools::ConstantFunction> rhs
);
std::shared_ptr<functools::ConstantFunction> operator/(
	std::shared_ptr<functools::IdentityFunction> lhs,
	std::shared_ptr<functools::IdentityFunction> rhs
);
std::shared_ptr<functools::ComplexFunction> operator/(
	std::shared_ptr<functools::IdentityFunction> lhs,
	std::shared_ptr<functools::TrigonometryFunction> rhs
);
std::shared_ptr<functools::ComplexFunction> operator/(
	std::shared_ptr<functools::IdentityFunction> lhs,
	std::shared_ptr<functools::ComplexFunction> rhs
);

std::shared_ptr<functools::ComplexFunction> operator/(
	std::shared_ptr<functools::TrigonometryFunction> lhs,
	Type rhs
);
std::shared_ptr<functools::ComplexFunction> operator/(
	std::shared_ptr<functools::TrigonometryFunction> lhs,
	std::shared_ptr<functools::ConstantFunction> rhs
);
std::shared_ptr<functools::ComplexFunction> operator/(
	std::shared_ptr<functools::TrigonometryFunction> lhs,
	std::shared_ptr<functools::IdentityFunction> rhs
);
std::shared_ptr<functools::ComplexFunction> operator/(
	std::shared_ptr<functools::TrigonometryFunction> lhs,
	std::shared_ptr<functools::TrigonometryFunction> rhs
);
std::shared_ptr<functools::ComplexFunction> operator/(
	std::shared_ptr<functools::TrigonometryFunction> lhs,
	std::shared_ptr<functools::ComplexFunction> rhs
);

std::shared_ptr<functools::ComplexFunction> operator/(
	std::shared_ptr<functools::TrigonometryFunction> lhs,
	Type rhs
);
std::shared_ptr<functools::ComplexFunction> operator/(
	std::shared_ptr<functools::ComplexFunction> lhs,
	std::shared_ptr<functools::ConstantFunction> rhs
);
std::shared_ptr<functools::ComplexFunction> operator/(
	std::shared_ptr<functools::ComplexFunction> lhs,
	std::shared_ptr<functools::IdentityFunction> rhs
);
std::shared_ptr<functools::ComplexFunction> operator/(
	std::shared_ptr<functools::ComplexFunction> lhs,
	std::shared_ptr<functools::TrigonometryFunction> rhs
);
std::shared_ptr<functools::ComplexFunction> operator/(
	std::shared_ptr<functools::ComplexFunction> lhs,
	std::shared_ptr<functools::ComplexFunction> rhs
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
std::shared_ptr<functools::ComplexFunction> operator+(
	std::shared_ptr<functools::ConstantFunction> lhs,
	std::shared_ptr<functools::IdentityFunction> rhs
) {
	return std::make_shared<functools::ComplexFunction>(
		lhs,
		functools::FunctionOperator::PLUS,
		rhs
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
std::shared_ptr<functools::ComplexFunction> operator+(
	std::shared_ptr<functools::ConstantFunction> lhs,
	std::shared_ptr<functools::ComplexFunction> rhs
) {
	return std::make_shared<functools::ComplexFunction>(
		lhs,
		functools::FunctionOperator::PLUS,
		rhs
	);
}

std::shared_ptr<functools::ComplexFunction> operator+(
	std::shared_ptr<functools::IdentityFunction> lhs,
	Type rhs
) {
	return std::make_shared<functools::ConstantFunction>(rhs) + lhs;
}
std::shared_ptr<functools::ComplexFunction> operator+(
	std::shared_ptr<functools::IdentityFunction> lhs,
	std::shared_ptr<functools::ConstantFunction> rhs
) {
	return rhs + lhs;
}
std::shared_ptr<functools::ComplexFunction> operator+(
	std::shared_ptr<functools::IdentityFunction> lhs,
	std::shared_ptr<functools::IdentityFunction> rhs
) {
	return std::make_shared<functools::ComplexFunction>(
		std::make_shared<functools::ConstantFunction>(2),
		functools::FunctionOperator::TIMES,
		lhs
	);
}
std::shared_ptr<functools::ComplexFunction> operator+(
	std::shared_ptr<functools::IdentityFunction> lhs,
	std::shared_ptr<functools::TrigonometryFunction> rhs
) {
	return std::make_shared<functools::ComplexFunction>(
		lhs,
		functools::FunctionOperator::PLUS,
		rhs
	);
}
std::shared_ptr<functools::ComplexFunction> operator+(
	std::shared_ptr<functools::IdentityFunction> lhs,
	std::shared_ptr<functools::ComplexFunction> rhs
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
		std::make_shared<functools::ConstantFunction>(rhs)
	);
}
std::shared_ptr<functools::ComplexFunction> operator+(
	std::shared_ptr<functools::TrigonometryFunction> lhs,
	std::shared_ptr<functools::ConstantFunction> rhs
) {
	return std::make_shared<functools::ComplexFunction>(
		lhs,
		functools::FunctionOperator::PLUS,
		rhs
	);
}
std::shared_ptr<functools::ComplexFunction> operator+(
	std::shared_ptr<functools::TrigonometryFunction> lhs,
	std::shared_ptr<functools::IdentityFunction> rhs
) {
	return rhs + lhs;
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
std::shared_ptr<functools::ComplexFunction> operator+(
	std::shared_ptr<functools::TrigonometryFunction> lhs,
	std::shared_ptr<functools::ComplexFunction> rhs
) {
	return std::make_shared<functools::ComplexFunction>(
		lhs,
		functools::FunctionOperator::PLUS,
		rhs
	);
}

std::shared_ptr<functools::ComplexFunction> operator+(
	std::shared_ptr<functools::ComplexFunction> lhs,
	Type rhs
) {
	return std::make_shared<functools::ComplexFunction>(
		lhs,
		functools::FunctionOperator::PLUS,
		std::make_shared<functools::ConstantFunction>(rhs)
	);
}
std::shared_ptr<functools::ComplexFunction> operator+(
	std::shared_ptr<functools::ComplexFunction> lhs,
	std::shared_ptr<functools::ConstantFunction> rhs
) {
	return std::make_shared<functools::ComplexFunction>(
		lhs,
		functools::FunctionOperator::PLUS,
		rhs
	);
}
std::shared_ptr<functools::ComplexFunction> operator+(
	std::shared_ptr<functools::ComplexFunction> lhs,
	std::shared_ptr<functools::IdentityFunction> rhs
) {
	return rhs + lhs;
}
std::shared_ptr<functools::ComplexFunction> operator+(
	std::shared_ptr<functools::ComplexFunction> lhs,
	std::shared_ptr<functools::TrigonometryFunction> rhs
) {
	return std::make_shared<functools::ComplexFunction>(
		lhs,
		functools::FunctionOperator::PLUS,
		rhs
	);
}
std::shared_ptr<functools::ComplexFunction> operator+(
	std::shared_ptr<functools::ComplexFunction> lhs,
	std::shared_ptr<functools::ComplexFunction> rhs
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
std::shared_ptr<functools::ComplexFunction> operator-(
	std::shared_ptr<functools::ConstantFunction> lhs,
	std::shared_ptr<functools::IdentityFunction> rhs
) {
	return lhs + rhs * (-1);
}
std::shared_ptr<functools::ComplexFunction> operator-(
	std::shared_ptr<functools::ConstantFunction> lhs,
	std::shared_ptr<functools::TrigonometryFunction> rhs
) {
	return lhs + rhs * (-1);
}
std::shared_ptr<functools::ComplexFunction> operator-(
	std::shared_ptr<functools::ConstantFunction> lhs,
	std::shared_ptr<functools::ComplexFunction> rhs
) {
	return lhs + rhs * (-1);
}

std::shared_ptr<functools::ComplexFunction> operator-(
	std::shared_ptr<functools::IdentityFunction> lhs,
	Type rhs
) {
	return lhs + (-rhs);
}
std::shared_ptr<functools::ComplexFunction> operator-(
	std::shared_ptr<functools::IdentityFunction> lhs,
	std::shared_ptr<functools::ConstantFunction> rhs
) {
	return lhs - rhs->GetValue();
}
std::shared_ptr<functools::ConstantFunction> operator-(
	std::shared_ptr<functools::IdentityFunction> lhs,
	std::shared_ptr<functools::IdentityFunction> rhs
) {
	return std::make_shared<functools::ConstantFunction>(0);
}
std::shared_ptr<functools::ComplexFunction> operator-(
	std::shared_ptr<functools::IdentityFunction> lhs,
	std::shared_ptr<functools::TrigonometryFunction> rhs
) {
	return lhs + rhs * (-1);
}
std::shared_ptr<functools::ComplexFunction> operator-(
	std::shared_ptr<functools::IdentityFunction> lhs,
	std::shared_ptr<functools::ComplexFunction> rhs
) {
	return lhs + rhs * (-1);
}

std::shared_ptr<functools::ComplexFunction> operator-(
	std::shared_ptr<functools::TrigonometryFunction> lhs,
	Type rhs
) {
	return lhs + (-rhs);
}
std::shared_ptr<functools::ComplexFunction> operator-(
	std::shared_ptr<functools::TrigonometryFunction> lhs,
	std::shared_ptr<functools::ConstantFunction> rhs
) {
	return lhs + rhs * (-1);
}
std::shared_ptr<functools::ComplexFunction> operator-(
	std::shared_ptr<functools::TrigonometryFunction> lhs,
	std::shared_ptr<functools::IdentityFunction> rhs
) {
	return lhs + rhs * (-1);
}
std::shared_ptr<functools::ComplexFunction> operator-(
	std::shared_ptr<functools::TrigonometryFunction> lhs,
	std::shared_ptr<functools::TrigonometryFunction> rhs
) {
	return lhs + rhs * (-1);
}
std::shared_ptr<functools::ComplexFunction> operator-(
	std::shared_ptr<functools::TrigonometryFunction> lhs,
	std::shared_ptr<functools::ComplexFunction> rhs
) {
	return lhs + rhs * (-1);
}

std::shared_ptr<functools::ComplexFunction> operator-(
	std::shared_ptr<functools::ComplexFunction> lhs,
	Type rhs
) {
	return lhs + (-rhs);
}
std::shared_ptr<functools::ComplexFunction> operator-(
	std::shared_ptr<functools::ComplexFunction> lhs,
	std::shared_ptr<functools::ConstantFunction> rhs
) {
	return lhs + rhs * (-1);
}
std::shared_ptr<functools::ComplexFunction> operator-(
	std::shared_ptr<functools::ComplexFunction> lhs,
	std::shared_ptr<functools::IdentityFunction> rhs
) {
	return lhs + rhs * (-1);
}
std::shared_ptr<functools::ComplexFunction> operator-(
	std::shared_ptr<functools::ComplexFunction> lhs,
	std::shared_ptr<functools::TrigonometryFunction> rhs
) {
	return lhs + rhs * (-1);
}
std::shared_ptr<functools::ComplexFunction> operator-(
	std::shared_ptr<functools::ComplexFunction> lhs,
	std::shared_ptr<functools::ComplexFunction> rhs
) {
	return lhs + rhs * (-1);
}

// *

std::shared_ptr<functools::ConstantFunction> operator*(
	std::shared_ptr<functools::ConstantFunction> lhs,
	Type rhs
) {
	return std::make_shared<functools::ConstantFunction>(
		lhs->GetValue() * rhs
	);	
}
std::shared_ptr<functools::ConstantFunction> operator*(
	std::shared_ptr<functools::ConstantFunction> lhs,
	std::shared_ptr<functools::ConstantFunction> rhs
) {
	return std::make_shared<functools::ConstantFunction>(
		lhs->GetValue() * rhs->GetValue()
	);
}
std::shared_ptr<functools::ComplexFunction> operator*(
	std::shared_ptr<functools::ConstantFunction> lhs,
	std::shared_ptr<functools::IdentityFunction> rhs
) {
	return std::make_shared<functools::ComplexFunction>(
		lhs,
		functools::FunctionOperator::TIMES,
		rhs
	);
}
std::shared_ptr<functools::ComplexFunction> operator*(
	std::shared_ptr<functools::ConstantFunction> lhs,
	std::shared_ptr<functools::TrigonometryFunction> rhs
) {
	return std::make_shared<functools::ComplexFunction>(
		lhs,
		functools::FunctionOperator::TIMES,
		rhs
	);
}
std::shared_ptr<functools::ComplexFunction> operator*(
	std::shared_ptr<functools::ConstantFunction> lhs,
	std::shared_ptr<functools::ComplexFunction> rhs
) {
	return std::make_shared<functools::ComplexFunction>(
		lhs,
		functools::FunctionOperator::TIMES,
		rhs
	);
}

std::shared_ptr<functools::ComplexFunction> operator*(
	std::shared_ptr<functools::IdentityFunction> lhs,
	Type rhs
) {
	return std::make_shared<functools::ConstantFunction>(rhs) * lhs;
}
std::shared_ptr<functools::ComplexFunction> operator*(
	std::shared_ptr<functools::IdentityFunction> lhs,
	std::shared_ptr<functools::ConstantFunction> rhs
) {
	return rhs * lhs;
}
std::shared_ptr<functools::ComplexFunction> operator*(
	std::shared_ptr<functools::IdentityFunction> lhs,
	std::shared_ptr<functools::IdentityFunction> rhs
) {
	return std::make_shared<functools::ComplexFunction>(
		lhs,
		functools::FunctionOperator::POWER,
		std::make_shared<functools::ConstantFunction>(2)
	);
}
std::shared_ptr<functools::ComplexFunction> operator*(
	std::shared_ptr<functools::IdentityFunction> lhs,
	std::shared_ptr<functools::TrigonometryFunction> rhs
) {
	return std::make_shared<functools::ComplexFunction>(
		lhs,
		functools::FunctionOperator::TIMES,
		rhs
	);
}
std::shared_ptr<functools::ComplexFunction> operator*(
	std::shared_ptr<functools::IdentityFunction> lhs,
	std::shared_ptr<functools::ComplexFunction> rhs
) {
	return std::make_shared<functools::ComplexFunction>(
		lhs,
		functools::FunctionOperator::TIMES,
		rhs
	);
}

std::shared_ptr<functools::ComplexFunction> operator*(
	std::shared_ptr<functools::TrigonometryFunction> lhs,
	Type rhs
) {
	return std::make_shared<functools::ComplexFunction>(
		lhs,
		functools::FunctionOperator::TIMES,
		std::make_shared<functools::ConstantFunction>(rhs)
	);
}
std::shared_ptr<functools::ComplexFunction> operator*(
	std::shared_ptr<functools::TrigonometryFunction> lhs,
	std::shared_ptr<functools::ConstantFunction> rhs
) {
	return rhs * lhs;
}
std::shared_ptr<functools::ComplexFunction> operator*(
	std::shared_ptr<functools::TrigonometryFunction> lhs,
	std::shared_ptr<functools::IdentityFunction> rhs
) {
	return rhs * lhs;
}
std::shared_ptr<functools::ComplexFunction> operator*(
	std::shared_ptr<functools::TrigonometryFunction> lhs,
	std::shared_ptr<functools::TrigonometryFunction> rhs
) {
	return std::make_shared<functools::ComplexFunction>(
		lhs,
		functools::FunctionOperator::TIMES,
		rhs
	);
}
std::shared_ptr<functools::ComplexFunction> operator*(
	std::shared_ptr<functools::TrigonometryFunction> lhs,
	std::shared_ptr<functools::ComplexFunction> rhs
) {
	return std::make_shared<functools::ComplexFunction>(
		lhs,
		functools::FunctionOperator::TIMES,
		rhs
	);
}

std::shared_ptr<functools::ComplexFunction> operator*(
	std::shared_ptr<functools::ComplexFunction> lhs,
	Type rhs
) {
	return std::make_shared<functools::ComplexFunction>(
		lhs,
		functools::FunctionOperator::TIMES,
		std::make_shared<functools::ConstantFunction>(rhs)
	);
}
std::shared_ptr<functools::ComplexFunction> operator*(
	std::shared_ptr<functools::ComplexFunction> lhs,
	std::shared_ptr<functools::ConstantFunction> rhs
) {
	return rhs * lhs;
}
std::shared_ptr<functools::ComplexFunction> operator*(
	std::shared_ptr<functools::ComplexFunction> lhs,
	std::shared_ptr<functools::IdentityFunction> rhs
) {
	return rhs * lhs;
}
std::shared_ptr<functools::ComplexFunction> operator*(
	std::shared_ptr<functools::ComplexFunction> lhs,
	std::shared_ptr<functools::TrigonometryFunction> rhs
) {
	return rhs * lhs;
}
std::shared_ptr<functools::ComplexFunction> operator*(
	std::shared_ptr<functools::ComplexFunction> lhs,
	std::shared_ptr<functools::ComplexFunction> rhs
) {
	return std::make_shared<functools::ComplexFunction>(
		lhs,
		functools::FunctionOperator::TIMES,
		rhs
	);
}

// /

std::shared_ptr<functools::ConstantFunction> operator/(
	std::shared_ptr<functools::ConstantFunction> lhs,
	Type rhs
) {
	return std::make_shared<functools::ConstantFunction>(
		lhs->GetValue() / rhs
	);
}
std::shared_ptr<functools::ConstantFunction> operator/(
	std::shared_ptr<functools::ConstantFunction> lhs,
	std::shared_ptr<functools::ConstantFunction> rhs
) {
	return lhs / rhs->GetValue();
}
std::shared_ptr<functools::ComplexFunction> operator/(
	std::shared_ptr<functools::ConstantFunction> lhs,
	std::shared_ptr<functools::IdentityFunction> rhs
) {
	return std::make_shared<functools::ComplexFunction>(
		lhs,
		functools::FunctionOperator::DIVIDED,
		rhs
	);
}
std::shared_ptr<functools::ComplexFunction> operator/(
	std::shared_ptr<functools::ConstantFunction> lhs,
	std::shared_ptr<functools::TrigonometryFunction> rhs
) {
	return std::make_shared<functools::ComplexFunction>(
		lhs,
		functools::FunctionOperator::DIVIDED,
		rhs
	);
}
std::shared_ptr<functools::ComplexFunction> operator/(
	std::shared_ptr<functools::ConstantFunction> lhs,
	std::shared_ptr<functools::ComplexFunction> rhs
) {
	return std::make_shared<functools::ComplexFunction>(
		lhs,
		functools::FunctionOperator::DIVIDED,
		rhs
	);
}

std::shared_ptr<functools::ComplexFunction> operator/(
	std::shared_ptr<functools::IdentityFunction> lhs,
	Type rhs
) {
	return std::make_shared<functools::ComplexFunction>(
		lhs,
		functools::FunctionOperator::DIVIDED,
		std::make_shared<functools::ConstantFunction>(rhs)
	);
}
std::shared_ptr<functools::ComplexFunction> operator/(
	std::shared_ptr<functools::IdentityFunction> lhs,
	std::shared_ptr<functools::ConstantFunction> rhs
) {
	return lhs / rhs->GetValue();
}
std::shared_ptr<functools::ConstantFunction> operator/(
	std::shared_ptr<functools::IdentityFunction> lhs,
	std::shared_ptr<functools::IdentityFunction> rhs
) {
	return std::make_shared<functools::ConstantFunction>(1);
}
std::shared_ptr<functools::ComplexFunction> operator/(
	std::shared_ptr<functools::IdentityFunction> lhs,
	std::shared_ptr<functools::TrigonometryFunction> rhs
) {
	return std::make_shared<functools::ComplexFunction>(
		lhs,
		functools::FunctionOperator::DIVIDED,
		rhs
	);
}
std::shared_ptr<functools::ComplexFunction> operator/(
	std::shared_ptr<functools::IdentityFunction> lhs,
	std::shared_ptr<functools::ComplexFunction> rhs
) {
	return std::make_shared<functools::ComplexFunction>(
		lhs,
		functools::FunctionOperator::DIVIDED,
		rhs
	);
}

std::shared_ptr<functools::ComplexFunction> operator/(
	std::shared_ptr<functools::TrigonometryFunction> lhs,
	Type rhs
) {
	return std::make_shared<functools::ComplexFunction>(
		lhs,
		functools::FunctionOperator::DIVIDED,
		std::make_shared<functools::ConstantFunction>(rhs)
	);
}
std::shared_ptr<functools::ComplexFunction> operator/(
	std::shared_ptr<functools::TrigonometryFunction> lhs,
	std::shared_ptr<functools::ConstantFunction> rhs
) {
	return lhs / rhs->GetValue();
}
std::shared_ptr<functools::ComplexFunction> operator/(
	std::shared_ptr<functools::TrigonometryFunction> lhs,
	std::shared_ptr<functools::IdentityFunction> rhs
) {
	return std::make_shared<functools::ComplexFunction>(
		lhs,
		functools::FunctionOperator::DIVIDED,
		rhs
	);
}
std::shared_ptr<functools::ComplexFunction> operator/(
	std::shared_ptr<functools::TrigonometryFunction> lhs,
	std::shared_ptr<functools::TrigonometryFunction> rhs
) {
	return std::make_shared<functools::ComplexFunction>(
		lhs,
		functools::FunctionOperator::DIVIDED,
		rhs
	);
}
std::shared_ptr<functools::ComplexFunction> operator/(
	std::shared_ptr<functools::TrigonometryFunction> lhs,
	std::shared_ptr<functools::ComplexFunction> rhs
) {
	return std::make_shared<functools::ComplexFunction>(
		lhs,
		functools::FunctionOperator::DIVIDED,
		rhs
	);
}

std::shared_ptr<functools::ComplexFunction> operator/(
	std::shared_ptr<functools::ComplexFunction> lhs,
	Type rhs
) {
	return std::make_shared<functools::ComplexFunction>(
		lhs,
		functools::FunctionOperator::DIVIDED,
		std::make_shared<functools::ConstantFunction>(rhs)
	);
}
std::shared_ptr<functools::ComplexFunction> operator/(
	std::shared_ptr<functools::ComplexFunction> lhs,
	std::shared_ptr<functools::ConstantFunction> rhs
) {
	return lhs / rhs->GetValue();
}
std::shared_ptr<functools::ComplexFunction> operator/(
	std::shared_ptr<functools::ComplexFunction> lhs,
	std::shared_ptr<functools::IdentityFunction> rhs
) {
	return std::make_shared<functools::ComplexFunction>(
		lhs,
		functools::FunctionOperator::DIVIDED,
		rhs
	);
}
std::shared_ptr<functools::ComplexFunction> operator/(
	std::shared_ptr<functools::ComplexFunction> lhs,
	std::shared_ptr<functools::TrigonometryFunction> rhs
) {
	return std::make_shared<functools::ComplexFunction>(
		lhs,
		functools::FunctionOperator::DIVIDED,
		rhs
	);
}
std::shared_ptr<functools::ComplexFunction> operator/(
	std::shared_ptr<functools::ComplexFunction> lhs,
	std::shared_ptr<functools::ComplexFunction> rhs
) {
	return std::make_shared<functools::ComplexFunction>(
		lhs,
		functools::FunctionOperator::DIVIDED,
		rhs
	);
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
		if(auto rhsCast = std::dynamic_pointer_cast<functools::IdentityFunction>(rhs)) {
			return lhsCast + rhsCast;
		}
		if(auto rhsCast = std::dynamic_pointer_cast<functools::TrigonometryFunction>(rhs)) {
			return lhsCast + rhsCast;
		}
		if(auto rhsCast = std::dynamic_pointer_cast<functools::ComplexFunction>(rhs)) {
			return lhsCast + rhsCast;
		}
	}
	if(auto lhsCast = std::dynamic_pointer_cast<functools::IdentityFunction>(lhs)) {
		if(auto rhsCast = std::dynamic_pointer_cast<functools::ConstantFunction>(rhs)) {
			return lhsCast + rhsCast;
		}
		if(auto rhsCast = std::dynamic_pointer_cast<functools::IdentityFunction>(rhs)) {
			return lhsCast + rhsCast;
		}
		if(auto rhsCast = std::dynamic_pointer_cast<functools::TrigonometryFunction>(rhs)) {
			return lhsCast + rhsCast;
		}
		if(auto rhsCast = std::dynamic_pointer_cast<functools::ComplexFunction>(rhs)) {
			return lhsCast + rhsCast;
		}
	}
	if(auto lhsCast = std::dynamic_pointer_cast<functools::TrigonometryFunction>(lhs)) {
		if(auto rhsCast = std::dynamic_pointer_cast<functools::ConstantFunction>(rhs)) {
			return lhsCast + rhsCast;
		}
		if(auto rhsCast = std::dynamic_pointer_cast<functools::IdentityFunction>(rhs)) {
			return lhsCast + rhsCast;
		}
		if(auto rhsCast = std::dynamic_pointer_cast<functools::TrigonometryFunction>(rhs)) {
			return lhsCast + rhsCast;
		}
		if(auto rhsCast = std::dynamic_pointer_cast<functools::ComplexFunction>(rhs)) {
			return lhsCast + rhsCast;
		}
	}
	if(auto lhsCast = std::dynamic_pointer_cast<functools::ComplexFunction>(lhs)) {
		if(auto rhsCast = std::dynamic_pointer_cast<functools::ConstantFunction>(rhs)) {
			return lhsCast + rhsCast;
		}
		if(auto rhsCast = std::dynamic_pointer_cast<functools::IdentityFunction>(rhs)) {
			return lhsCast + rhsCast;
		}
		if(auto rhsCast = std::dynamic_pointer_cast<functools::TrigonometryFunction>(rhs)) {
			return lhsCast + rhsCast;
		}
		if(auto rhsCast = std::dynamic_pointer_cast<functools::ComplexFunction>(rhs)) {
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
	if(auto lhsCast = std::dynamic_pointer_cast<functools::IdentityFunction>(lhs)) {
		return lhsCast + rhs;
	}
	if(auto lhsCast = std::dynamic_pointer_cast<functools::TrigonometryFunction>(lhs)) {
		return lhsCast + rhs;
	}
	if(auto lhsCast = std::dynamic_pointer_cast<functools::ComplexFunction>(lhs)) {
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
		if(auto rhsCast = std::dynamic_pointer_cast<functools::IdentityFunction>(rhs)) {
			return lhsCast - rhsCast;
		}
		if(auto rhsCast = std::dynamic_pointer_cast<functools::TrigonometryFunction>(lhs)) {
			return lhsCast - rhsCast;
		}
		if(auto rhsCast = std::dynamic_pointer_cast<functools::ComplexFunction>(lhs)) {
			return lhsCast - rhsCast;
		}
	}
	if(auto lhsCast = std::dynamic_pointer_cast<functools::IdentityFunction>(lhs)) {
		if(auto rhsCast = std::dynamic_pointer_cast<functools::ConstantFunction>(rhs)) {
			return lhsCast - rhsCast;
		}
		if(auto rhsCast = std::dynamic_pointer_cast<functools::IdentityFunction>(rhs)) {
			return lhsCast - rhsCast;
		}
		if(auto rhsCast = std::dynamic_pointer_cast<functools::TrigonometryFunction>(lhs)) {
			return lhsCast - rhsCast;
		}
		if(auto rhsCast = std::dynamic_pointer_cast<functools::ComplexFunction>(lhs)) {
			return lhsCast - rhsCast;
		}
	}
	if(auto lhsCast = std::dynamic_pointer_cast<functools::TrigonometryFunction>(lhs)) {
		if(auto rhsCast = std::dynamic_pointer_cast<functools::ConstantFunction>(rhs)) {
			return lhsCast - rhsCast;
		}
		if(auto rhsCast = std::dynamic_pointer_cast<functools::IdentityFunction>(rhs)) {
			return lhsCast - rhsCast;
		}
		if(auto rhsCast = std::dynamic_pointer_cast<functools::TrigonometryFunction>(lhs)) {
			return lhsCast - rhsCast;
		}
		if(auto rhsCast = std::dynamic_pointer_cast<functools::ComplexFunction>(lhs)) {
			return lhsCast - rhsCast;
		}
	}
	if(auto lhsCast = std::dynamic_pointer_cast<functools::ComplexFunction>(lhs)) {
		if(auto rhsCast = std::dynamic_pointer_cast<functools::ConstantFunction>(rhs)) {
			return lhsCast - rhsCast;
		}
		if(auto rhsCast = std::dynamic_pointer_cast<functools::IdentityFunction>(rhs)) {
			return lhsCast - rhsCast;
		}
		if(auto rhsCast = std::dynamic_pointer_cast<functools::TrigonometryFunction>(lhs)) {
			return lhsCast - rhsCast;
		}
		if(auto rhsCast = std::dynamic_pointer_cast<functools::ComplexFunction>(lhs)) {
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
	if(auto lhsCast = std::dynamic_pointer_cast<functools::IdentityFunction>(lhs)) {
		return lhsCast - rhs;
	}
	if(auto lhsCast = std::dynamic_pointer_cast<functools::TrigonometryFunction>(lhs)) {
		return lhsCast - rhs;
	}
	if(auto lhsCast = std::dynamic_pointer_cast<functools::ComplexFunction>(lhs)) {
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
		if(auto rhsCast = std::dynamic_pointer_cast<functools::IdentityFunction>(rhs)) {
			return lhsCast * rhsCast;
		}
		if(auto rhsCast = std::dynamic_pointer_cast<functools::TrigonometryFunction>(rhs)) {
			return lhsCast * rhsCast;
		}
		if(auto rhsCast = std::dynamic_pointer_cast<functools::ComplexFunction>(rhs)) {
			return lhsCast * rhsCast;
		}
	}
	if(auto lhsCast = std::dynamic_pointer_cast<functools::IdentityFunction>(lhs)) {
		if(auto rhsCast = std::dynamic_pointer_cast<functools::ConstantFunction>(rhs)) {
			return lhsCast * rhsCast;
		}
		if(auto rhsCast = std::dynamic_pointer_cast<functools::IdentityFunction>(rhs)) {
			return lhsCast * rhsCast;
		}
		if(auto rhsCast = std::dynamic_pointer_cast<functools::TrigonometryFunction>(rhs)) {
			return lhsCast * rhsCast;
		}
		if(auto rhsCast = std::dynamic_pointer_cast<functools::ComplexFunction>(rhs)) {
			return lhsCast * rhsCast;
		}
	}
	if(auto lhsCast = std::dynamic_pointer_cast<functools::TrigonometryFunction>(lhs)) {
		if(auto rhsCast = std::dynamic_pointer_cast<functools::ConstantFunction>(rhs)) {
			return lhsCast * rhsCast;
		}
		if(auto rhsCast = std::dynamic_pointer_cast<functools::IdentityFunction>(rhs)) {
			return lhsCast * rhsCast;
		}
		if(auto rhsCast = std::dynamic_pointer_cast<functools::TrigonometryFunction>(rhs)) {
			return lhsCast * rhsCast;
		}
		if(auto rhsCast = std::dynamic_pointer_cast<functools::ComplexFunction>(rhs)) {
			return lhsCast * rhsCast;
		}
	}
	if(auto lhsCast = std::dynamic_pointer_cast<functools::ComplexFunction>(lhs)) {
		if(auto rhsCast = std::dynamic_pointer_cast<functools::ConstantFunction>(rhs)) {
			return lhsCast * rhsCast;
		}
		if(auto rhsCast = std::dynamic_pointer_cast<functools::IdentityFunction>(rhs)) {
			return lhsCast * rhsCast;
		}
		if(auto rhsCast = std::dynamic_pointer_cast<functools::TrigonometryFunction>(rhs)) {
			return lhsCast * rhsCast;
		}
		if(auto rhsCast = std::dynamic_pointer_cast<functools::ComplexFunction>(rhs)) {
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
	if(auto lhsCast = std::dynamic_pointer_cast<functools::IdentityFunction>(lhs)) {
		return lhsCast * rhs;
	}
	if(auto lhsCast = std::dynamic_pointer_cast<functools::TrigonometryFunction>(lhs)) {
		return lhsCast * rhs;
	}
	if(auto lhsCast = std::dynamic_pointer_cast<functools::ComplexFunction>(lhs)) {
		return lhsCast * rhs;
	}

	throw std::runtime_error("Invalid function");
}

std::shared_ptr<functools::Function> operator/(
	std::shared_ptr<functools::Function> lhs,
	std::shared_ptr<functools::Function> rhs
) {
	if(auto lhsCast = std::dynamic_pointer_cast<functools::ConstantFunction>(lhs)) {
		if(auto rhsCast = std::dynamic_pointer_cast<functools::ConstantFunction>(rhs)) {
			return lhsCast / rhsCast;
		}
		if(auto rhsCast = std::dynamic_pointer_cast<functools::IdentityFunction>(rhs)) {
			return lhsCast / rhsCast;
		}
		if(auto rhsCast = std::dynamic_pointer_cast<functools::TrigonometryFunction>(rhs)) {
			return lhsCast / rhsCast;
		}
		if(auto rhsCast = std::dynamic_pointer_cast<functools::ComplexFunction>(rhs)) {
			return lhsCast / rhsCast;
		}
	}
	if(auto lhsCast = std::dynamic_pointer_cast<functools::IdentityFunction>(lhs)) {
		if(auto rhsCast = std::dynamic_pointer_cast<functools::ConstantFunction>(rhs)) {
			return lhsCast / rhsCast;
		}
		if(auto rhsCast = std::dynamic_pointer_cast<functools::IdentityFunction>(rhs)) {
			return lhsCast / rhsCast;
		}
		if(auto rhsCast = std::dynamic_pointer_cast<functools::TrigonometryFunction>(rhs)) {
			return lhsCast / rhsCast;
		}
		if(auto rhsCast = std::dynamic_pointer_cast<functools::ComplexFunction>(rhs)) {
			return lhsCast / rhsCast;
		}
	}
	if(auto lhsCast = std::dynamic_pointer_cast<functools::TrigonometryFunction>(lhs)) {
		if(auto rhsCast = std::dynamic_pointer_cast<functools::ConstantFunction>(rhs)) {
			return lhsCast / rhsCast;
		}
		if(auto rhsCast = std::dynamic_pointer_cast<functools::IdentityFunction>(rhs)) {
			return lhsCast / rhsCast;
		}
		if(auto rhsCast = std::dynamic_pointer_cast<functools::TrigonometryFunction>(rhs)) {
			return lhsCast / rhsCast;
		}
		if(auto rhsCast = std::dynamic_pointer_cast<functools::ComplexFunction>(rhs)) {
			return lhsCast / rhsCast;
		}
	}
	if(auto lhsCast = std::dynamic_pointer_cast<functools::ComplexFunction>(lhs)) {
		if(auto rhsCast = std::dynamic_pointer_cast<functools::ConstantFunction>(rhs)) {
			return lhsCast / rhsCast;
		}
		if(auto rhsCast = std::dynamic_pointer_cast<functools::IdentityFunction>(rhs)) {
			return lhsCast / rhsCast;
		}
		if(auto rhsCast = std::dynamic_pointer_cast<functools::TrigonometryFunction>(rhs)) {
			return lhsCast / rhsCast;
		}
		if(auto rhsCast = std::dynamic_pointer_cast<functools::ComplexFunction>(rhs)) {
			return lhsCast / rhsCast;
		}
	}
	
	throw std::runtime_error("Invalid function");
}

std::shared_ptr<functools::Function> operator/(
	std::shared_ptr<functools::Function> lhs,
	Type rhs
) {
	if(auto lhsCast = std::dynamic_pointer_cast<functools::ConstantFunction>(lhs)) {
		return lhsCast / rhs;
	}
	if(auto lhsCast = std::dynamic_pointer_cast<functools::IdentityFunction>(lhs)) {
		return lhsCast / rhs;
	}
	if(auto lhsCast = std::dynamic_pointer_cast<functools::TrigonometryFunction>(lhs)) {
		return lhsCast / rhs;
	}
	if(auto lhsCast = std::dynamic_pointer_cast<functools::ComplexFunction>(lhs)) {
		return lhsCast / rhs;
	}

	throw std::runtime_error("Invalid function");
}
