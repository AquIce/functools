#pragma once

#include "debug.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>

#define DegreeType char
#define Type double

namespace functools {

	class Function {
	
	public:

		virtual Type Evaluate(Type x) { throw std::runtime_error("Invalid function"); }

		virtual std::unique_ptr<Function> GetDerivative() const = 0; 
		virtual std::unique_ptr<Function> GetPrimitive() const = 0; 

		virtual std::string Repr() const = 0;

	};

	// --- 
	// Polynomial Functions
	// ---

	template <DegreeType Degree = 0>
	class PolynomialFunction : public Function {
	
	public:

		PolynomialFunction(
			const std::array<Type, Degree + 1>& coefficients
		) :
			coefficients(coefficients)
		{}

		Type Evaluate(Type x) override {
			Type result = 0;
			DegreeType currentDegree = Degree;

			for(const auto& coefficient : coefficients) {
				result += coefficient * std::pow(x, currentDegree--);
			}
			return result;
		}

		std::unique_ptr<Function> GetDerivative() const override {

			std::array<Type, Degree> derivedCoefficients;

			for(DegreeType i = 0; i < Degree; i++) {
				derivedCoefficients.at(i) = coefficients.at(i) * (Degree - i);
			}

			return std::make_unique<PolynomialFunction<Degree - 1>>(derivedCoefficients);
		}

		std::unique_ptr<Function> GetPrimitive() const override {

			std::array<Type, Degree + 2> primitiveCoefficients;

			for(DegreeType i = 0; i <= Degree; i++) {
				primitiveCoefficients.at(i) = coefficients.at(i) / (Degree - i + 1);
			}

			return std::make_unique<PolynomialFunction<Degree + 1>>(primitiveCoefficients);
		}

		std::string Repr() const override {
			std::string repr = "";
			DegreeType currentDegree = Degree;
			for(const auto& coefficient : coefficients) {
				
				if(coefficient == 0) {
					currentDegree--;
					continue;
				}

				std::string coeffString = std::to_string(coefficient);

				if(coefficient == 1) {
					coeffString = currentDegree == 0 ? "1" : "";
				} else if(coefficient == -1) {
					coeffString = currentDegree == 0 ? "-1" : "-";
				}

				std::string xForm = "x^" + std::to_string(currentDegree);
				
				if(currentDegree == 0) {
					xForm = "";
				} else if(currentDegree == 1) {
					xForm = "x";
				}

				repr +=
					((coefficient > 0 && currentDegree != Degree) ? "+ " : "")
					+ coeffString + xForm + " "; 

				currentDegree--;
			}

			return repr;
		}

	private:
		std::array<Type, Degree + 1> coefficients;
	};

	// For Degree = 0

	template <>
	Type PolynomialFunction<0>::Evaluate(Type x) {
		return coefficients.at(0);
	}

	template <>
	std::unique_ptr<Function> PolynomialFunction<0>::GetDerivative() const {
		return std::make_unique<PolynomialFunction<0>>(
			std::array<Type, 1>({
				static_cast<Type>(0)
			})
		);
	}

	template <>
	std::string PolynomialFunction<0>::Repr() const {
		return std::to_string(coefficients.at(0));
	}

	// For Degree = std::numeric_limits<DegreeType>::max()

	template <>
	std::unique_ptr<Function> PolynomialFunction<std::numeric_limits<DegreeType>::max()>::GetPrimitive() const {
		throw std::runtime_error("Exceeded maximum allowable degree for polynomial functions.");
	}

	// --- 
	// Trigonometry Functions
	// ---

	enum class TrigonometryFunctionType {
		SIN,
		COS,
		TAN,
		COT
	};
	
	class TrigonometryFunction : public Function {
	
	public:

		TrigonometryFunction(
			const TrigonometryFunctionType type,
			std::unique_ptr<Function> inner
		) :
			type(type),
			inner(inner)
		{}

		Type Evaluate(Type x) override {

		}

		std::unique_ptr<Function> GetDerivative() const override {

		}

		std::unique_ptr<Function> GetPrimitive() const override {

		}

		std::string Repr() const override {

		}

	private:
		TrigonometryFunctionType type;
		std::unique_ptr<Function> inner;
	};

}
