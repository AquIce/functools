#pragma once

#include "debug.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <memory>
#include <stdexcept>
#include <string>

#define DEGREE_TYPE unsigned int
#define Type double

namespace functools {

	class Function {
	
	public:

		virtual Type Evaluate(Type x) { throw std::runtime_error("Invalid function"); }

		virtual std::unique_ptr<Function> GetDerivative() const = 0; 
		virtual std::unique_ptr<Function> GetPrimitive() const = 0; 

		virtual std::string Repr() const = 0;

	};

	// Polynomial Functions

	template <DEGREE_TYPE Degree = 0>
	class PolynomialFunction : public Function {
	
	public:

		PolynomialFunction(
			const std::array<Type, Degree + 1>& coefficients
		) :
			coefficients(coefficients)
		{}

		virtual Type Evaluate(Type x) override {
			Type result = 0;
			DEGREE_TYPE currentDegree = Degree;

			for(const auto& coefficient : coefficients) {
				result += coefficient * std::pow(x, currentDegree--);
			}
			return result;
		}

		std::unique_ptr<Function> GetDerivative() const override {

			std::array<Type, Degree> derivedCoefficients;

			for(DEGREE_TYPE i = 0; i < Degree; i++) {
				derivedCoefficients.at(i) = coefficients.at(i) * (Degree - i);
			}

			return std::make_unique<PolynomialFunction<Degree - 1>>(derivedCoefficients);
		}

		virtual std::unique_ptr<Function> GetPrimitive() const override {
			return nullptr;
		}

		virtual std::string Repr() const override {
			std::string repr = "";
			DEGREE_TYPE currentDegree = Degree;
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

	template <>
	class PolynomialFunction<0> : public Function {

	public:

		PolynomialFunction(
			const std::array<Type, 1>& coefficients
		) :
			coefficients(coefficients)
		{}

		Type Evaluate(Type x) override {
			return coefficients.at(0);
		}

		std::unique_ptr<Function> GetDerivative() const override {
			return std::make_unique<PolynomialFunction<0>>(
				std::array<Type, 1>({
					static_cast<Type>(0)
				})
			);
		}

		std::unique_ptr<Function> GetPrimitive() const override {
			return nullptr;
		}

		std::string Repr() const override {
			std::string repr = "";
			repr += std::to_string(coefficients.at(0));
			return repr;
		}

	private:
		std::array<Type, 1> coefficients;
	};


	// Trigonometry Functions

	enum class TrigonometryFunctionType {
		SIN,
		COS,
		TAN,
		COT
	};
}
