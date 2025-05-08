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

// ---
// Macros
// ---

#define HANDLE_POLYNOMIAL_CASE(T, code) \
    if(auto child = std::dynamic_pointer_cast<functools::PolynomialFunction<T>>(lhs)) code

#define HANDLE_POLYNOMIAL(lhs, code) \
    HANDLE_POLYNOMIAL_CASE(0, code) \
    HANDLE_POLYNOMIAL_CASE(1, code) HANDLE_POLYNOMIAL_CASE(2, code) HANDLE_POLYNOMIAL_CASE(3, code) HANDLE_POLYNOMIAL_CASE(4, code) \
    HANDLE_POLYNOMIAL_CASE(5, code) HANDLE_POLYNOMIAL_CASE(6, code) HANDLE_POLYNOMIAL_CASE(7, code) HANDLE_POLYNOMIAL_CASE(8, code) \
    HANDLE_POLYNOMIAL_CASE(9, code) HANDLE_POLYNOMIAL_CASE(10, code) HANDLE_POLYNOMIAL_CASE(11, code) HANDLE_POLYNOMIAL_CASE(12, code) \
    HANDLE_POLYNOMIAL_CASE(13, code) HANDLE_POLYNOMIAL_CASE(14, code) HANDLE_POLYNOMIAL_CASE(15, code) HANDLE_POLYNOMIAL_CASE(16, code) \
    HANDLE_POLYNOMIAL_CASE(17, code) HANDLE_POLYNOMIAL_CASE(18, code) HANDLE_POLYNOMIAL_CASE(19, code) HANDLE_POLYNOMIAL_CASE(20, code) \
    HANDLE_POLYNOMIAL_CASE(21, code) HANDLE_POLYNOMIAL_CASE(22, code) HANDLE_POLYNOMIAL_CASE(23, code) HANDLE_POLYNOMIAL_CASE(24, code) \
    HANDLE_POLYNOMIAL_CASE(25, code) HANDLE_POLYNOMIAL_CASE(26, code) HANDLE_POLYNOMIAL_CASE(27, code) HANDLE_POLYNOMIAL_CASE(28, code) \
    HANDLE_POLYNOMIAL_CASE(29, code) HANDLE_POLYNOMIAL_CASE(30, code) HANDLE_POLYNOMIAL_CASE(31, code) HANDLE_POLYNOMIAL_CASE(32, code) \
    HANDLE_POLYNOMIAL_CASE(33, code) HANDLE_POLYNOMIAL_CASE(34, code) HANDLE_POLYNOMIAL_CASE(35, code) HANDLE_POLYNOMIAL_CASE(36, code) \
    HANDLE_POLYNOMIAL_CASE(37, code) HANDLE_POLYNOMIAL_CASE(38, code) HANDLE_POLYNOMIAL_CASE(39, code) HANDLE_POLYNOMIAL_CASE(40, code) \
    HANDLE_POLYNOMIAL_CASE(41, code) HANDLE_POLYNOMIAL_CASE(42, code) HANDLE_POLYNOMIAL_CASE(43, code) HANDLE_POLYNOMIAL_CASE(44, code) \
    HANDLE_POLYNOMIAL_CASE(45, code) HANDLE_POLYNOMIAL_CASE(46, code) HANDLE_POLYNOMIAL_CASE(47, code) HANDLE_POLYNOMIAL_CASE(48, code) \
    HANDLE_POLYNOMIAL_CASE(49, code) HANDLE_POLYNOMIAL_CASE(50, code) HANDLE_POLYNOMIAL_CASE(51, code) HANDLE_POLYNOMIAL_CASE(52, code) \
    HANDLE_POLYNOMIAL_CASE(53, code) HANDLE_POLYNOMIAL_CASE(54, code) HANDLE_POLYNOMIAL_CASE(55, code) HANDLE_POLYNOMIAL_CASE(56, code) \
    HANDLE_POLYNOMIAL_CASE(57, code) HANDLE_POLYNOMIAL_CASE(58, code) HANDLE_POLYNOMIAL_CASE(59, code) HANDLE_POLYNOMIAL_CASE(60, code) \
    HANDLE_POLYNOMIAL_CASE(61, code) HANDLE_POLYNOMIAL_CASE(62, code) HANDLE_POLYNOMIAL_CASE(63, code) HANDLE_POLYNOMIAL_CASE(64, code) \
    HANDLE_POLYNOMIAL_CASE(65, code) HANDLE_POLYNOMIAL_CASE(66, code) HANDLE_POLYNOMIAL_CASE(67, code) HANDLE_POLYNOMIAL_CASE(68, code) \
    HANDLE_POLYNOMIAL_CASE(69, code) HANDLE_POLYNOMIAL_CASE(70, code) HANDLE_POLYNOMIAL_CASE(71, code) HANDLE_POLYNOMIAL_CASE(72, code) \
    HANDLE_POLYNOMIAL_CASE(73, code) HANDLE_POLYNOMIAL_CASE(74, code) HANDLE_POLYNOMIAL_CASE(75, code) HANDLE_POLYNOMIAL_CASE(76, code) \
    HANDLE_POLYNOMIAL_CASE(77, code) HANDLE_POLYNOMIAL_CASE(78, code) HANDLE_POLYNOMIAL_CASE(79, code) HANDLE_POLYNOMIAL_CASE(80, code) \
    HANDLE_POLYNOMIAL_CASE(81, code) HANDLE_POLYNOMIAL_CASE(82, code) HANDLE_POLYNOMIAL_CASE(83, code) HANDLE_POLYNOMIAL_CASE(84, code) \
    HANDLE_POLYNOMIAL_CASE(85, code) HANDLE_POLYNOMIAL_CASE(86, code) HANDLE_POLYNOMIAL_CASE(87, code) HANDLE_POLYNOMIAL_CASE(88, code) \
    HANDLE_POLYNOMIAL_CASE(89, code) HANDLE_POLYNOMIAL_CASE(90, code) HANDLE_POLYNOMIAL_CASE(91, code) HANDLE_POLYNOMIAL_CASE(92, code) \
    HANDLE_POLYNOMIAL_CASE(93, code) HANDLE_POLYNOMIAL_CASE(94, code) HANDLE_POLYNOMIAL_CASE(95, code) HANDLE_POLYNOMIAL_CASE(96, code) \
    HANDLE_POLYNOMIAL_CASE(97, code) HANDLE_POLYNOMIAL_CASE(98, code) HANDLE_POLYNOMIAL_CASE(99, code) HANDLE_POLYNOMIAL_CASE(100, code) \
    HANDLE_POLYNOMIAL_CASE(101, code) HANDLE_POLYNOMIAL_CASE(102, code) HANDLE_POLYNOMIAL_CASE(103, code) HANDLE_POLYNOMIAL_CASE(104, code) \
    HANDLE_POLYNOMIAL_CASE(105, code) HANDLE_POLYNOMIAL_CASE(106, code) HANDLE_POLYNOMIAL_CASE(107, code) HANDLE_POLYNOMIAL_CASE(108, code) \
    HANDLE_POLYNOMIAL_CASE(109, code) HANDLE_POLYNOMIAL_CASE(110, code) HANDLE_POLYNOMIAL_CASE(111, code) HANDLE_POLYNOMIAL_CASE(112, code) \
    HANDLE_POLYNOMIAL_CASE(113, code) HANDLE_POLYNOMIAL_CASE(114, code) HANDLE_POLYNOMIAL_CASE(115, code) HANDLE_POLYNOMIAL_CASE(116, code) \
    HANDLE_POLYNOMIAL_CASE(117, code) HANDLE_POLYNOMIAL_CASE(118, code) HANDLE_POLYNOMIAL_CASE(119, code) HANDLE_POLYNOMIAL_CASE(120, code) \
    HANDLE_POLYNOMIAL_CASE(121, code) HANDLE_POLYNOMIAL_CASE(122, code) HANDLE_POLYNOMIAL_CASE(123, code) HANDLE_POLYNOMIAL_CASE(124, code) \
    HANDLE_POLYNOMIAL_CASE(125, code) HANDLE_POLYNOMIAL_CASE(126, code) HANDLE_POLYNOMIAL_CASE(127, code) HANDLE_POLYNOMIAL_CASE(128, code) \
    HANDLE_POLYNOMIAL_CASE(129, code) HANDLE_POLYNOMIAL_CASE(130, code) HANDLE_POLYNOMIAL_CASE(131, code) HANDLE_POLYNOMIAL_CASE(132, code) \
    HANDLE_POLYNOMIAL_CASE(133, code) HANDLE_POLYNOMIAL_CASE(134, code) HANDLE_POLYNOMIAL_CASE(135, code) HANDLE_POLYNOMIAL_CASE(136, code) \
    HANDLE_POLYNOMIAL_CASE(137, code) HANDLE_POLYNOMIAL_CASE(138, code) HANDLE_POLYNOMIAL_CASE(139, code) HANDLE_POLYNOMIAL_CASE(140, code) \
    HANDLE_POLYNOMIAL_CASE(141, code) HANDLE_POLYNOMIAL_CASE(142, code) HANDLE_POLYNOMIAL_CASE(143, code) HANDLE_POLYNOMIAL_CASE(144, code) \
    HANDLE_POLYNOMIAL_CASE(145, code) HANDLE_POLYNOMIAL_CASE(146, code) HANDLE_POLYNOMIAL_CASE(147, code) HANDLE_POLYNOMIAL_CASE(148, code) \
    HANDLE_POLYNOMIAL_CASE(149, code) HANDLE_POLYNOMIAL_CASE(150, code) HANDLE_POLYNOMIAL_CASE(151, code) HANDLE_POLYNOMIAL_CASE(152, code) \
    HANDLE_POLYNOMIAL_CASE(153, code) HANDLE_POLYNOMIAL_CASE(154, code) HANDLE_POLYNOMIAL_CASE(155, code) HANDLE_POLYNOMIAL_CASE(156, code) \
    HANDLE_POLYNOMIAL_CASE(157, code) HANDLE_POLYNOMIAL_CASE(158, code) HANDLE_POLYNOMIAL_CASE(159, code) HANDLE_POLYNOMIAL_CASE(160, code) \
    HANDLE_POLYNOMIAL_CASE(161, code) HANDLE_POLYNOMIAL_CASE(162, code) HANDLE_POLYNOMIAL_CASE(163, code) HANDLE_POLYNOMIAL_CASE(164, code) \
    HANDLE_POLYNOMIAL_CASE(165, code) HANDLE_POLYNOMIAL_CASE(166, code) HANDLE_POLYNOMIAL_CASE(167, code) HANDLE_POLYNOMIAL_CASE(168, code) \
    HANDLE_POLYNOMIAL_CASE(169, code) HANDLE_POLYNOMIAL_CASE(170, code) HANDLE_POLYNOMIAL_CASE(171, code) HANDLE_POLYNOMIAL_CASE(172, code) \
    HANDLE_POLYNOMIAL_CASE(173, code) HANDLE_POLYNOMIAL_CASE(174, code) HANDLE_POLYNOMIAL_CASE(175, code) HANDLE_POLYNOMIAL_CASE(176, code) \
    HANDLE_POLYNOMIAL_CASE(177, code) HANDLE_POLYNOMIAL_CASE(178, code) HANDLE_POLYNOMIAL_CASE(179, code) HANDLE_POLYNOMIAL_CASE(180, code) \
    HANDLE_POLYNOMIAL_CASE(181, code) HANDLE_POLYNOMIAL_CASE(182, code) HANDLE_POLYNOMIAL_CASE(183, code) HANDLE_POLYNOMIAL_CASE(184, code) \
    HANDLE_POLYNOMIAL_CASE(185, code) HANDLE_POLYNOMIAL_CASE(186, code) HANDLE_POLYNOMIAL_CASE(187, code) HANDLE_POLYNOMIAL_CASE(188, code) \
    HANDLE_POLYNOMIAL_CASE(189, code) HANDLE_POLYNOMIAL_CASE(190, code) HANDLE_POLYNOMIAL_CASE(191, code) HANDLE_POLYNOMIAL_CASE(192, code) \
    HANDLE_POLYNOMIAL_CASE(193, code) HANDLE_POLYNOMIAL_CASE(194, code) HANDLE_POLYNOMIAL_CASE(195, code) HANDLE_POLYNOMIAL_CASE(196, code) \
    HANDLE_POLYNOMIAL_CASE(197, code) HANDLE_POLYNOMIAL_CASE(198, code) HANDLE_POLYNOMIAL_CASE(199, code) HANDLE_POLYNOMIAL_CASE(200, code) \
    HANDLE_POLYNOMIAL_CASE(201, code) HANDLE_POLYNOMIAL_CASE(202, code) HANDLE_POLYNOMIAL_CASE(203, code) HANDLE_POLYNOMIAL_CASE(204, code) \
    HANDLE_POLYNOMIAL_CASE(205, code) HANDLE_POLYNOMIAL_CASE(206, code) HANDLE_POLYNOMIAL_CASE(207, code) HANDLE_POLYNOMIAL_CASE(208, code) \
    HANDLE_POLYNOMIAL_CASE(209, code) HANDLE_POLYNOMIAL_CASE(210, code) HANDLE_POLYNOMIAL_CASE(211, code) HANDLE_POLYNOMIAL_CASE(212, code) \
    HANDLE_POLYNOMIAL_CASE(213, code) HANDLE_POLYNOMIAL_CASE(214, code) HANDLE_POLYNOMIAL_CASE(215, code) HANDLE_POLYNOMIAL_CASE(216, code) \
    HANDLE_POLYNOMIAL_CASE(217, code) HANDLE_POLYNOMIAL_CASE(218, code) HANDLE_POLYNOMIAL_CASE(219, code) HANDLE_POLYNOMIAL_CASE(220, code) \
    HANDLE_POLYNOMIAL_CASE(221, code) HANDLE_POLYNOMIAL_CASE(222, code) HANDLE_POLYNOMIAL_CASE(223, code) HANDLE_POLYNOMIAL_CASE(224, code) \
    HANDLE_POLYNOMIAL_CASE(225, code) HANDLE_POLYNOMIAL_CASE(226, code) HANDLE_POLYNOMIAL_CASE(227, code) HANDLE_POLYNOMIAL_CASE(228, code) \
    HANDLE_POLYNOMIAL_CASE(229, code) HANDLE_POLYNOMIAL_CASE(230, code) HANDLE_POLYNOMIAL_CASE(231, code) HANDLE_POLYNOMIAL_CASE(232, code) \
    HANDLE_POLYNOMIAL_CASE(233, code) HANDLE_POLYNOMIAL_CASE(234, code) HANDLE_POLYNOMIAL_CASE(235, code) HANDLE_POLYNOMIAL_CASE(236, code) \
    HANDLE_POLYNOMIAL_CASE(237, code) HANDLE_POLYNOMIAL_CASE(238, code) HANDLE_POLYNOMIAL_CASE(239, code) HANDLE_POLYNOMIAL_CASE(240, code) \
    HANDLE_POLYNOMIAL_CASE(241, code) HANDLE_POLYNOMIAL_CASE(242, code) HANDLE_POLYNOMIAL_CASE(243, code) HANDLE_POLYNOMIAL_CASE(244, code) \
    HANDLE_POLYNOMIAL_CASE(245, code) HANDLE_POLYNOMIAL_CASE(246, code) HANDLE_POLYNOMIAL_CASE(247, code) HANDLE_POLYNOMIAL_CASE(248, code) \
    HANDLE_POLYNOMIAL_CASE(249, code) HANDLE_POLYNOMIAL_CASE(250, code) HANDLE_POLYNOMIAL_CASE(251, code) HANDLE_POLYNOMIAL_CASE(252, code) \
    HANDLE_POLYNOMIAL_CASE(253, code) HANDLE_POLYNOMIAL_CASE(254, code) HANDLE_POLYNOMIAL_CASE(255, code)


template <DegreeType lhs, DegreeType rhs>
struct MaxDegree {
public:
	static const DegreeType value = (lhs > rhs) ? lhs : rhs;
};
template <DegreeType lhs, DegreeType rhs>
struct MinDegree {
public:
	static const DegreeType value = (lhs > rhs) ? rhs : lhs;
};


template<typename T, size_t S, size_t S_>
constexpr std::array<T, MaxDegree<S, S_>::value> GetMaxArray(
	std::array<T, S> lhs,
	std::array<T, S_> rhs
) {
	if constexpr (S > S_) {
		return lhs;
	}
	return rhs;
}
template<typename T, size_t S, size_t S_>
constexpr std::array<T, MinDegree<S, S_>::value> GetMinArray(
	std::array<T, S> lhs,
	std::array<T, S_> rhs
) {
	if constexpr (S > S_) {
		return rhs;
	}
	return lhs;
}

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

std::shared_ptr<functools::Function> operator+(
	std::shared_ptr<functools::Function> lhs,
	std::shared_ptr<functools::Function> rhs
);
std::shared_ptr<functools::Function> operator+(
	std::shared_ptr<functools::Function> lhs,
	Type rhs
);

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
	std::array<std::shared_ptr<functools::Function>, Degree + 1> res = rhs->GetCoefficients();
	res.at(Degree) = res.at(Degree) + lhs->GetValue();
	return std::make_shared<functools::PolynomialFunction<Degree>>(
		res
	);
}
template <DegreeType Degree>
std::shared_ptr<functools::PolynomialFunction<Degree>> operator+(
	std::shared_ptr<functools::PolynomialFunction<Degree>> lhs,
	Type rhs
) {
	std::array<std::shared_ptr<functools::Function>, Degree + 1> res = lhs->GetCoefficients();
	res.at(Degree) = res.at(Degree) + rhs;
	return std::make_shared<functools::PolynomialFunction<Degree>>(
		res
	);
}
template <DegreeType Degree>
std::shared_ptr<functools::PolynomialFunction<Degree>> operator+(
	std::shared_ptr<functools::PolynomialFunction<Degree>> lhs,
	std::shared_ptr<functools::ConstantFunction> rhs
) {
	return rhs + lhs;
}
template <DegreeType lhsDegree, DegreeType rhsDegree>
std::shared_ptr<functools::PolynomialFunction<MaxDegree<lhsDegree, rhsDegree>::value>> operator+(
	std::shared_ptr<functools::PolynomialFunction<lhsDegree>> lhs,
	std::shared_ptr<functools::PolynomialFunction<rhsDegree>> rhs
) {

	if constexpr (lhsDegree == rhsDegree) {

		std::array<std::shared_ptr<functools::Function>, lhsDegree + 1> res = lhs->GetCoefficients() + rhs->GetCoefficients();

		return std::make_shared<functools::PolynomialFunction<lhsDegree>>(res);
	}

	constexpr DegreeType maxDegree = MaxDegree<lhsDegree, rhsDegree>::value;
	constexpr DegreeType minDegree = MinDegree<lhsDegree, rhsDegree>::value;

	std::array<std::shared_ptr<functools::Function>, maxDegree + 1> maxCoeffs = GetMaxArray<
		std::shared_ptr<functools::Function>, lhsDegree + 1, rhsDegree + 1
	>(lhs->GetCoefficients(), rhs->GetCoefficients());
	std::array<std::shared_ptr<functools::Function>, minDegree + 1> minCoeffs = GetMinArray<
		std::shared_ptr<functools::Function>, lhsDegree + 1, rhsDegree + 1
	>(lhs->GetCoefficients(), rhs->GetCoefficients());

	std::array<std::shared_ptr<functools::Function>, maxDegree + 1> res;
	std::copy(maxCoeffs.begin(), maxCoeffs.end(), maxCoeffs.begin());

	for(DegreeType i = 0; i <= minDegree; i++) {
		res.at(maxDegree - minDegree + i) = res.at(maxDegree - minDegree + i) + minCoeffs.at(i);
	}

	return std::make_shared<functools::PolynomialFunction<maxDegree>>(
		res
	);
}

std::shared_ptr<functools::Function> operator+(
	std::shared_ptr<functools::Function> lhs,
	std::shared_ptr<functools::Function> rhs
) {
	if(auto lhsCast = std::dynamic_pointer_cast<functools::ConstantFunction>(lhs)) {
		if(auto rhsCast = std::dynamic_pointer_cast<functools::ConstantFunction>(rhs)) {
			return lhsCast + rhsCast;
		}
		HANDLE_POLYNOMIAL(rhs, {
			return lhsCast + child;
		});
	}
	// if(auto lhsCast = std::dynamic_pointer_cast<functools::PolynomialFunction<Degree>>(lhs)) {
	// 	if(auto rhsCast = std::dynamic_pointer_cast<functools::ConstantFunction>(rhs)) {
	// 		return lhsCast + rhsCast;
	// 	}
	// 	if(auto rhsCast = std::dynamic_pointer_cast<functools::PolynomialFunction<Degree>>(lhs)) {
	// 		return lhsCast + rhsCast;
	// 	}
	// }
	NOIMP;
}

std::shared_ptr<functools::Function> operator+(
	std::shared_ptr<functools::Function> lhs,
	Type rhs
) {
	if(auto lhsCast = std::dynamic_pointer_cast<functools::ConstantFunction>(lhs)) {
		return lhsCast + rhs;
	}
	HANDLE_POLYNOMIAL(lhs, {
		return child + rhs;
	})

	throw std::runtime_error("Invalid function");
}