#pragma once

#include <array>
#include <cstdlib>
#include <memory>

template<typename T, size_t S>
std::array<std::shared_ptr<T>, S> operator+(
	const std::array<std::shared_ptr<T>, S>& lhs,
	const std::array<std::shared_ptr<T>, S>& rhs
) {
	std::array<std::shared_ptr<T>, S> result;
	for(size_t i = 0; i < S; i++) {
		result.at(i) = *(lhs.at(i).get()) + rhs.at(i);
	}
	return result;
}
