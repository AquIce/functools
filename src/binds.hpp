#pragma once

#include <cstdlib>
#include <vector>

template<typename T>
std::vector<T> operator+(
	const std::vector<T>& lhs,
	const std::vector<T>& rhs
) {
	std::vector<T> result;
	for(size_t i = 0; i < lhs.size(); i++) {
		result.push_back(lhs.at(i) + rhs.at(i));
	}
	return result;
}