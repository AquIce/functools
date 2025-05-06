#pragma once

#include <iostream>

#define LOG(value) std::cout << "[LOG] " << value << "\n";

#define NOIMP throw std::runtime_error("Not implemented yet") 
