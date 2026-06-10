#pragma once

#include <type_traits>

template <typename T, bool is_const>
using maybe_const_t = std::conditional_t<is_const, const T, T>;
