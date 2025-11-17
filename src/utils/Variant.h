#pragma once

template <typename... Args>
struct Overload : Args... {
  explicit Overload(Args... args) : Args(args)... {}

  using Args::operator()...;
};

template <typename U, typename... Args>
U variant_cast(std::variant<Args...> variant) {
  return std::visit([](auto& value) { return U(std::move(value)); }, variant);
}
