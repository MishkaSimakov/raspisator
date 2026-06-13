#pragma once

#include <iostream>
#include <string>

#include "fwd.h"

template <typename Field>
class Variable {
  std::string name_;

  explicit Variable(std::string name) : name_(std::move(name)) {}

  friend class MILPProblem<Field>;
  friend class Expression<Field>;

 public:
  const std::string& get_name() const { return name_; }
};

template <typename Field>
bool operator==(const Variable<Field>& left, const Variable<Field>& right) {
  return left.get_name() == right.get_name();
}

template <typename Field>
struct std::hash<Variable<Field>> {
  size_t operator()(const Variable<Field>& variable) const {
    return std::hash<std::string>()(variable.get_name());
  }
};

template <typename Field>
std::ostream& operator<<(std::ostream& os, const Variable<Field>& variable) {
  os << variable.get_name();
  return os;
}
