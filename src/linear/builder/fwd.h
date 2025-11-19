#pragma once

template <typename Field>
class ProblemBuilder;

template <typename Field>
class Variable;

template <typename Field>
class Expression;

namespace details {
template <typename T>
struct expression_like : std::false_type {};

template <typename Field>
struct expression_like<Expression<Field>> : std::true_type {};

template <typename Field>
struct expression_like<Variable<Field>> : std::true_type {};

template <typename T>
concept ExpressionLike = expression_like<std::decay_t<T>>::value;

}  // namespace details
