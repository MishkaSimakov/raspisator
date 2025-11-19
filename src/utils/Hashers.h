#pragma once

#include <utility>

inline auto hash_fn = []<typename T>(const T& value) {
  return std::hash<T>()(value);
};

struct StreamHasher {
 private:
  size_t current_ = 0;

 public:
  template <typename T>
  StreamHasher& operator<<(const T& value) {
    current_ ^=
        std::hash<T>()(value) + 0xeeffccdd + (current_ << 5) + (current_ >> 3);

    return *this;
  }

  size_t get_hash() const { return current_; }
};

template <typename... Args>
struct TupleHasher {
  size_t operator()(const Args&... args) const {
    StreamHasher hasher{};
    (hasher << ... << args);
    return hasher.get_hash();
  }

  size_t operator()(const std::tuple<Args...>& tuple) const {
    return [this, &tuple]<size_t... I>(std::index_sequence<I...>) {
      return operator()(std::get<I>(tuple)...);
    }(std::index_sequence_for<Args...>{});
  }
};

inline auto tuple_hasher_fn = []<typename... Args>(const Args&... args) {
  return TupleHasher<Args...>()(args...);
};

template <typename U, typename V>
struct std::hash<std::pair<U, V>> {
  size_t operator()(const std::pair<U, V>& pair) const noexcept {
    return tuple_hasher_fn(pair.first, pair.second);
  }
};

template <typename... Args>
struct std::hash<std::tuple<Args...>> {
  size_t operator()(const std::tuple<Args...>& tuple) const noexcept {
    return [&tuple]<size_t... Is>(std::index_sequence<Is...>) {
      return tuple_hasher_fn(std::get<Is>(tuple)...);
    }(std::make_index_sequence<sizeof...(Args)>());
  }
};
