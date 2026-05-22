#pragma once

#include <format>
#include <optional>
#include <stdexcept>
#include <string>

namespace mps {

class ParseError final : public std::exception {
  std::optional<size_t> line_;
  std::string message_;
  std::string what_;

 public:
  explicit ParseError(std::string message)
      : message_(std::move(message)), what_(message_) {}

  void set_line(size_t line) {
    line_ = line;
    what_ = std::format("MPS parse error on line {}: {}", line, message_);
  }

  const char* what() const noexcept override { return what_.c_str(); }

  std::optional<size_t> line() const noexcept { return line_; }
  const std::string& message() const noexcept { return message_; }
};

}  // namespace mps
