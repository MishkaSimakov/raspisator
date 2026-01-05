#pragma once

#include <chrono>
#include <fstream>

#include "linear/matrix/Matrix.h"
#include "linear/model/LP.h"
#include "linear/sparse/CSCMatrix.h"

namespace simplex {


// double to_double(std::initializer_list<uint8_t> bytes) {
//   double result;
//   memcpy(&result, bytes.begin(), sizeof(double));
//
//   return result;
// }
//
// void print_double(double value) {
//   auto bytes = reinterpret_cast<uint8_t*>(&value);
//
//   std::cout << "{";
//
//   for (size_t i = 0; i < sizeof(value); ++i) {
//     std::cout << static_cast<int>(bytes[i]) << ", ";
//   }
//
//   std::cout << "}";
// }


template <typename Field>
class SimplexCoreDump {
  const CSCMatrix<Field>& A_;
  const Matrix<Field>& b_;
  const Matrix<Field>& c_;
  const Matrix<Field>& upper_;
  const Matrix<Field>& lower_;
  const Matrix<VariableState>& init_states_;
  const Matrix<VariableState>& last_states_;

  static size_t get_dump_id() {
    return std::chrono::system_clock::now().time_since_epoch() /
           std::chrono::milliseconds(1);
  }

  template <typename T>
  static std::string dump_initializer_list(std::span<T> values) {
    std::string result;

    for (size_t i = 0; i < values.size(); ++i) {
      result += std::to_string(values[i]);

      if (i + 1 != values.size()) {
        result += ", ";
      }
    }

    return "{" + result + "}";
  }

  static std::string dump_field(Field value) {
    auto bytes = reinterpret_cast<uint8_t*>(&value);
    return dump_initializer_list({bytes, sizeof(Field)});
  }

  static std::vector<std::string> states_to_string(
      const std::vector<VariableState>& states) {
    std::vector<std::string> result;

    for (auto state : states) {
      if (state == VariableState::BASIC) {
        result.emplace_back("VariableState::BASIC");
      } else if (state == VariableState::AT_LOWER) {
        result.emplace_back("VariableState::AT_LOWER");
      } else {
        result.emplace_back("VariableState::AT_UPPER");
      }
    }

    return result;
  }

 public:
  SimplexCoreDump(const CSCMatrix<Field>& a, const Matrix<Field>& b,
                  const Matrix<Field>& c, const Matrix<Field>& upper,
                  const Matrix<Field>& lower,
                  const Matrix<VariableState>& init_states,
                  const Matrix<VariableState>& last_states)
      : A_(a),
        b_(b),
        c_(c),
        upper_(upper),
        lower_(lower),
        init_states_(init_states),
        last_states_(last_states) {}

  void dump() const {
    size_t dump_id = get_dump_id();
    std::string dump_name = std::format("simplex_core_dump_{}.h", dump_id);

    std::ofstream os(dump_name);

    os << "namespace SimplexDump_" << dump_id << " {\n";

    os << "Matrix<Field> A = {" << linalg::to_dense(A_) << "};\n";
    os << "Matrix<Field> b = {" << b_ << "};\n";
    os << "Matrix<Field> c = {" << c_ << "};\n";

    os << "std::vector<Field> upper = " << dump_initializer_list(upper_)
       << ";\n";

    os << "std::vector<Field> lower = " << dump_initializer_list(lower_)
       << ";\n";

    os << "std::vector<VariableState> last_states = "
       << dump_initializer_list(states_to_string(last_states_)) << ";\n";

    os << "std::vector<VariableState> init_states = "
       << dump_initializer_list(states_to_string(init_states_)) << ";\n";

    os << "}\n";
    os << std::flush;

    std::println("Dumped simplex state into {}.", dump_name);
  }
};

}  // namespace simplex
