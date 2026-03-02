#pragma once

#include <filesystem>
#include <string_view>

#ifndef LOG_FILEPATH
#error \
    "Log filepath is not defined! Please, add LOG_FILEPATH compile definition.";
#endif

#ifndef RESOURCES_FILEPATH
#error \
    "Resources filepath is not defined! Please, add RESOURCES_FILEPATH compile definition.";
#endif

namespace paths {

inline std::filesystem::path log(std::string_view filename) {
  return std::filesystem::path(LOG_FILEPATH) / filename;
}

inline std::filesystem::path resource(std::string_view filename) {
  return std::filesystem::path(RESOURCES_FILEPATH) / filename;
}

}  // namespace paths
