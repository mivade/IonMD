#ifndef UTIL_HPP
#define UTIL_HPP

#include <iostream>
#include <ctime>
#include <string>
#include <cmath>
#include <armadillo>

using std::sqrt;
using arma::vec;

namespace ionmd {

inline auto timestamp_str() -> std::string {
    const auto now = std::time(nullptr);
    char buff[256];
    std::strftime(buff, sizeof(buff), "%FT%T%z", std::localtime(&now));
    return std::string(buff);
}

}  // namespace ionmd

#endif
