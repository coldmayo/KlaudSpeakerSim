#ifndef crossover_passive_h_INCLUDED
#define crossover_passive_h_INCLUDED

#include <vector>
#include <cmath>
#include <iostream>
#include <complex>
#include "json.hpp"

using json = nlohmann::json;

struct ABCD {
    std::complex<double> A, B, C, D;
};

std::vector<std::complex<double>> pass_cross_abcd(const std::vector<std::complex<double>>& s, const std::vector<std::complex<double>>& Z_driver_zma, const json& branch_network);

#endif // crossover_passive_h_INCLUDED
