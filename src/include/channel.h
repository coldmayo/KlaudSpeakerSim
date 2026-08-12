#ifndef channel_h_INCLUDED
#define channel_h_INCLUDED

#include <string>
#include <vector>
#include <complex>
#include "json.hpp"

struct Channel {
    std::string type;
    bool enclosure_loaded = false;
    bool is_ribbon = false;
    double V_in = 2.83;
    nlohmann::json crossover_cfg = nlohmann::json::object();
    std::vector<std::complex<double>> H_driver;
    std::vector<std::complex<double>> H_crossover;
    std::vector<std::complex<double>> H_branch;
    std::vector<std::complex<double>> Z_zma;
    std::vector<std::complex<double>> Z_in;
    std::vector<double> radii;
    std::vector<double> offsets;
    std::vector<double> power_dissipation;
};

#endif // channel_h_INCLUDED
