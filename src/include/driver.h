#ifndef driver_h_INCLUDED
#define driver_h_INCLUDED

#include <cmath>
#include <iostream>
#include <vector>
#include <complex>
#include <string>

std::vector<std::complex<double>> delay_transfer(std::vector<std::complex<double>> s, double d);
std::vector<std::complex<double>> driver_measured(const std::string frd_file, const std::vector<double> sim_freqs);

#endif // driver_h_INCLUDED
