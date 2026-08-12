#ifndef directivity_h_INCLUDED
#define directivity_h_INCLUDED

#include <cmath>
#include <iostream>
#include <vector>
#include <complex>
#include <string>
#include <fstream>
#include "channel.h"

double piston_dir(double f, double theta_rad, double a);
double ribbon_dir(double f, double theta_rad, double half_length);
std::complex<double> angular_delay(std::complex<double> s_val, double offset_m, double theta_rad);

std::vector<std::complex<double>> polar_response(const std::vector<std::complex<double>>& s, const std::vector<double>& freqs, const std::map<std::string, Channel>& channels, double theta_rad);

void export_sonogram_data(const std::string& filename, const std::vector<double>& freqs, const std::vector<double>& angles_deg, const std::vector<std::vector<std::complex<double>>>& sonogram_matrix);

#endif // directivity_h_INCLUDED
