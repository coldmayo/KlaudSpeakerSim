#ifndef directivity_h_INCLUDED
#define directivity_h_INCLUDED

#include <cmath>
#include <iostream>
#include <vector>
#include <complex>
#include <string>
#include <fstream>

std::vector<std::complex<double>> vented_transfer(std::vector<std::complex<double>> s, struct a_coe a, double omega_0);
void export_sonogram_data(const std::string& filename, const std::vector<double>& freqs, const std::vector<double>& angles_deg, const std::vector<std::vector<std::complex<double>>>& sonogram_matrix);
std::vector<std::complex<double>> polar_response(const std::vector<std::complex<double>>& s,
    const std::vector<double>& freqs, const std::vector<std::complex<double>>& H_woofer_branch, const std::vector<std::complex<double>>& H_tweeter_branch, const std::vector<double>& woofer_a, const std::vector<double>& tweeter_a, const std::vector<double>& woofer_offset_m, const std::vector<double>& tweeter_offset_m, double theta_rad);
#endif // directivity_h_INCLUDED
