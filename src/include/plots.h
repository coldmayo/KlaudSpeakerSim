#ifndef plots_h_INCLUDED
#define plots_h_INCLUDED

#include <vector>
#include <complex>
#include <string>

#include "channel.h"

void plot_phase(const std::vector<std::complex<double>>& H, const std::vector<std::complex<double>>& s, std::string title);
void plot_t_funcs(const std::vector<std::complex<double>>& H_system_total, const std::map<std::string, Channel>& channels, const std::vector<std::complex<double>>& s, const std::string& title);
void plot_imp(std::vector<std::complex<double>> Z, const std::vector<std::complex<double>>& s, std::string title);
void plot_filter(const std::map<std::string, Channel>& channels, const std::vector<std::complex<double>>& s, const std::string& title);
void plot_power(const std::map<std::string, Channel>& channels, const std::vector<std::complex<double>>& s, const std::string& title);

#endif // plots_h_INCLUDED
