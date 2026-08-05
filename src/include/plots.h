#ifndef plots_h_INCLUDED
#define plots_h_INCLUDED

#include <vector>
#include <complex>
#include <string>

void plot_t_funcs(const std::vector<std::complex<double>>& H, const std::vector<std::complex<double>>& H_woofer, const std::vector<std::complex<double>>& H_tweeter, const std::vector<std::complex<double>>& H_woofer_filt, const std::vector<std::complex<double>>& H_tweeter_filt, const std::vector<std::complex<double>>& s,const std::string& title);
void plot_imp(std::vector<std::complex<double>> Z, const std::vector<std::complex<double>>&s, std::string title);
void plot_filter(const std::vector<std::complex<double>>& H_HP, const std::vector<std::complex<double>>& H_LP, const std::vector<std::complex<double>>& s, std::string title);
void plot_power(const std::vector<double>& P_woofer, const std::vector<double>& P_tweeter, const std::vector<std::complex<double>>& s, std::string title);

#endif // plots_h_INCLUDED
