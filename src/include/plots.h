#ifndef plots_h_INCLUDED
#define plots_h_INCLUDED

#include <vector>
#include <complex>
#include <string>

void plot_t_funcs(const std::vector<std::complex<double>>& H, const std::vector<std::complex<double>>& H_woofer, const std::vector<std::complex<double>>& H_tweeter, const std::vector<std::complex<double>>& H_woofer_filt, const std::vector<std::complex<double>>& H_tweeter_filt, const std::vector<std::complex<double>>& s,const std::string& title);

#endif // plots_h_INCLUDED
