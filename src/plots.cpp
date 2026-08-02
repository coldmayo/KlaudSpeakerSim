#include "matplotlibcpp.h"
#include <complex>
#include <cmath>
#include <vector>
#include <string>

namespace plt = matplotlibcpp;

void plot_t_funcs(const std::vector<std::complex<double>>& H, const std::vector<std::complex<double>>& H_woofer, const std::vector<std::complex<double>>& H_tweeter, const std::vector<std::complex<double>>& H_woofer_filt, const std::vector<std::complex<double>>& H_tweeter_filt, const std::vector<std::complex<double>>& s,const std::string& title) {
    std::vector<double> freq(s.size());
    std::vector<double> mag_db_system(H.size());
    std::vector<double> mag_db_woofer(H.size());
    std::vector<double> mag_db_tweeter(H.size());
    std::vector<double> mag_db_woofer_f(H.size());
    std::vector<double> mag_db_tweeter_f(H.size());

    for (size_t i = 0; i < s.size(); ++i) {
        double f = std::imag(s[i]) / (2.0 * M_PI);
        freq.push_back(f);
        mag_db_system.push_back(20.0 * std::log10(std::abs(H[i])));
        mag_db_tweeter.push_back(20.0 * std::log10(std::abs(H_tweeter[i])));
        mag_db_tweeter_f.push_back(20.0 * std::log10(std::abs(H_tweeter_filt[i])));
        mag_db_woofer.push_back(20.0 * std::log10(std::abs(H_woofer[i])));
        mag_db_woofer_f.push_back(20.0 * std::log10(std::abs(H_woofer_filt[i])));
    }

    plt::figure_size(1200, 780);

    plt::named_semilogx("$H_{Tweeter(s)}$", freq, mag_db_tweeter);
    plt::named_semilogx("$H_{Woofer(s)}$", freq, mag_db_woofer);
    plt::named_semilogx("Filtered $H_{Tweeter(s)}$", freq, mag_db_tweeter_f);
    plt::named_semilogx("Filtered $H_{Woofer(s)}$", freq, mag_db_woofer_f);
    plt::named_semilogx("$H_{System}$", freq, mag_db_system);
    plt::grid(true);

    plt::xlabel("Frequency (Hz)");
    plt::ylabel("Magnitude (dB)");
    plt::title(title);

    plt::legend();

    plt::save(title + ".png");
    plt::close();
}
