#include "matplotlibcpp.h"
#include <complex>
#include <cmath>
#include <vector>
#include <string>

namespace plt = matplotlibcpp;

void plot_t_funcs(const std::vector<std::complex<double>>& H, const std::vector<std::complex<double>>& s,const std::string& title) {
    std::vector<double> freq;
    std::vector<double> mag_db;

    freq.reserve(s.size());
    mag_db.reserve(H.size());

    for (size_t i = 0; i < s.size(); ++i) {
        double f = std::imag(s[i]) / (2.0 * M_PI);
        freq.push_back(f);
        mag_db.push_back(20.0 * std::log10(std::abs(H[i])));
    }

    plt::figure_size(1200, 780);

    plt::named_semilogx(title, freq, mag_db);
    plt::grid(true);

    plt::xlabel("Frequency (Hz)");
    plt::ylabel("Magnitude (dB)");
    plt::title(title);

    plt::legend();

    plt::save(title + ".png");
    plt::close();
}
