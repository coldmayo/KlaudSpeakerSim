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

void plot_imp(std::vector<std::complex<double>> Z, const std::vector<std::complex<double>>&s, std::string title) {
    std::vector<double> freq(s.size());
    std::vector<double> mag_Z_system(Z.size());
    for (size_t i = 0; i < s.size(); ++i) {
        double f = std::imag(s[i]) / (2.0 * M_PI);
        freq.push_back(f);
        mag_Z_system.push_back(20.0 * std::log10(std::abs(Z[i])));
    }

    plt::figure_size(1200, 780);

    plt::named_semilogx("$Z_{System}$", freq, mag_Z_system);

    plt::grid(true);
    
    plt::xlabel("Frequency (Hz)");
    plt::ylabel("Impedence ($\Omega$)");

    plt::legend();

    plt::save(title + ".png");
    plt::close();
}

void plot_filter(const std::vector<std::complex<double>>& H_HP, const std::vector<std::complex<double>>& H_LP, const std::vector<std::complex<double>>& s, std::string title) {
    std::vector<double> freq(s.size());
    std::vector<double> mag_LP_db(s.size());
    std::vector<double> mag_HP_db(s.size());

    for (size_t i = 0; i < s.size(); ++i) {
        double f = std::imag(s[i]) / (2.0 * M_PI);
        freq[i] = f;

        mag_LP_db[i] = 20.0 * std::log10(std::abs(H_LP[i]) + 1e-12);
        mag_HP_db[i] = 20.0 * std::log10(std::abs(H_HP[i]) + 1e-12);
    }

    plt::figure_size(1200, 780);

    plt::named_semilogx("Woofer Filter (LP)", freq, mag_LP_db);
    plt::named_semilogx("Tweeter Filter (HP)", freq, mag_HP_db);

    plt::grid(true);

    plt::xlabel("Frequency (Hz)");
    plt::ylabel("Transfer Function (dB)");

    plt::legend();

    plt::save(title + ".png");
    plt::close();
}
