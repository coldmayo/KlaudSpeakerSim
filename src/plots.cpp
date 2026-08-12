#include "matplotlibcpp.h"
#include "include/channel.h"
#include "include/plots.h"
#include <complex>
#include <cmath>
#include <vector>
#include <string>
#include <map>

namespace plt = matplotlibcpp;

constexpr double PI = 3.14159265358979323846;

std::vector<double> phase_angle(const std::vector<std::complex<double>> & H) {
    std::vector<double> phi(H.size());

    for (size_t i = 0; i < H.size(); ++i) {
        phi[i] = std::atan2(H[i].imag(), H[i].real());
    }
    return phi;
}

std::vector<double> unwrap_phase(const std::vector<double>& phi) {
    std::vector<double> unwrapped(phi.size());
    if (phi.empty()) return unwrapped;

    unwrapped[0] = phi[0];
    double offset = 0.0;
    for (size_t i = 1; i < phi.size(); ++i) {
        double delta = phi[i] - phi[i - 1];
        if (delta > PI)       offset -= 2.0 * PI;
        else if (delta < -PI) offset += 2.0 * PI;
        unwrapped[i] = phi[i] + offset;
    }
    return unwrapped;
}

std::vector<double> phase_delay(std::vector<double> p_a, std::vector<double> omega) {
    std::vector<double> tau_p(p_a.size());
    for (size_t i = 0; i < p_a.size(); ++i) {
        tau_p[i] = (omega[i] > 0.0) ? -p_a[i] / omega[i] : 0.0;
    }
    return tau_p;
}

std::vector<double> group_delay(const std::vector<double>& phi_unwrapped, const std::vector<double>& omega) {
    size_t N = phi_unwrapped.size();
    std::vector<double> tau_g(N, 0.0);
    if (N < 2) return tau_g;

    tau_g[0] = -(phi_unwrapped[1] - phi_unwrapped[0]) / (omega[1] - omega[0]);

    for (size_t i = 1; i < N - 1; ++i) {
        tau_g[i] = -(phi_unwrapped[i + 1] - phi_unwrapped[i - 1])
                   / (omega[i + 1] - omega[i - 1]);
    }

    tau_g[N - 1] = -(phi_unwrapped[N - 1] - phi_unwrapped[N - 2])
                   / (omega[N - 1] - omega[N - 2]);

    return tau_g;
}

void plot_phase(const std::vector<std::complex<double>>& H, const std::vector<std::complex<double>>& s, std::string title) {
    auto phi = phase_angle(H);
    auto phi_uw = unwrap_phase(phi);
    std::vector<double> omega(s.size());

    for (size_t i = 0; i < s.size(); ++i) {
        omega[i] = s[i].imag();
    }

    auto tau_p = phase_delay(phi_uw, omega);
    auto tau_g = group_delay(phi_uw, omega);

    std::vector<double> freq(s.size());
    for (size_t i = 0; i < s.size(); ++i) {
        freq[i] = std::imag(s[i]) / (2.0 * PI);
    }

    plt::figure_size(1200, 780);

    plt::named_semilogx("Group Delay", freq, tau_g);
    plt::named_semilogx("Phase Delay", freq, tau_p);
    plt::grid(true);

    plt::xlabel("Frequency (Hz)");
    plt::ylabel("Samples");
    plt::title(title);

    plt::legend();

    plt::save(title + ".png");
    plt::close();
}

void plot_t_funcs(const std::vector<std::complex<double>>& H_system_total, const std::map<std::string, Channel>& channels, const std::vector<std::complex<double>>& s, const std::string& title) {
    const size_t N = s.size();
    std::vector<double> freq(N);
    std::vector<double> mag_db_system(N);

    for (size_t i = 0; i < N; ++i) {
        freq[i] = std::imag(s[i]) / (2.0 * PI);
        mag_db_system[i] = 20.0 * std::log10(std::abs(H_system_total[i]) + 1e-12);
    }

    plt::figure_size(1200, 780);

    for (const auto& [key, ch] : channels) {
        std::vector<double> mag_db_raw(N);
        std::vector<double> mag_db_branch(N);
        for (size_t i = 0; i < N; ++i) {
            mag_db_raw[i] = 20.0 * std::log10(std::abs(ch.H_driver[i]) + 1e-12);
            mag_db_branch[i] = 20.0 * std::log10(std::abs(ch.H_branch[i]) + 1e-12);
        }
        plt::named_semilogx("$H_{" + key + "}$", freq, mag_db_raw);
        plt::named_semilogx("Filtered $H_{" + key + "}$", freq, mag_db_branch);
    }

    plt::named_semilogx("$H_{System}$", freq, mag_db_system);
    plt::grid(true);

    plt::xlabel("Frequency (Hz)");
    plt::ylabel("Magnitude (dB)");
    plt::title(title);

    plt::legend();

    plt::save(title + ".png");
    plt::close();
}

void plot_imp(std::vector<std::complex<double>> Z, const std::vector<std::complex<double>>& s, std::string title) {
    const size_t N = s.size();
    std::vector<double> freq(N);
    std::vector<double> mag_Z_system(N);
    for (size_t i = 0; i < N; ++i) {
        freq[i] = std::imag(s[i]) / (2.0 * PI);
        mag_Z_system[i] = 20.0 * std::log10(std::abs(Z[i]) + 1e-12);
    }

    plt::figure_size(1200, 780);

    plt::named_semilogx("$Z_{System}$", freq, mag_Z_system);

    plt::grid(true);

    plt::xlabel("Frequency (Hz)");
    plt::ylabel("Impedence ($\\Omega$)");
    plt::title(title);

    plt::legend();

    plt::save(title + ".png");
    plt::close();
}

// One filter-response trace per channel.
void plot_filter(const std::map<std::string, Channel>& channels, const std::vector<std::complex<double>>& s, const std::string& title) {
    const size_t N = s.size();
    std::vector<double> freq(N);
    for (size_t i = 0; i < N; ++i) {
        freq[i] = std::imag(s[i]) / (2.0 * PI);
    }

    plt::figure_size(1200, 780);

    for (const auto& [key, ch] : channels) {
        std::vector<double> mag_db(N);
        for (size_t i = 0; i < N; ++i) {
            mag_db[i] = 20.0 * std::log10(std::abs(ch.H_crossover[i]) + 1e-12);
        }
        plt::named_semilogx(key + " Filter", freq, mag_db);
    }

    plt::grid(true);

    plt::xlabel("Frequency (Hz)");
    plt::ylabel("Transfer Function (dB)");
    plt::title(title);

    plt::legend();

    plt::save(title + ".png");
    plt::close();
}

// One power trace per channel, plus total system power.
void plot_power(const std::map<std::string, Channel>& channels, const std::vector<std::complex<double>>& s, const std::string& title) {
    const size_t N = s.size();
    std::vector<double> freq(N);
    std::vector<double> P_total(N, 0.0);

    for (size_t i = 0; i < N; ++i) {
        freq[i] = std::imag(s[i]) / (2.0 * PI);
    }

    plt::figure_size(1200, 780);

    for (const auto& [key, ch] : channels) {
        plt::named_semilogx(key + " Power", freq, ch.power_dissipation);
        for (size_t i = 0; i < N && i < ch.power_dissipation.size(); ++i) {
            P_total[i] += ch.power_dissipation[i];
        }
    }

    plt::named_semilogx("Total System Power", freq, P_total);

    plt::grid(true);

    plt::xlabel("Frequency (Hz)");
    plt::ylabel("Active Power (W)");
    plt::title(title);

    plt::legend();

    plt::save(title + ".png");
    plt::close();
}
