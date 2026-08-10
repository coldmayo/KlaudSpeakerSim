#include "matplotlibcpp.h"
#include <complex>
#include <cmath>
#include <vector>
#include <string>

namespace plt = matplotlibcpp;

constexpr double PI = 3.14159265358979323846;

std::vector<double> phase_angle(const std::vector<std::complex<double>> & H) {
    std::vector<double> phi(H.size());

    for (int i = 0; i < H.size(); ++i) {
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

	for (int i = 0; i < s.size(); ++i) {
    	omega[i] = s[i].imag();
	}
	
	auto tau_p = phase_delay(phi_uw, omega);
	auto tau_g = group_delay(phi_uw, omega);

	std::vector<double> freq(s.size());
	for (size_t i = 0; i < s.size(); ++i) {
		freq[i] = std::imag(s[i]) / (2.0 * M_PI);
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

void plot_power(const std::vector<double>& P_woofer, const std::vector<double>& P_tweeter, const std::vector<std::complex<double>>& s, std::string title) {
    std::vector<double> freq(s.size());
    std::vector<double> P_total(s.size());

    for (size_t i = 0; i < s.size(); ++i) {
        double f = std::imag(s[i]) / (2.0 * M_PI);
        freq[i] = f;
        P_total[i] = P_woofer[i] + P_tweeter[i];
    }

    plt::figure_size(1200, 780);

    plt::named_semilogx("Woofer Power", freq, P_woofer);
    plt::named_semilogx("Tweeter Power", freq, P_tweeter);
    plt::named_semilogx("Total System Power", freq, P_total);

    plt::grid(true);

    plt::xlabel("Frequency (Hz)");
    plt::ylabel("Active Power (W)");

    plt::legend();

    plt::save(title + ".png");
    plt::close();
}


