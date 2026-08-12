#include "include/directivity.h"

constexpr double PI = 3.14159265358979323846;

double piston_dir(double f, double theta_rad, double a) {
    if (std::abs(theta_rad) < 1e-6 || a <= 0.0 || f <= 0.0) {
        return 1.0;
    }
    const double c = 343.0;
    double k = 2.0 * PI * f / c;
    double x = std::abs(k * a * std::sin(theta_rad));
    if (std::abs(x) < 1e-9) return 1.0;   // avoid 0/0 at theta=0
    return 2.0 * std::cyl_bessel_j(1, x) / x;
}

double ribbon_dir(double f, double theta_rad, double half_length) {
    if (std::abs(theta_rad) < 1e-6 || half_length <= 0.0 || f <= 0.0) {
        return 1.0;
    }
    const double c = 343.0;
    double k = 2.0 * PI * f / c;
    double x = k * half_length * std::sin(theta_rad);
    if (std::abs(x) < 1e-9) return 1.0;   // avoid 0/0 at theta=0
    return std::sin(x) / x;

}

std::complex<double> angular_delay(std::complex<double> s_val, double offset_m, double theta_rad) {
    constexpr double c = 343.0; // speed of sound in m/s
    double delta_d = offset_m * std::sin(theta_rad);
    double delta_t = delta_d / c;
    return std::exp(-s_val * delta_t);
}

std::vector<std::complex<double>> polar_response(
    const std::vector<std::complex<double>>& s,
    const std::vector<double>& freqs,
    const std::map<std::string, Channel>& channels,
    double theta_rad) {

    std::vector<std::complex<double>> H_theta(s.size(), std::complex<double>(0.0, 0.0));

    for (size_t i = 0; i < s.size(); ++i) {
        std::complex<double> H_sum(0.0, 0.0);

        for (const auto& [key, ch] : channels) {
            std::complex<double> H_spatial(0.0, 0.0);

            if (!ch.radii.empty()) {
                for (size_t j = 0; j < ch.radii.size(); ++j) {
                    double directivity = ch.is_ribbon
                    ? ribbon_dir(freqs[i], theta_rad, ch.radii[j])
                    : piston_dir(freqs[i], theta_rad, ch.radii[j]);
                    std::complex<double> phase_rot = angular_delay(s[i], ch.offsets[j], theta_rad);
                    H_spatial += directivity * phase_rot;
                }
                H_spatial /= static_cast<double>(ch.radii.size());
            } else {
                H_spatial = 1.0;
            }

            H_sum += ch.H_branch[i] * H_spatial;
        }

        H_theta[i] = H_sum;
    }

    return H_theta;
}

void export_sonogram_data(const std::string& filename, const std::vector<double>& freqs, const std::vector<double>& angles_deg, const std::vector<std::vector<std::complex<double>>>& sonogram_matrix) {
    std::ofstream out(filename);

    // Header
    out << "Angle";
    for (double f : freqs) {
        out << "," << f;
    }
    out << "\n";

    // Rows
    for (size_t a = 0; a < angles_deg.size(); ++a) {
        out << angles_deg[a];
        for (size_t f = 0; f < freqs.size(); ++f) {
            double mag_db = 20.0 * std::log10(std::abs(sonogram_matrix[a][f]) + 1e-12);
            out << "," << mag_db;
        }
        out << "\n";
    }
}
