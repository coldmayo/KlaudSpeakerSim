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
    const std::vector<std::complex<double>>& H_woofer_branch,
    const std::vector<std::complex<double>>& H_tweeter_branch,
    const std::vector<double>& woofer_a,
    const std::vector<double>& tweeter_a,
    const std::vector<double>& woofer_offset_m,
    const std::vector<double>& tweeter_offset_m,
    double theta_rad) {
    std::vector<std::complex<double>> H_theta(s.size());

    for (size_t i = 0; i < s.size(); ++i) {
        std::complex<double> H_spatial_w(0.0, 0.0);

        if (!woofer_a.empty()) {
            for (size_t j = 0; j < woofer_a.size(); ++j) {
                double directivity = piston_dir(freqs[i], theta_rad, woofer_a[j]);
                std::complex<double> phase_rot = angular_delay(s[i], woofer_offset_m[j], theta_rad);
                H_spatial_w += directivity * phase_rot;
            }
            H_spatial_w /= static_cast<double>(woofer_a.size());
        } else {
            H_spatial_w = 1.0;
        }

        std::complex<double> H_spatial_t(0.0, 0.0);

        if (!tweeter_a.empty()) {
            for (size_t j = 0; j < tweeter_a.size(); ++j) {
                double directivity = ribbon_dir(freqs[i], theta_rad, tweeter_a[j]);
                std::complex<double> phase_rot = angular_delay(s[i], tweeter_offset_m[j], theta_rad);
                H_spatial_t += directivity * phase_rot;
            }
            H_spatial_t /= static_cast<double>(tweeter_a.size());
        } else {
            H_spatial_t = 1.0;
        }

        H_theta[i] = (H_woofer_branch[i] * H_spatial_w) + (H_tweeter_branch[i] * H_spatial_t);
    }

    std::cerr << "woofer a=" << woofer_a[0]
    << " f=500 D=" << piston_dir(500.0, PI/2, woofer_a[0])
    << " f=5000 D=" << piston_dir(5000.0, PI/2, woofer_a[0]) << "\n";
    std::cerr << "tweeter a=" << tweeter_a[0]
    << " f=2000 D=" << piston_dir(2000.0, PI/2, tweeter_a[0])
    << " f=18000 D=" << piston_dir(18000.0, PI/2, tweeter_a[0]) << "\n";

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
