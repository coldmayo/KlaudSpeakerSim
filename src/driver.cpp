#include "include/driver.h"
#include "include/parsing.h"

constexpr double DEG2RAD = 3.14159265358979323846 / 180.0;

double interp(double x, double x0, double x1, double y0, double y1) {
    return y0 + (x - x0) * (y1 - y0) / (x1 - x0);
}

//std::vector<std::complex<double>> driver_theoretical() {
    
//}

std::vector<std::complex<double>> driver_measured(const std::string frd_file, const std::vector<double> sim_freqs) {
    frd raw_data = parsing_frd(frd_file);
    std::vector<std::complex<double>> H_driver(sim_freqs.size());
    if (raw_data.freq.empty()) {
        std::cerr << "Warning: Empty FRD data loaded from " << frd_file << std::endl;
        return H_driver;
    }

    int raw_size = raw_data.freq.size();
    int raw_i = 0;
    for (size_t i = 0; i < sim_freqs.size(); ++i) {
        double f = sim_freqs[i];

        // Clamp boundaries
        if (f <= raw_data.freq.front()) {
            double mag_linear = std::pow(10.0, raw_data.spl.front() / 20.0);
            double rad = raw_data.phase.front() * DEG2RAD;
            H_driver[i] = std::polar(mag_linear, rad);
            continue;
        }
        if (f >= raw_data.freq.back()) {
            double mag_linear = std::pow(10.0, raw_data.spl.back() / 20.0);
            double rad = raw_data.phase.back() * DEG2RAD;
            H_driver[i] = std::polar(mag_linear, rad);
            continue;
        }

        // advance raw index to bracket current target frequency
        while (raw_i < raw_size - 1 && raw_data.freq[raw_i + 1] < f) {
            raw_i++;
        }

        // Interpolate
        double f0 = raw_data.freq[raw_i];
        double f1 = raw_data.freq[raw_i + 1];

        double spl = interp(f, f0, f1, raw_data.spl[raw_i], raw_data.spl[raw_i + 1]);
        double phase_deg = interp(f, f0, f1, raw_data.phase[raw_i], raw_data.phase[raw_i + 1]);

        // Convert dB to linear magnitude
        double mag_linear = std::pow(10.0, spl / 20.0);
        double phase_rad = phase_deg * DEG2RAD;
        H_driver[i] = std::polar(mag_linear, phase_rad);
    }
    return H_driver;
}

std::vector<std::complex<double>> delay_transfer(std::vector<std::complex<double>> s, double d) {
    const double c = 343.0;
    double delt = d/c;

    std::vector<std::complex<double>> H(s.size());
    for (size_t i = 0; i < s.size(); ++i) {
        H[i] = std::exp(-s[i]*delt);
    }

    return H;
}
