
#include "include/enclosure.h"

// put sequence into logspace
std::vector<double> logspace(double start, double end, int n, int base = 10) {
    std::vector<double> freqs(n); 
    const double step_size = (end-start) / (n-1);

    for (int i = 0; i < n; i++) {
        freqs[i] = std::pow(base, start + (i * step_size));
    }
    return freqs;
}

// make a logified s range
std::vector<std::complex<double>> init_s(std::vector<double> freqs) {
    std::vector<std::complex<double>> s(freqs.size());
    for (int i = 0; i < freqs.size(); ++i) {
        double omega = 2.0 * PI * freqs[i];
        s[i] = std::complex<double>(0.0, omega);
    }
    return s;
}

// second order high pass filter in s domain
std::vector<std::complex<double>> sealed_transfer(std::vector<std::complex<double>> s, double omega_0, double Qtc) {
    std::vector<std::complex<double>> H(s.size());
    const double omega2_0 = omega_0*omega_0;
    const double damping = omega_0/Qtc;
    for (int i = 0; i < s.size(); i++) {
        std::complex<double> s2 = s[i]*s[i];
        H[i] = s2/(s2 + (damping)*s[i] + omega2_0);
    }
    return H;
}

// find the denominator coefficents for vented transfer function
a_coe gen_a_coes(double Qtc, double QL, double h, double alpha) {
    a_coe a;
    a.a3 = (1/Qtc) + (h/QL);
    a.a2 = 1 + (h*h) + alpha + (h/(Qtc*QL));
    a.a1 = (h*h/Qtc) + (h/QL);
    return a;
}

std::vector<std::complex<double>> vented_transfer(std::vector<std::complex<double>> s, struct a_coe a, double omega_0) {
    std::vector<std::complex<double>> H(s.size());
    const double w0 = omega_0;
    const double w0_2 = w0 * w0;
    const double w0_3 = w0_2 * w0;
    const double w0_4 = w0_3 * w0;
    for (int i = 0; i < s.size(); ++i) {
        std::complex<double> si = s[i];
        std::complex<double> s2 = si * si;
        std::complex<double> s3 = s2 * si;
        std::complex<double> s4 = s2 * s2;

        std::complex<double> denom = s4 + (a.a3 * w0 * s3) + (a.a2 * w0_2 * s2) + (a.a1 * w0_3 * si) + w0_4;
        H[i] = s4 / denom;
    }
    return H;
}
