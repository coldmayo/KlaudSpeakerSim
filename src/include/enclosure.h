#ifndef dsp_h_INCLUDED
#define dsp_h_INCLUDED

#include <vector>
#include <cmath>
#include <iostream>
#include <complex>

struct a_coe {
    double a3;
    double a2;
    double a1;
};

constexpr double PI = 3.14159265358979323846;

std::vector<double> logspace(double start, double end, int n, int base);

std::vector<std::complex<double>> init_s(std::vector<double> freqs);
a_coe gen_a_coes(double Qtc, double QL, double h, double alpha);
std::vector<std::complex<double>> sealed_transfer(std::vector<std::complex<double>> s, double omega_0, double Qtc);

std::vector<std::complex<double>> vented_transfer(std::vector<std::complex<double>> s, struct a_coe a, double omega_0);
#endif // dsp_h_INCLUDED
