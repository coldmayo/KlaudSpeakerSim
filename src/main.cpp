#include <math.h>
#include "include/enclosure.h"
#include "include/plots.h"

int main(int argc, char * argv[]) {

    std::vector<double> freqs = logspace(1, 3, 1000, 10);
    std::cout << "logspace done\n";
    auto s = init_s(freqs);
    std::cout << "init_s done\n";

    // Driver base resonance (e.g., 35 Hz)
    double fs = 35.0;
    double omega_0 = 2.0 * PI * fs;

    // Sealed Box (Butterworth B2)
    double Qtc_sealed = 0.707;
    auto H_sealed = sealed_transfer(s, omega_0, Qtc_sealed);
    std::cout << "sealed done\n";
    plot_t_funcs(H_sealed, s, "Sealed Box");

    // Vented Box (Butterworth B4 alignment)
    double Qts = 0.40;
    double h = 1.00;   // tuning ratio fb/fs
    double alpha = 1.00;   // volume ratio Vas/Vb
    double QL = 7.0;   // enclosure leakage

    a_coe a = gen_a_coes(Qts, QL, h, alpha);
    auto H_vented = vented_transfer(s, a, omega_0);
    std::cout << "vented done\n";
    plot_t_funcs(H_vented, s, "Vented Box");
}
