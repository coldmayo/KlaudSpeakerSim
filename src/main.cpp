#include <math.h>
#include "include/enclosure.h"
#include "include/plots.h"
#include "include/driver.h"

int main(int argc, char * argv[]) {

    std::vector<double> freqs = logspace(1, 3, 1000, 10);
    std::cout << "logspace done\n";
    auto s = init_s(freqs);
    std::cout << "init_s done\n";

    auto H_driver_bass = driver_measured("driverFiles/Audax AM100G2.frd", freqs);
    auto H_driver_tweeter = driver_measured("driverFiles/Fountek JP-3.frd", freqs);
    auto H_delay = delay_transfer(s, 0.02);

    double fs = 35.0;
    double omega_0 = 2.0 * PI * fs;

    // Sealed Box (Butterworth B2)
    double Qtc_sealed = 0.707;
    auto H_sealed = sealed_transfer(s, omega_0, Qtc_sealed);
    std::cout << "sealed done\n";
    // plot_t_funcs(H_sealed, s, "Sealed Box");

    std::vector<std::complex<double>> H_system_total_sealed(freqs.size());
    for (int i = 0; i < freqs.size(); ++i) {
        H_system_total_sealed[i] = (H_driver_bass[i] + H_driver_tweeter[i]) * H_sealed[i] * H_delay[i];
    }
    plot_t_funcs(H_system_total_sealed, s, "Sealed Box");

    // Vented Box (Butterworth B4 alignment)
    double Qts = 0.40;
    double h = 1.00;   // tuning ratio fb/fs
    double alpha = 1.00;   // volume ratio Vas/Vb
    double QL = 7.0;   // enclosure leakage

    a_coe a = gen_a_coes(Qts, QL, h, alpha);
    auto H_vented = vented_transfer(s, a, omega_0);
    std::cout << "vented done\n";
    //plot_t_funcs(H_vented, s, "Vented Box");

    std::vector<std::complex<double>> H_system_total_vented(freqs.size());
    for (int i = 0; i < freqs.size(); ++i) {
        H_system_total_vented[i] = (H_driver_bass[i] + H_driver_tweeter[i]) * H_vented[i] * H_delay[i];
    }
    plot_t_funcs(H_system_total_vented, s, "Vented Box");
}
