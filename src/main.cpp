#include <iostream>
#include <fstream>
#include <vector>
#include <complex>
#include <cmath>
#include "json.hpp"

#include "include/enclosure.h"
#include "include/plots.h"
#include "include/driver.h"
#include "include/crossover.h"

using json = nlohmann::json;

int main(int argc, char * argv[]) {
    std::string config_file = "configs/speaker1.json";
    if (argc >= 2) {
        config_file = argv[1];
    }
    std::ifstream file(config_file);
    if (!file.is_open()) {
        throw std::runtime_error("Failed to open config file: " + config_file);
    }

    json json_f;
    file >> json_f;
    
    std::vector<double> freqs = logspace(1, 3, 1000, 10);
    auto s = init_s(freqs);

    std::vector<std::complex<double>> H_driver_woofer(s.size(), 0.0);
    std::vector<std::complex<double>> H_driver_tweeter(s.size(), 0.0);

    std::string type;
    std::string frd_file;
    double d;
    for (const auto& d_json : json_f["drivers"]) {
        type = d_json["type"];
        frd_file = d_json["frd_path"];
        d = d_json["physical_delay_m"];

        auto H_raw = driver_measured(frd_file, freqs);
        auto H_delay = delay_transfer(s, d);

        for (size_t j = 0; j < s.size(); ++j) {
            std::complex<double> H_delayed_driver = H_raw[j] * H_delay[j];
            if (type == "woofer") {
                H_driver_woofer[j] += H_delayed_driver;
            } else if (type == "tweeter") {
                H_driver_tweeter[j] += H_delayed_driver;
            }
        }
    }

    double fs = json_f["fs_hz"];
    double fc = json_f["crossover"]["crossover_freq_hz"];
    double omega_0 = 2.0 * PI * fs;
    double omega_c = 2.0 * PI * fc;
    std::vector<std::complex<double>> H_box(s.size(), 1.0);

    // Sealed Box (Butterworth B2)
    if (json_f["box_type"] == "sealed") {
        double Qtc = json_f["sealed_params"]["Qtc"];
        H_box = sealed_transfer(s, omega_0, Qtc);
    } else if (json_f["box_type"] == "vented") {
        double Qts = json_f["vented_params"]["Qts"];
        double h = json_f["vented_params"]["tuning_ratio_h"];
        double alpha = json_f["vented_params"]["volume_ratio_alpha"];
        double QL = json_f["vented_params"]["leakage_QL"];

        a_coe a = gen_a_coes(Qts, QL, h, alpha);
        H_box = vented_transfer(s, a, omega_0);
    }

    auto H_HP = select_crossover(json_f["crossover"]["type"], json_f["crossover"]["order"], s, omega_c, 0);
    auto H_LP = select_crossover(json_f["crossover"]["type"], json_f["crossover"]["order"], s, omega_c, 1);

    std::vector<std::complex<double>> H_system_total(freqs.size());
    for (size_t i = 0; i < freqs.size(); ++i) {
        std::complex<double> woofer_branch  = H_driver_woofer[i] * H_box[i] * H_LP[i];
        std::complex<double> tweeter_branch = H_driver_tweeter[i] * H_HP[i];

        H_system_total[i] = woofer_branch + tweeter_branch;
    }

    plot_t_funcs(H_system_total, s, json_f["box_type"] == "sealed" ? "Sealed Box System" : "Vented Box System");
    return 0;

}
