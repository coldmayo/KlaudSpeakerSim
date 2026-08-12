#include <iostream>
#include <fstream>
#include <vector>
#include <complex>
#include <cmath>
#include <format>
#include "json.hpp"

#include "include/channel.h"
#include "include/enclosure.h"
#include "include/plots.h"
#include "include/driver.h"
#include "include/crossover.h"
#include "include/directivity.h"
#include "include/crossover_passive.h"
#include "include/baffle.h"
#include "include/parsing.h"
#include "include/Thiele_Small.h"

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
    
    std::vector<double> freqs = logspace(1, 4.301, 1000, 10);
    auto s = init_s(freqs);

    std::map<std::string, Channel> channels;

    double baffle_w = json_f.value("baffle", json::object()).value("width_m", 0.20);
    double baffle_h = json_f.value("baffle", json::object()).value("height_m", 0.35);
    double listener_z = json_f.value("baffle", json::object()).value("listener_dist_m", 1.0);
    double dl_step = json_f.value("baffle", json::object()).value("dl_step_m", 0.002);
    double fs_baffle = json_f.value("baffle", json::object()).value("fs_hz", 96000.0);
    int num_samples = json_f.value("baffle", json::object()).value("num_samples", 2048);

	// Gather data from zma and frd files
    std::string type;
    std::string frd_file;
    std::string zma_file;
    double d;

    for (const auto& d_json : json_f["drivers"]) {
        type = d_json["type"];
        frd_file = d_json["frd_path"];
        zma_file = d_json["zma_path"];
        d = d_json["physical_delay_m"];
        double Sd = d_json.value("eff_surface_area", -1.0);
        double w  = d_json.value("w_ribbon_strip", -1.0);
        double y_offset = d;
        double radius_m = eff_piston_r(Sd, w);

        std::string channel_key = d_json.value("channel", type);
        Channel& ch = channels[channel_key];
        if (ch.H_driver.empty()) {
            ch.H_driver.assign(s.size(), std::complex<double>(0.0, 0.0));
            ch.Z_zma.assign(s.size(), std::complex<double>(8.0, 0.0));
            ch.type = type;
            ch.enclosure_loaded = d_json.value("enclosure_loaded", false);
            ch.V_in = d_json.value("V_in", 2.83);
            ch.crossover_cfg = d_json.value("crossover", json::object());
        }


		source_coords s_c {
            d_json.value("pos_x_m", baffle_w / 2.0),
            d_json.value("pos_y_m", 0.10)
        };

        listener_coords l_c{s_c.x, s_c.y, listener_z};

        auto H_raw = driver_measured(frd_file, freqs);
        auto H_delay = delay_transfer(s, d);
        auto H_raw_fft = b_response(baffle_w, baffle_h, radius_m, dl_step, fs_baffle, num_samples, s_c, l_c);

        double R_direct = std::sqrt((l_c.x - s_c.x)*(l_c.x - s_c.x) + (l_c.y - s_c.y)*(l_c.y - s_c.y) + (l_c.z)*(l_c.z));
        double p_direct_amp = 2.0 / R_direct;
        std::vector<std::complex<double>> H_baffle(s.size(), 1.0);
        for (size_t j = 0; j < freqs.size(); ++j) {
            double f = freqs[j];
            double bin_float = f * num_samples / fs_baffle;
            int bin0 = std::clamp(static_cast<int>(std::floor(bin_float)), 0, num_samples - 1);
            int bin1 = std::clamp(bin0 + 1, 0, num_samples - 1);

            double frac = bin_float - bin0;
            std::complex<double> H_fft_interp = H_raw_fft[bin0] * (1.0 - frac) + H_raw_fft[bin1] * frac;

            double w_k = (2.0 * PI * f) / 343.0;
            double kp = w_k * radius_m;
            double driver_filter = 1.0;
            if (kp > 1e-6) {
                driver_filter = 2.0 * j1(kp) / kp;
            }

            std::complex<double> H_direct(p_direct_amp, 0.0);
            std::complex<double> H_diffracted = (H_fft_interp - H_direct) * driver_filter; // Directivity scaled[cite: 1]
            H_baffle[j] = (H_direct + H_diffracted) / p_direct_amp; // Relative transfer factor
        }

        for (size_t j = 0; j < s.size(); ++j) {
            std::complex<double> H_delayed_driver = H_raw[j] * H_delay[j] * H_baffle[j];
            ch.H_driver[j] += H_delayed_driver;
        }

        ch.radii.push_back(radius_m);
        ch.offsets.push_back(y_offset);

        if (!zma_file.empty()) {
            auto Z_interp = zma_measured(zma_file, freqs);
            ch.Z_zma = Z_interp;
        }

        // Find T/S Parameters
        if (d_json.contains("fund_TS_params")) {
            const auto& ts = d_json["fund_TS_params"];
            if (ts.contains("V_as")) {
                zma_data z_ = zma_interp(parsing_zma(zma_file), 250);
                double M_ms = ts["M_ms"];
                double V_as = ts["V_as"];
                double Cms = C_ms(V_as, Sd);
                double fs_d = f_s(M_ms, Cms);
                double Qts = Q_ts(z_.Q_es, z_.Q_ms);
                double Bl = B_l(fs_d, M_ms, z_.R_e, z_.Q_es);
                double n0 = n_0(fs_d, V_as, z_.Q_es);
                double SPL_1W = SPL1W(n0);
                double SPL_283 = SPL283(SPL_1W, z_.R_e);
                double XL = X_L(z_.Z_high, z_.R_e);
                double Le = L_e(XL, z_.f_high);
                double EBP_ = EBP(fs_d, z_.Q_es);
                std::cout << "T/S Parameters of the " << d_json["name"] << " Driver:\nC_ms: " << Cms << "\nf_s: " << fs_d << "\nQts: " << Qts << "\nBl: " << Bl << "\nQ_es: " << z_.Q_es << "\nQ_m: " << z_.Q_ms << "\nn_0: " << n0 << "\nSPL_1W: " << SPL_1W << "\nSPL_2.83W: " << SPL_283 << "\nX_L: " << XL << "\nL_e: " << Le << "\nEBP: " << EBP_ << "\n\n" << std::endl;
            } else {
                zma_data z_ = zma_interp(parsing_zma(zma_file), 20000);
                double S_d = d_json["eff_surface_area"];
                double XL = X_L(z_.Z_high, z_.R_e);
                double Le = L_e(XL, z_.f_high);
                std::cout << "T/S Parameters of the " << d_json["name"] << " Driver:" << "\nX_L: " << XL << "\nL_e: " << Le << "\n\n" << std::endl;
            }
        }
    }

    double fs = json_f["fs_hz"];
    double omega_0 = 2.0 * PI * fs;
    std::vector<std::complex<double>> H_box(s.size(), 1.0);

    // Get transfer functions describing the enclosure of a loudspeaker
    std::string box_type = json_f.value("box_type", "none");
    if (box_type == "sealed") {
        double Qtc = json_f["sealed_params"]["Qtc"];
        H_box = sealed_transfer(s, omega_0, Qtc);
    } else if (box_type == "vented") {
        double Qts = json_f["vented_params"]["Qts"];
        double h = json_f["vented_params"]["tuning_ratio_h"];
        double alpha = json_f["vented_params"]["volume_ratio_alpha"];
        double QL = json_f["vented_params"]["leakage_QL"];

        a_coe a = gen_a_coes(Qts, QL, h, alpha);
        H_box = vented_transfer(s, a, omega_0);
    }

	// finding the transfer function for active or passive crossover circuits
	for (auto& [key, ch] : channels) {
        const json& cx = ch.crossover_cfg;
        std::string mode = cx.value("mode", "active");
        ch.H_crossover.assign(s.size(), std::complex<double>(1.0, 0.0));
        ch.Z_in = ch.Z_zma; // default (no crossover network): amp/branch sees the driver's own impedance

        if (mode == "active") {
            std::vector<std::complex<double>> H_hp(s.size(), std::complex<double>(1.0, 0.0));
            std::vector<std::complex<double>> H_lp(s.size(), std::complex<double>(1.0, 0.0));
            std::string xtype = cx.value("type", "linkwitz_riley");

            if (cx.contains("hp_freq_hz")) {
                double fc = cx["hp_freq_hz"];
                int order = cx.value("hp_order", 2);
                double omega_c = 2.0 * PI * fc;
                H_hp = select_crossover(xtype, order, s, omega_c, 0);
            }
            if (cx.contains("lp_freq_hz")) {
                double fc = cx["lp_freq_hz"];
                int order = cx.value("lp_order", 2);
                double omega_c = 2.0 * PI * fc;
                H_lp = select_crossover(xtype, order, s, omega_c, 1);
            }
            for (size_t j = 0; j < s.size(); ++j) {
                ch.H_crossover[j] = H_hp[j] * H_lp[j];
            }
            ch.Z_in = ch.Z_zma;
        } else if (mode == "passive") {
            if (cx.contains("branch")) {
                std::tie(ch.H_crossover, ch.Z_in) = pass_cross_abcd(s, ch.Z_zma, cx["branch"]);
            }
        }
    }

	// Finding the systems impulse response
	std::vector<std::complex<double>> Z_system(freqs.size(), std::complex<double>(0.0, 0.0));
    for (size_t i = 0; i < freqs.size(); ++i) {
        std::complex<double> Y_sum(0.0, 0.0);
        for (auto& [key, ch] : channels) {
            if (std::abs(ch.Z_in[i]) > 1e-12) {
                Y_sum += 1.0 / ch.Z_in[i];
            }
        }
        Z_system[i] = (std::abs(Y_sum) > 1e-12) ? (1.0 / Y_sum) : std::complex<double>(0.0, 0.0);
    }


	// finding the total system frequency response
	std::vector<std::complex<double>> H_system_total(freqs.size(), std::complex<double>(0.0, 0.0));
    for (auto& [key, ch] : channels) {
        ch.H_branch.assign(s.size(), std::complex<double>(0.0, 0.0));
        for (size_t i = 0; i < s.size(); ++i) {
            std::complex<double> enclosure_factor = ch.enclosure_loaded ? H_box[i] : std::complex<double>(1.0, 0.0);
            ch.H_branch[i] = ch.H_driver[i] * enclosure_factor * ch.H_crossover[i];
            H_system_total[i] += ch.H_branch[i];
        }
    }


    // Plotting
    std::string title = std::format("{}-way {} System", channels.size(), box_type);

    plot_t_funcs(H_system_total, channels, s, title);
    plot_imp(Z_system, s, "System Impedence");
    plot_filter(channels, s, "Filter Response");

    // Power
    for (auto& [key, ch] : channels) {
        ch.power_dissipation.assign(s.size(), 0.0);
        for (size_t i = 0; i < s.size(); ++i) {
            double re_z = std::real(ch.Z_in[i]);
            ch.power_dissipation[i] = (re_z > 1e-9)
            ? (ch.V_in * ch.V_in * std::norm(ch.H_crossover[i])) / re_z
            : 0.0;
        }
    }
    plot_power(channels, s, "Power Dissipation");

    plot_phase(H_system_total, s, "Group and Phase Delay");

    // Directivity
    std::vector<double> angles_deg;
    for (double deg = -90.0; deg <= 90.0; deg += 2.0) {
        angles_deg.push_back(deg);
    }

    std::vector<std::vector<std::complex<double>>> sonogram_matrix;
    for (double deg : angles_deg) {
        double rad = deg * M_PI / 180.0;
        auto H_theta = polar_response(s, freqs, channels, rad);
        sonogram_matrix.push_back(H_theta);
    }

    export_sonogram_data("sonogram_data.csv", freqs, angles_deg, sonogram_matrix);
    int result = std::system("python3 sonogram.py");
    if (result != 0) {
        std::cerr << "Could not run Python Sonogram code";
    } else {
        std::cout << "Sonogram matrix exported to sonogram_data.csv successfully.\n";
    }
    return 0;


}
