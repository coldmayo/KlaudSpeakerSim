#include <iostream>
#include <fstream>
#include <vector>
#include <complex>
#include <cmath>
#include <format>
#include "json.hpp"

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

    std::vector<std::complex<double>> H_driver_woofer(s.size(), 0.0);
    std::vector<std::complex<double>> H_driver_tweeter(s.size(), 0.0);

    std::vector<double> woofer_a, tweeter_a;
    std::vector<double> woofer_offset_m, tweeter_offset_m;

    std::vector<std::complex<double>> Z_woofer_zma(s.size(), std::complex<double>(8.0, 0.0));
    std::vector<std::complex<double>> Z_tweeter_zma(s.size(), std::complex<double>(8.0, 0.0));

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

		source_coords s_c{
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
            std::complex<double> H_delayed_driver = H_raw[j] * H_delay[j]* H_baffle[j];
            if (type == "woofer") {
                H_driver_woofer[j] += H_delayed_driver;
            } else if (type == "tweeter") {
                H_driver_tweeter[j] += H_delayed_driver;
            }
        }

        if (type == "woofer") {
            woofer_a.push_back(radius_m);
            woofer_offset_m.push_back(y_offset);
        } else if (type == "tweeter") {
            tweeter_a.push_back(radius_m);
            tweeter_offset_m.push_back(y_offset);
        }

        if (!zma_file.empty()) {
            auto Z_interp = zma_measured(zma_file, freqs);
            if (type == "woofer") {
                Z_woofer_zma = Z_interp;
            } else if (type == "tweeter") {
                Z_tweeter_zma = Z_interp;
            }
        }

        // Find T/S Parameters
        if (d_json.contains("fund_TS_params") && d_json["t_type"] != "ribbon") {
            zma_data z_ = zma_interp(parsing_zma(zma_file), 250);
            double M_ms = d_json["fund_TS_params"]["M_ms"]; 
            double V_as = d_json["fund_TS_params"]["V_as"];
            double S_d = d_json["eff_surface_area"];
        	double Cms = C_ms(V_as, S_d);
        	double fs = f_s(M_ms, Cms);
        	double Qts = Q_ts(z_.Q_es, z_.Q_ms);
        	double Bl = B_l(fs, M_ms, z_.R_e, z_.Q_es);
            double n0 = n_0(fs, V_as, z_.Q_es);
            double SPL_1W = SPL1W(n0);
            double SPL_283 = SPL283(SPL_1W, z_.R_e);
            double XL = X_L(z_.Z_high, z_.R_e);
            double Le = L_e(XL, z_.f_high);
            double EBP_ = EBP(fs, z_.Q_es);
        	std::cout << "T/S Parameters of the" << d_json["name"] << " Driver:\nC_ms: " << Cms << "\nf_s: " << fs << "\nQts: " << Qts << "\nBl: " << Bl << "\nQ_es: " << z_.Q_es << "\nQ_m: " << z_.Q_ms << "\nn_0: " << n0 << "\nSPL_1W: " << SPL_1W << "\nSPL_2.83W: " << SPL_283 << "\nX_L: " << XL << "\nL_e: " << Le << "\nEBP: " << EBP_ << "\n\n" << std::endl;
        } else if (d_json["t_type"] == "ribbon") {
            zma_data z_ = zma_interp(parsing_zma(zma_file), 20000);
            double M_ms = d_json["fund_TS_params"]["M_ms"];
            double S_d = d_json["eff_surface_area"];
            double XL = X_L(z_.Z_high, z_.R_e);
            double Le = L_e(XL, z_.f_high);
            std::cout << "T/S Parameters of the" << d_json["name"] << " Driver:" << "\nX_L: " << XL << "\nL_e: " << Le << "\n\n" << std::endl;
        }
    }

    double fs = json_f["fs_hz"];
    double fc = json_f["crossover"]["crossover_freq_hz"];
    double omega_0 = 2.0 * PI * fs;
    double omega_c = 2.0 * PI * fc;
    std::vector<std::complex<double>> H_box(s.size(), 1.0);

    // Get transfer functions describing the enclosure of a loudspeaker
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
    std::vector<std::complex<double>> H_HP;
    std::vector<std::complex<double>> H_LP;
    std::vector<std::complex<double>> Z_system(freqs.size());
    std::vector<std::complex<double>> Z_woofer_in(freqs.size());
    std::vector<std::complex<double>> Z_tweeter_in(freqs.size());

	// finding the transfer function for active or passive crossover circuits
	if (json_f["crossover"]["mode"] == "active") {
    	H_HP = select_crossover(json_f["crossover"]["type"], json_f["crossover"]["order"], s, omega_c, 0);
    	H_LP = select_crossover(json_f["crossover"]["type"], json_f["crossover"]["order"], s, omega_c, 1);
    	Z_woofer_in = Z_woofer_zma;
    	Z_tweeter_in = Z_tweeter_zma;
	} else if (json_f["crossover"]["mode"] == "passive") {
        if (json_f["crossover"].contains("tweeter_branch")) {
            std::tie(H_HP, Z_tweeter_in) = pass_cross_abcd(s, Z_tweeter_zma, json_f["crossover"]["tweeter_branch"]);
        }
        if (json_f["crossover"].contains("woofer_branch")) {
            std::tie(H_LP, Z_woofer_in) = pass_cross_abcd(s, Z_woofer_zma, json_f["crossover"]["woofer_branch"]);
        }
	}

	// Finding the systems impulse response
	for (size_t i = 0; i < s.size(); ++i) {
        Z_system[i] = (Z_woofer_in[i] * Z_tweeter_in[i]) / (Z_woofer_in[i] + Z_tweeter_in[i]);
    }

	// finding the total system frequency response
    std::vector<std::complex<double>> H_system_total(freqs.size());
    std::vector<std::complex<double>> H_woofer(freqs.size());
    std::vector<std::complex<double>> H_tweeter(freqs.size());
    for (size_t i = 0; i < freqs.size(); ++i) {
        std::complex<double> woofer_branch  = H_driver_woofer[i] * H_box[i] * H_LP[i];
        std::complex<double> tweeter_branch = H_driver_tweeter[i] * H_HP[i];

        H_woofer[i] = woofer_branch;
        H_tweeter[i] = tweeter_branch;
        H_system_total[i] = woofer_branch + tweeter_branch;
    }

    // Plotting
    std::string box_type = json_f["box_type"];
    std::string cross_mode = json_f["crossover"]["mode"];
	std::string title = std::format("{} {} {}-order {}-drive System", box_type, cross_mode, json_f["crossover"]["order"], json_f["drivers"].size());
	
    plot_t_funcs(H_system_total, H_driver_woofer, H_driver_tweeter, H_woofer, H_tweeter, s, title);
    plot_imp(Z_system, s, "System Impedence");
    plot_filter(H_HP, H_LP, s, "Filter Response");

	// Power
	std::vector<double> P_woofer(s.size()), P_tweeter(s.size());
	double V_in = json_f.value("V_in_tweeter", 2.83);
	for (int i = 0; i < s.size(); ++i) {
    	P_woofer[i]  = (V_in * V_in * std::abs(H_LP[i]) * std::abs(H_LP[i])) / std::real(Z_woofer_in[i]);
        P_tweeter[i] = (V_in * V_in * std::abs(H_HP[i]) * std::abs(H_HP[i])) / std::real(Z_tweeter_in[i]);
	}
	plot_power(P_woofer, P_tweeter, s, "Power Dissipation");

    // Directivity
    std::vector<double> angles_deg;
    for (double deg = -90.0; deg <= 90.0; deg += 2.0) {
        angles_deg.push_back(deg);
    }

    std::vector<std::vector<std::complex<double>>> sonogram_matrix;
    for (double deg : angles_deg) {
        double rad = deg * M_PI / 180.0;

        auto H_theta = polar_response(s, freqs, H_woofer, H_tweeter, woofer_a, tweeter_a, woofer_offset_m, tweeter_offset_m, rad);
        
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
