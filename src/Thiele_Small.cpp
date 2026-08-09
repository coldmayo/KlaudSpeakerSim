#include "include/Thiele_Small.h"
#include "include/parsing.h"

// https://en.wikipedia.org/wiki/Thiele/Small_parameters#Small_signal_parameters

constexpr double PI = 3.14159265358979323846;

zma_data zma_interp(const zma& data, double max_search_freq) {
    zma_data d = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    double Z_max = -1.0;
    int peak_index = 0;

    for (size_t i = 0; i < data.freq.size(); ++i) {
        if (data.freq[i] > max_search_freq) break;

        if (data.imp[i] > Z_max) {
            Z_max = data.imp[i];
            d.F_s = data.freq[i];
            d.Zmax = Z_max;
            peak_index = static_cast<int>(i);
        }
    }

    double f_high = 10000.0;
    double Z_at_high_freq = Z_max;

    for (size_t i = 0; i < data.freq.size(); ++i) {
        if (data.freq[i] >= f_high) {
            f_high = data.freq[i];
            Z_at_high_freq = data.imp[i];
            break;
        }
    }
    d.f_high = f_high;
    d.Z_high = Z_at_high_freq;

    double imp_min = 1000000000.0;
    for (int j = peak_index; j >= 0; --j) {
        if (data.imp[j] < imp_min) {
            imp_min = data.imp[j];
        }
    }
    d.R_e = imp_min;

    double Z_c = std::sqrt(d.Zmax * d.R_e);
    double f1 = 0.0;
    double f2 = 0.0;

    for (int i = 0; i < peak_index; ++i) {
        if (data.imp[i] >= Z_c) {
            f1 = data.freq[i]; // Simple approximation (or linear interpolate between i-1 and i)
            break;
        }
    }

    // Find f2 (above Fs) by sweeping past the peak index
    for (size_t i = peak_index; i < data.freq.size(); ++i) {
        if (data.imp[i] <= Z_c) {
            f2 = data.freq[i];
            break;
        }
    }

    if (f2 - f1 > 0.001) {
        d.Q_ms = (d.F_s * std::sqrt(d.Zmax / d.R_e)) / (f2 - f1);
        d.Q_es = d.Q_ms / ((d.Zmax / d.R_e) - 1.0);
    } else {
        d.Q_ms = 0.0;
        d.Q_es = 0.0;
    }

    return d;
}

double f_s(double M_ms, double C_ms) {
    return 1/(2.0*PI*std::sqrt(C_ms*M_ms));
}

double Q_es(double f_s, double M_ms, double C_ms, double R_e, double B_l) {
    return (2*PI*f_s*M_ms*R_e)/(std::pow(B_l, 2));
}

double Q_ms(double f_s, double M_ms, double R_ms) {
    return (2*PI*f_s*M_ms)/R_ms;
}

double Q_ts(double Q_es, double Q_ms) {
    return (Q_ms*Q_es)/(Q_ms+Q_es);
}

double V_as(double S_d, double C_ms) {
    double rho = 1.184;
    double c = 346.1;
    return rho * (std::pow(c, 2)) * (std::pow(S_d, 2)) * C_ms;
}

double B_l(double f_s, double M_ms, double R_e, double Q_es) {
    return std::sqrt((2*PI*f_s*M_ms*R_e)/Q_es);
}

double R_ms(double S_d, double Q_ms, double V_as, double f_s) {
    double rho = 1.184;
    double c = 346.1;
    return (rho*(std::pow(c, 2))*(std::pow(S_d, 2)))/(2*PI*f_s*Q_ms*V_as);
}

double C_ms(double V_as, double S_d) {
    double rho = 1.184;
    double c = 346.1;
    return (V_as)/((rho*(std::pow(c, 2)))*(std::pow(S_d, 2)));
}

double n_0(double f_s, double V_as, double Q_es) {
    double c = 346.1;
    return (4*std::pow(PI, 2)*std::pow(f_s, 3)*V_as)/(std::pow(c, 3)*Q_es);
}

double SPL1W(double n0) {
    return 112.2 + 10*std::log10(n0);
}

double SPL283(double SPL1W, double R_e) {
    return SPL1W + 20*std::log10((std::sqrt(8.0))/(std::sqrt(R_e)));
}

double X_L(double Z, double R_e) {
    return std::sqrt(std::pow(Z, 2)-std::pow(R_e, 2));
}

double L_e(double X_L, double f_s) {
    return X_L/(2*PI*f_s);
}

double EBP(double f_s, double Q_es) {
    return f_s/Q_es;
}
