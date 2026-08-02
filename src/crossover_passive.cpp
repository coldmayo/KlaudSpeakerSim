#include "include/crossover_passive.h"

ABCD series_element(std::complex<double> Z) {
    return {1.0, Z, 0.0, 1.0};
}

ABCD shunt_element(std::complex<double> Z) {
    return {1.0, 0.0, 1.0 / Z, 1.0};
}

ABCD multiply(const ABCD & x, const ABCD & y) {
    return {
		x.A * y.A + x.B * y.C,  x.A * y.B + x.B * y.D,
        x.C * y.A + x.D * y.C,  x.C * y.B + x.D * y.D
    };
}

std::complex<double> element_imp(const json& elem, std::complex<double> s) {
    std::string type = elem["type"];
    if (type == "series_inductor") {
        double L = elem["L_h"];
        double dcr = elem.value("dcr_ohm", 0.0);
        return dcr + s * L;
    } else if (type == "shunt_capacitor" || type == "series_capacitor") {
        double C = elem["C_f"];
        return 1.0/(s*C);
    } else if (type == "shunt_inductor") {
        double L = elem["L_h"];
        double dcr = elem.value("dcr_ohm", 0.0);
        return dcr + s * L;
    } else if (type == "shunt_resistor" || type == "series_resistor") {
        return std::complex<double>(elem["R_ohm"], 0.0);
    } else if (type == "zobel_network") {
        double R = elem["R_ohm"];
        double C = elem["C_f"];
        return R + (1.0 / (s * C));
    }
    
    return std::complex<double>(1e-9, 0.0);
}

std::pair<std::vector<std::complex<double>>, std::vector<std::complex<double>>> pass_cross_abcd(const std::vector<std::complex<double>>& s, const std::vector<std::complex<double>>& Z_driver_zma, const json& branch_network) {
    std::vector<std::complex<double>> H_branch(s.size());
    std::vector<std::complex<double>> Z_in_branch(s.size());

    for (int i = 0; i < s.size(); ++i) {
        std::complex<double> Z_L = Z_driver_zma[i];
        ABCD M_total = {1.0, 0.0, 0.0, 1.0};   // identity matrix
        for (const auto& elem : branch_network) {
            std::complex<double> Z = element_imp(elem, s[i]);
            std::string type = elem["type"];
            ABCD M_elem = {1.0, 0.0, 0.0, 1.0};
            if (type.rfind("series", 0) == 0 || type.substr(0, 6) == "series") {
                M_elem = series_element(Z);
            } else if (type.rfind("shunt", 0) == 0 || type.substr(0, 5) == "shunt" || type == "zobel_network") {
                M_elem = shunt_element(Z);
            }

            M_total = multiply(M_total, M_elem);
        }

        H_branch[i] = Z_L / (M_total.A * Z_L + M_total.B);
        Z_in_branch[i] = (M_total.A * Z_L + M_total.B) / (M_total.C * Z_L + M_total.D);
    }
    return {H_branch, Z_in_branch};
}
