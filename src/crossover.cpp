#include "include/crossover.h"

/*
	First order:
    	Single component per driver. 
    	Tweeters use a capacitor for a high pass filer
    	Subs use a single inductor (coil) for a low pass filter
    	6dB per octave slope
*/

std::vector<std::complex<double>> first_order_HP(std::vector<std::complex<double>> s, double omega_c) {
    std::vector<std::complex<double>> H_HP(s.size());
    for (int i = 0; i < s.size(); ++i) {
        H_HP[i] = s[i] / (s[i] + omega_c);
    }
    return H_HP;
}

std::vector<std::complex<double>> first_order_LP(std::vector<std::complex<double>> s, double omega_c) {
    std::vector<std::complex<double>> H_HP(s.size());
    for (int i = 0; i < s.size(); ++i) {
        H_HP[i] = omega_c / (s[i] + omega_c);
    }
    return H_HP;
}

/*
	Second order:
    	Two components per driver
    	Butterworth: Q = 0.707; Linkwitz-Riley 2nd order: Q = 0.5; Bessel: Q = 0.577
    	12dB per octave slope
*/

std::vector<std::complex<double>> second_order_HP(std::vector<std::complex<double>> s, double omega_c, double Q) {
    std::vector<std::complex<double>> H(s.size());
    const double omega2_c = omega_c*omega_c;
    const double damping = omega_c/Q;
    for (int i = 0; i < s.size(); i++) {
        std::complex<double> s2 = s[i]*s[i];
        H[i] = s2/(s2 + (damping)*s[i] + omega2_c);
    }
    return H;
}

std::vector<std::complex<double>> second_order_LP(std::vector<std::complex<double>> s, double omega_c, double Q) {
    std::vector<std::complex<double>> H(s.size());
    const double omega2_c = omega_c*omega_c;
    const double damping = omega_c/Q;
    for (int i = 0; i < s.size(); i++) {
        std::complex<double> s2 = s[i]*s[i];
        H[i] = omega2_c/(s2 + (damping)*s[i] + omega2_c);
    }
    return H;
}

/*
	Third Order:
    	Not super common since you could use the 4th order version
    	18dB per octave slope
*/

std::vector<std::complex<double>> third_order_LP(std::vector<std::complex<double>> s, double omega_c, double Q) {
    std::vector<std::complex<double>> H(s.size());
    const double omega2_c = omega_c*omega_c;
    const double omega3_c = omega_c*omega_c*omega_c;
    for (int i = 0; i < s.size(); ++i) {
        std::complex<double> s2 = s[i]*s[i];
        std::complex<double> s3 = s[i]*s[i]*s[i];
        H[i] = omega3_c/(s3 + 2*omega_c*s2 + 2*omega2_c*s[i] + omega3_c);
    }
    return H;
}

std::vector<std::complex<double>> third_order_HP(std::vector<std::complex<double>> s, double omega_c, double Q) {
    std::vector<std::complex<double>> H(s.size());
    const double omega2_c = omega_c*omega_c;
    const double omega3_c = omega_c*omega_c*omega_c;
    for (int i = 0; i < s.size(); ++i) {
        std::complex<double> s2 = s[i]*s[i];
        std::complex<double> s3 = s[i]*s[i]*s[i];
        H[i] = s3/(s3 + 2*omega_c*s2 + 2*omega2_c*s[i] + omega3_c);
    }
    return H;
}

/*
	Fourth Order:
    	Only the Linkwitz-Riley
    	Default for industry
*/

std::vector<std::complex<double>> fourth_order_HP(std::vector<std::complex<double>> s, double omega_c, double Q) {
    std::vector<std::complex<double>> H(s.size());
    auto H_2 = second_order_HP(s, omega_c, Q);
    for (int i = 0; i < s.size(); ++i) {
        H[i] = H_2[i]*H_2[i];
    }
    return H;
}

std::vector<std::complex<double>> fourth_order_LP(std::vector<std::complex<double>> s, double omega_c, double Q) {
    std::vector<std::complex<double>> H(s.size());
    auto H_2 = second_order_LP(s, omega_c, Q);
    for (int i = 0; i < s.size(); ++i) {
        H[i] = H_2[i]*H_2[i];
    }
    return H;
}

// HPLP_toggle, 0: HP, 1: LP
std::vector<std::complex<double>> select_crossover(std::string type, int order, std::vector<std::complex<double>> s, double omega_c, int HPLP_toggle) {
	double Q = 0;
	if (type == "linkwitz_riley") {
    	Q = 0.5;
	} else if (type == "butterworth") {
    	Q = 0.707;
	} else if (type == "bessel") {
    	Q = 0.577;
	}

	switch(order) {
    	case 1:
        	if (HPLP_toggle == 0) {
            	return first_order_HP(s, omega_c);
        	} else {
            	return first_order_LP(s, omega_c);
        	}
        case 2:
            if (HPLP_toggle == 0) {
            	return second_order_HP(s, omega_c, Q);
        	} else {
            	return second_order_LP(s, omega_c, Q);
        	}
        case 3:
            if (HPLP_toggle == 0) {
            	return third_order_HP(s, omega_c, Q);
        	} else {
            	return third_order_LP(s, omega_c, Q);
        	}
        case 4:
            if (HPLP_toggle == 0) {
            	return fourth_order_HP(s, omega_c, Q);
        	} else {
            	return fourth_order_LP(s, omega_c, Q);
        	}
        default:
            return first_order_HP(s, omega_c);
	}
}
