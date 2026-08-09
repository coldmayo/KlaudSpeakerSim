#ifndef Thiele_Small_h_INCLUDED
#define Thiele_Small_h_INCLUDED

#include "include/parsing.h"

#include <vector>
#include <cmath>
#include <iostream>
#include <complex>
#include <algorithm>

struct zma_data {
    double R_e;
    double F_s;
    double Zmax;
    double Q_es;
    double Q_ms;
    double f_high;
    double Z_high;
};

zma_data zma_interp(const zma& data, double max_search_freq);
double f_s(double M_ms, double C_ms);
double Q_es(double f_s, double M_ms, double C_ms, double R_e);
double Q_ms(double f_s, double M_ms, double R_ms);
double Q_ts(double Q_es, double Q_ms);
double V_as(double S_d, double C_ms);
double B_l(double f_s, double M_ms, double R_e, double Q_es);
double R_ms(double S_d, double Q_ms, double V_as, double f_s);
double C_ms(double V_as, double S_d);
double n_0(double f_s, double V_as, double Q_es);
double SPL1W(double n0);
double SPL283(double SPL1W, double R_e);
double X_L(double Z, double R_e);
double L_e(double X_L, double f_s);
double EBP(double f_s, double Q_es);

#endif // Thiele_Small_h_INCLUDED
