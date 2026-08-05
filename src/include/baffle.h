#ifndef baffle_h_INCLUDED
#define baffle_h_INCLUDED

#include <vector>
#include <cmath>
#include <iostream>
#include <complex>
#include <cassert>
#include <cstdint>

struct source_coords {
    double x;
    double y;
};

struct listener_coords {
    double x;
    double y;
    double z;
};

struct EdgeSegment {
    double x, y, dl;
};

std::vector<std::complex<double>> b_response(double b_w, double b_h, double driver_r, double dl_step, double fs, int num_samples, source_coords& s_c, listener_coords& l_c);

#endif // baffle_h_INCLUDED
