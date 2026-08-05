#include "include/baffle.h"
// Bit-reversal

constexpr double PI = 3.14159265358979323846;

uint32_t reverse_bits(uint32_t x, int log2N) {
    uint32_t n = 0;
    for (int i = 0; i < log2N; ++i) {
        n = (n << 1) | (x & 1);
        x >>= 1;
    }
    return n;
}

void FFT(std::vector<std::complex<double>>& x) {
    int N = static_cast<int>(x.size());
    assert((N & (N - 1)) == 0 && "FFT size N must be a power of 2");

    int log2N = static_cast<int>(std::log2(N));

    // 1. Bit-reversal permutation
    for (int i = 0; i < N; ++i) {
        int j = static_cast<int>(reverse_bits(static_cast<uint32_t>(i), log2N));
        if (j > i) {
            std::swap(x[i], x[j]);
        }
    }

    for (int size = 2; size <= N; size *= 2) {
        int half = size / 2;
        std::complex<double> w_m = std::polar(1.0, -2.0 * PI / size);

        for (int i = 0; i < N; i += size) {
            std::complex<double> w(1.0, 0.0);

            for (int j = 0; j < half; ++j) {
                std::complex<double> u = x[i + j];
                std::complex<double> t = w * x[i + j + half];

                x[i + j] = u + t;
                x[i + j + half] = u - t;

                w *= w_m;
            }
        }
    }
}
/*
Edge Sources: When the sound wave traveling across the cabinet surface reaches an edge where the geometry suddenly changes (e.g., a $90^\circ$ corner dropping off into free space), the edge acts as a secondary line source radiating a scattered diffracted wave.
*/

double vanderkooy_F_theta(double theta, double wedge_angle = PI / 2.0) {
    double nu = (2.0 * PI - wedge_angle) / PI;
    double numerator = (2.0 / nu) * std::sin(PI / nu);
    double denominator = std::cos(PI / nu) - std::cos(theta / nu);
    constexpr double kMinDenom = 1e-3;
    if (std::abs(denominator) < kMinDenom) {
        denominator = (denominator < 0.0 ? -1.0 : 1.0) * kMinDenom;
    }
    return numerator / denominator;
}

std::vector<std::complex<double>> b_response(double b_w, double b_h, double driver_r, double dl_step, double fs, int num_samples, source_coords& s_c, listener_coords& l_c) {
    double c = 343;
    // Define Baffle Perimeter
    std::vector<std::vector<double>> corners = {{0, 0}, {b_w, 0}, {b_w, b_h}, {0, b_h}, {0,0}};

	// Discretize Edges into Line Segments
    std::vector<EdgeSegment> segments;
    for (int i = 0; i < 4; ++i) {
        std::vector<double> p1 = corners[i];
        std::vector<double> p2 = corners[i + 1];
        double length = std::hypot(p2[0] - p1[0], p2[1] - p1[1]);
        int num_segs = std::max(2, static_cast<int>(std::ceil(length / dl_step)));
        double segment_dl = length / num_segs;
        for (int j = 0; j < num_segs; ++j) {
            double t = (j + 0.5) / num_segs; // Segment midpoint
            segments.push_back({p1[0] + t * (p2[0] - p1[0]), p1[1] + t * (p2[1] - p1[1]), segment_dl});
        }
    }

	// IR buffer
	assert((num_samples & (num_samples - 1)) == 0 && "num_samples must be a power of 2");
    std::vector<double> ir(num_samples, 0.0);

    double Sx = s_c.x, Sy = s_c.y, Sz = 0.0;
    double Lx = l_c.x, Ly = l_c.y, Lz = l_c.z;

    // Direct Sound Path
    double R_direct = std::sqrt((Lx - Sx)*(Lx - Sx) + (Ly - Sy)*(Ly - Sy) + (Lz - Sz)*(Lz - Sz));
    double p_direct_amp = 2.0 / R_direct;
    ir[0] += p_direct_amp;

    // Accumulate 1st-order edge diffractions
    for (const auto& seg : segments) {
		double Qx = seg.x, Qy = seg.y, Qz = 0;

		double r_S = std::sqrt((Qx - Sx)*(Qx - Sx) + (Qy - Sy)*(Qy - Sy));   // source to edge
		double r_P = std::sqrt((Lx - Qx)*(Lx - Qx) + (Ly - Qy)*(Ly - Qy) + (Lz - Qz)*(Lz - Qz)); // Edge to Listener
		
		double t_delay = (r_S + r_P - R_direct) / c;
		double sample_idx = t_delay * fs;

		// Solve for observation angle theta
		double xy_dist = std::hypot(Lx - Qx, Ly - Qy);
        double theta = std::atan2(Lz, xy_dist) + (PI / 2.0);

        double F_th = vanderkooy_F_theta(theta);

        // Impulse strength increment dp
        double dp = (F_th / (2.0 * PI)) * (seg.dl / (r_S * r_P));

        size_t idx0 = static_cast<size_t>(std::floor(sample_idx));
        size_t idx1 = idx0 + 1;
        double frac = sample_idx - idx0;
 
        if (idx0 < static_cast<size_t>(num_samples)) ir[idx0] += dp * (1.0 - frac);
        if (idx1 < static_cast<size_t>(num_samples)) ir[idx1] += dp * frac;
        
    }

    std::vector<std::complex<double>> H_raw(ir.begin(), ir.end());
    FFT(H_raw);

    return H_raw;
}
