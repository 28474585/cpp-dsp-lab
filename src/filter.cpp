#include "dsp/filter.h"

#include <cmath>
#include <stdexcept>
#include <string>

namespace dsp {

namespace {

constexpr double kPi = 3.14159265358979323846;

// 归一化 sinc:sin(πx)/(πx)，x = 0 时按极限取 1。
double sinc(double x) {
    if (std::abs(x) < 1e-12) return 1.0;
    const double px = kPi * x;
    return std::sin(px) / px;
}

}  // namespace

// ---------------------------------------------------------------- FIR

std::vector<double> designFirLowpass(double cutoff, std::size_t numTaps, WindowType win) {
    if (cutoff <= 0.0 || cutoff >= 0.5) {
        throw std::invalid_argument("cutoff 必须在 (0, 0.5) 内，当前为 " + std::to_string(cutoff));
    }
    if (numTaps == 0) throw std::invalid_argument("numTaps 必须大于 0");

    const std::vector<double> w = window(win, numTaps);
    // 线性相位 FIR 的冲激响应关于中心对称，中心在 (N-1)/2。奇数 N 使中心落在采样点上。
    const double center = (static_cast<double>(numTaps) - 1.0) / 2.0;

    std::vector<double> h(numTaps);
    for (std::size_t i = 0; i < numTaps; ++i) {
        // 理想低通的冲激响应是一条无限长的 sinc；乘窗即把它截断成有限长。
        const double n = static_cast<double>(i) - center;
        h[i] = 2.0 * cutoff * sinc(2.0 * cutoff * n) * w[i];
    }

    // 归一化到直流增益 1，保证通带增益为 1。
    double sum = 0.0;
    for (double v : h) sum += v;
    if (sum != 0.0) {
        for (double& v : h) v /= sum;
    }
    return h;
}

Signal firFilter(const std::vector<double>& h, const Signal& x) {
    if (h.empty() || x.empty()) return {};

    Signal y(x.size() + h.size() - 1, Complex(0.0, 0.0));
    for (std::size_t n = 0; n < x.size(); ++n) {
        // 上采样后的信号里有大量精确的 0，跳过可直接省掉 (sps-1)/sps 的乘法。
        if (x[n] == Complex(0.0, 0.0)) continue;
        for (std::size_t k = 0; k < h.size(); ++k) {
            y[n + k] += x[n] * h[k];
        }
    }
    return y;
}

std::vector<Complex> freqResponse(const std::vector<double>& h, std::size_t points) {
    std::vector<Complex> out(points);
    for (std::size_t p = 0; p < points; ++p) {
        const double w = (points == 1)
                             ? 0.0
                             : kPi * static_cast<double>(p) / static_cast<double>(points - 1);
        Complex acc(0.0, 0.0);
        for (std::size_t k = 0; k < h.size(); ++k) {
            acc += h[k] * std::exp(Complex(0.0, -w * static_cast<double>(k)));
        }
        out[p] = acc;
    }
    return out;
}

double groupDelay(std::size_t numTaps) { return (static_cast<double>(numTaps) - 1.0) / 2.0; }

// ---------------------------------------------------------------- IIR

Biquad biquadLowpass(double cutoff, double q) {
    if (cutoff <= 0.0 || cutoff >= 0.5) {
        throw std::invalid_argument("cutoff 必须在 (0, 0.5) 内，当前为 " + std::to_string(cutoff));
    }
    if (q <= 0.0) throw std::invalid_argument("q 必须大于 0");

    const double w0 = 2.0 * kPi * cutoff;
    const double cw = std::cos(w0);
    const double alpha = std::sin(w0) / (2.0 * q);

    const double a0 = 1.0 + alpha;
    Biquad bq;
    bq.b0 = ((1.0 - cw) / 2.0) / a0;
    bq.b1 = (1.0 - cw) / a0;
    bq.b2 = ((1.0 - cw) / 2.0) / a0;
    bq.a1 = (-2.0 * cw) / a0;
    bq.a2 = (1.0 - alpha) / a0;
    return bq;
}

Signal biquadFilter(const Biquad& bq, const Signal& x) {
    Signal y(x.size(), Complex(0.0, 0.0));
    Complex x1(0.0, 0.0), x2(0.0, 0.0), y1(0.0, 0.0), y2(0.0, 0.0);

    for (std::size_t n = 0; n < x.size(); ++n) {
        const Complex xn = x[n];
        const Complex yn = bq.b0 * xn + bq.b1 * x1 + bq.b2 * x2 - bq.a1 * y1 - bq.a2 * y2;
        y[n] = yn;
        x2 = x1;
        x1 = xn;
        y2 = y1;
        y1 = yn;
    }
    return y;
}

std::vector<Complex> biquadResponse(const Biquad& bq, std::size_t points) {
    std::vector<Complex> out(points);
    for (std::size_t p = 0; p < points; ++p) {
        const double w = (points == 1)
                             ? 0.0
                             : kPi * static_cast<double>(p) / static_cast<double>(points - 1);
        const Complex z1 = std::exp(Complex(0.0, -w));
        const Complex z2 = z1 * z1;
        const Complex num = bq.b0 + bq.b1 * z1 + bq.b2 * z2;
        const Complex den = Complex(1.0, 0.0) + bq.a1 * z1 + bq.a2 * z2;
        out[p] = num / den;
    }
    return out;
}

}  // namespace dsp
