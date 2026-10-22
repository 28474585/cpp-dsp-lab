#include "dsp/window.h"

#include <cmath>

namespace dsp {

namespace {

constexpr double kTwoPi = 6.283185307179586476925286766559;

}  // namespace

std::vector<double> window(WindowType type, std::size_t n) {
    if (n == 0) return {};
    if (n == 1) return {1.0};

    std::vector<double> w(n);
    // 分母取 N-1：让首尾两端正好落在余弦的 0 与 2π 处，窗序列关于中心对称。
    const double denom = static_cast<double>(n - 1);

    for (std::size_t i = 0; i < n; ++i) {
        const double x = kTwoPi * static_cast<double>(i) / denom;
        switch (type) {
            case WindowType::Rectangular:
                w[i] = 1.0;
                break;
            case WindowType::Hann:
                w[i] = 0.5 - 0.5 * std::cos(x);
                break;
            case WindowType::Hamming:
                w[i] = 0.54 - 0.46 * std::cos(x);
                break;
            case WindowType::Blackman:
                w[i] = 0.42 - 0.5 * std::cos(x) + 0.08 * std::cos(2.0 * x);
                break;
        }
    }
    return w;
}

const char* windowName(WindowType type) {
    switch (type) {
        case WindowType::Rectangular: return "矩形窗";
        case WindowType::Hann:        return "汉宁窗";
        case WindowType::Hamming:     return "汉明窗";
        case WindowType::Blackman:    return "布莱克曼窗";
    }
    return "未知窗";
}

}  // namespace dsp
