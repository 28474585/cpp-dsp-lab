// 阶段二演示：FIR / IIR 滤波器设计 + 窗函数 + 频谱泄漏。
//
// 输出：
//   - 终端：窗函数旁瓣对比表、滤波器指标（直流增益、群延迟、阻带衰减）
//   - output/filter_freq_response.csv：三种窗设计的 FIR 与一个 IIR 的幅频响应(dB)
//   - output/filter_time_domain.csv  ：同一信号滤波前后的时域波形
//   - output/filter_impulse_response.csv：IIR 的冲激响应（看衰减与稳定性）
//   - output/spectral_leakage.csv     ：非整数周期正弦用矩形窗 / 汉明窗的频谱对比

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "dsp/console.h"
#include "dsp/csv.h"
#include "dsp/fft.h"
#include "dsp/filter.h"
#include "dsp/window.h"

using namespace dsp;

namespace {

constexpr double kTwoPi = 6.283185307179586476925286766559;

// 复数的幅度换算成 dB（并夹一个下限，避免 log(0)）
double toDb(const Complex& v) {
    const double m = std::abs(v);
    return 20.0 * std::log10(std::max(m, 1e-12));
}

// 第一旁瓣电平（dB，相对主瓣峰值）。
// 注意：加窗 sinc 滤波器在通带里还有很小的"纹波"（Gibbs 现象），深度不到 1 dB，
// 直接用"第一个局部极小"会误把这种纹波当成主瓣后的第一零点。因此要求极小值必须
// 比峰值低 3 dB 以上，才是真正的主瓣-旁瓣分界零点。
double firstSidelobeDb(const std::vector<double>& h, std::size_t points) {
    const std::vector<Complex> H = freqResponse(h, points);
    std::vector<double> db(points);
    for (std::size_t i = 0; i < points; ++i) db[i] = toDb(H[i]);

    const double peak = *std::max_element(db.begin(), db.end());

    std::size_t nullIdx = points - 1;
    for (std::size_t i = 1; i + 1 < points; ++i) {
        const bool localMin = (db[i] <= db[i - 1]) && (db[i] <= db[i + 1]);
        if (localMin && db[i] < peak - 3.0) {
            nullIdx = i;
            break;
        }
    }

    // 从第一零点向右爬到第一个局部极大，即第一旁瓣峰。
    std::size_t j = nullIdx;
    while (j + 1 < points && db[j + 1] >= db[j]) ++j;
    return db[j] - peak;
}

void compareWindows() {
    std::cout << "===== 窗函数对比（低通 FIR，截止 0.1，65 抽头）=====\n";
    std::cout << "窗口            第一旁瓣(dB)   直流增益\n";

    const double cutoff = 0.1;
    const std::size_t taps = 65;

    std::vector<CsvRow> rows;
    const std::vector<std::pair<WindowType, const char*>> tried = {
        {WindowType::Rectangular, "rect_db"},
        {WindowType::Hann, "hann_db"},
        {WindowType::Hamming, "hamming_db"},
        {WindowType::Blackman, "blackman_db"},
    };

    // 频响点一样，横向拼成一张表，方便绘图直接对比四条曲线。
    const std::size_t points = 512;
    std::vector<std::vector<double>> curves(tried.size());

    for (std::size_t t = 0; t < tried.size(); ++t) {
        const auto h = designFirLowpass(cutoff, taps, tried[t].first);
        double dc = 0.0;
        for (double v : h) dc += v;
        const double side = firstSidelobeDb(h, 4096);

        std::cout << std::left << std::setw(14) << windowName(tried[t].first) << std::right
                  << std::setw(12) << fmt(side, 1) << std::setw(12) << fmt(dc, 4) << "\n";

        const auto H = freqResponse(h, points);
        curves[t].resize(points);
        for (std::size_t i = 0; i < points; ++i) curves[t][i] = toDb(H[i]);
    }
    std::cout << "结论：旁瓣越低 → 阻带衰减越好（泄漏/串扰越小）；代价是主瓣更宽、过渡带更缓。\n\n";

    // 再算一条 IIR 曲线放同一张图，直观对比"低阶 IIR 也能给出陡峭的过渡带"。
    const Biquad bq = biquadLowpass(cutoff, 0.70710678118654752);
    const auto Hi = biquadResponse(bq, points);
    std::vector<double> iirDb(points);
    for (std::size_t i = 0; i < points; ++i) iirDb[i] = toDb(Hi[i]);

    for (std::size_t i = 0; i < points; ++i) {
        const double f = static_cast<double>(i) / static_cast<double>(points - 1) * 0.5;  // 归一化频率
        rows.push_back({fmt(f, 6), fmt(curves[0][i], 4), fmt(curves[1][i], 4),
                        fmt(curves[2][i], 4), fmt(curves[3][i], 4), fmt(iirDb[i], 4)});
    }
    writeCsv("output/filter_freq_response.csv",
             {"freq_norm", "rect_db", "hann_db", "hamming_db", "blackman_db", "iir_db"}, rows);
    std::cout << "幅频响应已写入 output/filter_freq_response.csv\n\n";
}

void filterDemo() {
    // 构造一个"低频有用 + 高频干扰"的信号，用低通滤波器把高频削掉。
    // 两个频率都取成整数个周期（k/n），这样 DFT 投影不泄漏，幅度读数才干净。
    const std::size_t n = 256;
    const std::size_t kLow = 8;      // 8/256 个周期/采样点
    const std::size_t kHigh = 77;    // 77/256，高于截止 0.1，应被滤除
    const double fLow = static_cast<double>(kLow) / static_cast<double>(n);
    const double fHigh = static_cast<double>(kHigh) / static_cast<double>(n);

    Signal x(n);
    for (std::size_t i = 0; i < n; ++i) {
        const double v = std::cos(kTwoPi * fLow * i) + 0.6 * std::cos(kTwoPi * fHigh * i);
        x[i] = Complex(v, 0.0);
    }

    const auto h = designFirLowpass(0.1, 65, WindowType::Hamming);
    const Signal y = firFilter(h, x);
    const std::size_t delay = static_cast<std::size_t>(groupDelay(h.size()));

    std::vector<CsvRow> rows;
    for (std::size_t i = 0; i < n; ++i) {
        // 输出右移了 delay 个点，这里对齐后输出，方便与输入直接比较。
        const std::size_t j = i + delay;
        const double out = (j < y.size()) ? y[j].real() : 0.0;
        rows.push_back({fmtInt(static_cast<long long>(i)), fmt(x[i].real(), 6), fmt(out, 6)});
    }
    writeCsv("output/filter_time_domain.csv", {"n", "input", "output"}, rows);

    // 量化一下：滤波前后各频率分量的幅度。
    // 对实余弦做单频投影，只取到正频率分量，所以乘 2 还原出真实幅度。
    // 输出要从 delay 处开始取，避开卷积前沿那段还没起来的暂态。
    auto ampAt = [&](const Signal& s, double f, std::size_t len, std::size_t start) {
        Complex acc(0.0, 0.0);
        for (std::size_t i = 0; i < len && start + i < s.size(); ++i) {
            acc += s[start + i] * std::exp(Complex(0.0, -kTwoPi * f * static_cast<double>(i)));
        }
        return 2.0 * std::abs(acc) / static_cast<double>(len);
    };
    std::cout << "===== 低通滤波（截止 0.1，汉明窗 FIR）=====\n";
    std::cout << "滤波前 高频分量 (f=0.30) 幅度: " << fmt(ampAt(x, fHigh, n, 0), 4) << "\n";
    std::cout << "滤波后 高频分量 (f=0.30) 幅度: " << fmt(ampAt(y, fHigh, n, delay), 4) << "\n";
    std::cout << "滤波后 低频分量 (f=0.03) 幅度: " << fmt(ampAt(y, fLow, n, delay), 4)
              << "（应保持 ≈1）\n";
    std::cout << "时域对比已写入 output/filter_time_domain.csv\n\n";

    // IIR 冲激响应：看它是否衰减到 0（稳定性的直观证据）
    const Biquad bq = biquadLowpass(0.1, 0.70710678118654752);
    Signal imp(n, Complex(0.0, 0.0));
    imp[0] = Complex(1.0, 0.0);
    const Signal ir = biquadFilter(bq, imp);

    std::vector<CsvRow> irows;
    for (std::size_t i = 0; i < n; ++i) {
        irows.push_back({fmtInt(static_cast<long long>(i)), fmt(ir[i].real(), 8)});
    }
    writeCsv("output/filter_impulse_response.csv", {"n", "iir_impulse"}, irows);

    double tail = 0.0;
    for (std::size_t i = n / 2; i < n; ++i) tail = std::max(tail, std::abs(ir[i]));
    std::cout << "IIR 冲激响应（截止 0.1，Q=0.7071）后段最大幅度: " << tail
              << "（趋于 0 说明极点都在单位圆内，滤波器稳定）\n";
    std::cout << "冲激响应已写入 output/filter_impulse_response.csv\n\n";
}

void spectralLeakage() {
    // 泄漏的根源：截断时长不是信号周期的整数倍。
    // 这里让频率落在两个频点之间（10.3 个周期），矩形窗截断会在谱上产生很高的旁瓣，
    // 把能量"漏"到邻近频点；换成汉明窗可以明显压低这些旁瓣。
    const std::size_t n = 256;
    const double cycles = 10.3;

    Signal x(n);
    for (std::size_t i = 0; i < n; ++i) {
        x[i] = Complex(std::cos(kTwoPi * cycles * static_cast<double>(i) / n), 0.0);
    }

    const auto wRect = window(WindowType::Rectangular, n);
    const auto wHam = window(WindowType::Hamming, n);

    Signal xRect(n), xHam(n);
    for (std::size_t i = 0; i < n; ++i) {
        xRect[i] = x[i] * wRect[i];
        xHam[i] = x[i] * wHam[i];
    }
    const Signal XRect = fft(xRect);
    const Signal XHam = fft(xHam);

    // 归一化到各自的主瓣峰值，这样比较的是"旁瓣相对主瓣有多高"。
    std::vector<double> ar(n / 2 + 1), ah(n / 2 + 1);
    for (std::size_t k = 0; k <= n / 2; ++k) {
        ar[k] = std::abs(XRect[k]);
        ah[k] = std::abs(XHam[k]);
    }
    const double pr = *std::max_element(ar.begin(), ar.end());
    const double ph = *std::max_element(ah.begin(), ah.end());

    // 主瓣所在频点（用未加窗的谱找），旁瓣统计时跳过主瓣附近的 guard 个频点。
    std::size_t peakBin = 0;
    for (std::size_t k = 1; k <= n / 2; ++k) {
        if (ar[k] > ar[peakBin]) peakBin = k;
    }
    const std::size_t guard = 4;

    std::vector<CsvRow> rows;
    double sidelobeRect = -1e9, sidelobeHam = -1e9;
    for (std::size_t k = 0; k <= n / 2; ++k) {
        const double dr = 20.0 * std::log10(std::max(ar[k] / pr, 1e-12));
        const double dh = 20.0 * std::log10(std::max(ah[k] / ph, 1e-12));
        rows.push_back({fmtInt(static_cast<long long>(k)), fmt(dr, 4), fmt(dh, 4)});

        // 主瓣之外（跳过 guard 个频点）的最高电平，即约等于第一旁瓣。
        if (k > peakBin + guard) {
            sidelobeRect = std::max(sidelobeRect, dr);
            sidelobeHam = std::max(sidelobeHam, dh);
        }
    }
    writeCsv("output/spectral_leakage.csv", {"bin", "rect_db", "hamming_db"}, rows);

    std::cout << "===== 频谱泄漏（10.3 个周期的正弦，非整数截断）=====\n";
    std::cout << "矩形窗  主瓣外最高电平: " << fmt(sidelobeRect, 1) << " dB\n";
    std::cout << "汉明窗  主瓣外最高电平: " << fmt(sidelobeHam, 1) << " dB\n";
    std::cout << "结论：矩形窗旁瓣高、泄漏严重；汉明窗把泄漏压下去，代价是主瓣变宽、分辨率下降。\n";
    std::cout << "泄漏对比已写入 output/spectral_leakage.csv\n\n";
}

}  // namespace

int main() {
    enableUtf8Console();
    try {
        std::cout << "FIR / IIR 滤波器与窗函数演示（阶段二）\n\n";
        compareWindows();
        filterDemo();
        spectralLeakage();
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "错误: " << e.what() << "\n";
        return 1;
    }
}
