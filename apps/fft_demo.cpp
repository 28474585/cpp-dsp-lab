// 阶段一演示：朴素 DFT 与基 2 FFT 的正确性对照 + 性能对比。
//
// 输出：
//   - 终端：正确性验证结论 + 性能对比表
//   - output/fft_bench.csv：各长度的耗时与加速比
//   - output/fft_spectrum.csv：一个示例信号的频谱（便于绘图）

#include <algorithm>
#include <chrono>
#include <cmath>
#include <complex>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "dsp/console.h"
#include "dsp/csv.h"
#include "dsp/fft.h"

using namespace dsp;

namespace {

constexpr double kTwoPi = 6.283185307179586476925286766559;

// 可复现的伪随机数（线性同余），不依赖 rand()，保证每次运行结果一致。
class Lcg {
public:
    explicit Lcg(std::uint64_t seed) : s_(seed) {}
    double next() {
        s_ = s_ * 6364136223846793005ULL + 1442695040888963407ULL;
        // 取高 53 位映射到 [-1, 1)
        const std::uint64_t hi = s_ >> 11;
        return static_cast<double>(hi) / static_cast<double>(1ULL << 53) * 2.0 - 1.0;
    }

private:
    std::uint64_t s_;
};

Signal makeTestSignal(std::size_t n, std::uint64_t seed = 20260920ULL) {
    Lcg rng(seed);
    Signal x(n);
    for (std::size_t i = 0; i < n; ++i) {
        x[i] = Complex(rng.next(), rng.next());
    }
    return x;
}

double maxAbsDiff(const Signal& a, const Signal& b) {
    double m = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i) m = std::max(m, std::abs(a[i] - b[i]));
    return m;
}

double energy(const Signal& x) {
    double e = 0.0;
    for (const auto& v : x) e += std::norm(v);
    return e;
}

template <typename Fn>
double timeMs(Fn&& fn, int reps) {
    const auto t0 = std::chrono::steady_clock::now();
    for (int i = 0; i < reps; ++i) {
        // 防止编译器把整个调用优化掉
        volatile double sink = std::abs(fn()[0]);
        (void)sink;
    }
    const auto t1 = std::chrono::steady_clock::now();
    return std::chrono::duration<double, std::milli>(t1 - t0).count() / reps;
}

void verifyCorrectness() {
    std::cout << "===== 正确性验证 =====\n";

    // 1. FFT vs 朴素 DFT
    {
        const std::size_t n = 1024;
        const Signal x = makeTestSignal(n);
        const double err = maxAbsDiff(dft(x), fft(x));
        std::cout << "[1] FFT 与朴素 DFT 的最大绝对误差 (N=" << n << "): " << err << "\n";
    }
    // 2. 迭代 DIT 与递归 DIF 两条独立路径互相印证
    {
        const std::size_t n = 512;
        const Signal x = makeTestSignal(n);
        const double err = maxAbsDiff(fft(x), fftRecursiveDIF(x));
        std::cout << "[2] 迭代 DIT 与递归 DIF 的最大绝对误差 (N=" << n << "): " << err << "\n";
    }
    // 3. Parseval 定理：时域能量 == 频域能量 / N
    {
        const std::size_t n = 256;
        const Signal x = makeTestSignal(n);
        const Signal X = fft(x);
        const double et = energy(x);
        const double ef = energy(X) / static_cast<double>(n);
        std::cout << "[3] Parseval: 时域能量=" << et << "  频域能量/N=" << ef << "\n";
    }
    // 4. IFFT 往返
    {
        const std::size_t n = 128;
        const Signal x = makeTestSignal(n);
        const double err = maxAbsDiff(x, ifft(fft(x)));
        std::cout << "[4] FFT -> IFFT 往返最大误差 (N=" << n << "): " << err << "\n";
    }
    // 5. 已知闭式频谱：单位冲激 -> 全 1
    {
        const std::size_t n = 64;
        Signal x(n, Complex(0.0, 0.0));
        x[0] = Complex(1.0, 0.0);
        const Signal X = fft(x);
        double err = 0.0;
        for (const auto& v : X) err = std::max(err, std::abs(v - Complex(1.0, 0.0)));
        std::cout << "[5] 单位冲激 -> 全 1 频响，最大误差: " << err << "\n";
    }
    // 6. 已知闭式频谱：直流 -> 只有第 0 点非零，值为 N
    {
        const std::size_t n = 64;
        Signal x(n, Complex(2.0, 0.0));
        const Signal X = fft(x);
        double offCenter = 0.0;
        for (std::size_t k = 1; k < n; ++k) offCenter = std::max(offCenter, std::abs(X[k]));
        std::cout << "[6] 直流信号 -> 第 0 点为 " << X[0].real() << " (期望 " << 2.0 * n
                  << ")，其余频点最大幅度: " << offCenter << "\n";
    }
    // 7. 单频正弦 -> 对应频点出现峰值
    {
        const std::size_t n = 256;
        const std::size_t k0 = 19;  // 目标频点
        Signal x(n);
        for (std::size_t i = 0; i < n; ++i) {
            x[i] = Complex(std::cos(kTwoPi * static_cast<double>(k0) * i / n), 0.0);
        }
        const Signal X = fft(x);
        std::size_t peak = 0;
        for (std::size_t k = 0; k < n; ++k) {
            if (std::abs(X[k]) > std::abs(X[peak])) peak = k;
        }
        std::cout << "[7] 单频信号 (k0=" << k0 << ") 的频谱峰值出现在 k=" << peak
                  << "，峰值幅度 " << std::abs(X[peak]) << " (期望 " << n / 2.0 << ")\n";
    }
    std::cout << "\n";
}

void benchmark() {
    std::cout << "===== 性能对比 =====\n";
    std::cout << "        N      DFT(ms)      FFT(ms)      加速比\n";

    std::vector<CsvRow> rows;

    // 朴素 DFT 是 O(N^2)，到 2^13 已经要几秒，再大就失去意义，故只测到 2^13
    const std::size_t dftMaxPow = 13;
    const std::size_t fftMaxPow = 20;

    for (std::size_t p = 8; p <= fftMaxPow; ++p) {
        const std::size_t n = static_cast<std::size_t>(1) << p;
        const Signal x = makeTestSignal(n, 12345ULL + p);

        double dftMs = -1.0;
        if (p <= dftMaxPow) {
            const int reps = (p <= 11) ? 3 : 1;
            dftMs = timeMs([&] { return dft(x); }, reps);
        }
        const int fftReps = (p <= 16) ? 20 : 5;
        const double fftMs = timeMs([&] { return fft(x); }, fftReps);

        const double speedup = dftMs > 0 ? dftMs / fftMs : 0.0;

        if (dftMs > 0) {
            std::cout << std::setw(9) << n << std::setw(12) << fmt(dftMs, 3) << std::setw(13)
                      << fmt(fftMs, 4) << std::setw(12) << fmt(speedup, 1) << "x\n";
        } else {
            std::cout << std::setw(9) << n << std::setw(12) << "-" << std::setw(13)
                      << fmt(fftMs, 4) << std::setw(12) << "-" << "\n";
        }

        rows.push_back({fmtInt(static_cast<long long>(n)),
                        dftMs > 0 ? fmt(dftMs, 3) : "",
                        fmt(fftMs, 4),
                        dftMs > 0 ? fmt(speedup, 1) : ""});
    }

    writeCsv("output/fft_bench.csv", {"n", "dft_ms", "fft_ms", "speedup"}, rows);
    std::cout << "\n结果已写入 output/fft_bench.csv\n\n";
}

void dumpSpectrum() {
    // 一个便于讲解的示例：两个正弦 + 一点噪声，输出幅度谱
    const std::size_t n = 256;
    Signal x(n);
    Lcg rng(777ULL);
    for (std::size_t i = 0; i < n; ++i) {
        const double t = static_cast<double>(i);
        const double v = std::cos(kTwoPi * 10.0 * t / n) + 0.5 * std::cos(kTwoPi * 33.0 * t / n) +
                         0.1 * rng.next();
        x[i] = Complex(v, 0.0);
    }
    const Signal X = fft(x);

    std::vector<CsvRow> rows;
    for (std::size_t k = 0; k <= n / 2; ++k) {
        rows.push_back({fmtInt(static_cast<long long>(k)), fmt(std::abs(X[k]), 6)});
    }
    writeCsv("output/fft_spectrum.csv", {"bin", "magnitude"}, rows);
    std::cout << "示例信号幅度谱已写入 output/fft_spectrum.csv（峰值应在 bin 10 与 33）\n";
}

}  // namespace

int main() {
    enableUtf8Console();
    try {
        std::cout << "基 2 FFT 演示（阶段一）\n\n";
        verifyCorrectness();
        benchmark();
        dumpSpectrum();
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "错误: " << e.what() << "\n";
        return 1;
    }
}
