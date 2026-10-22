// 阶段一测试：DFT / FFT 的数值正确性。
//
// 设计思路：不依赖任何外部参考实现，用三种自证手段验证：
//   1. 慢而不易错的朴素 DFT 作为基准真值，FFT 必须与它一致；
//   2. 两条互相独立的 FFT 路径（迭代 DIT / 递归 DIF）交叉印证；
//   3. 已知闭式频谱的解析结果（冲激、直流、单频、Parseval 定理）。

#define DOCTEST_CONFIG_IMPLEMENT
#include <doctest.h>

#include <cmath>
#include <complex>
#include <cstdint>
#include <stdexcept>
#include <vector>

#include "dsp/console.h"
#include "dsp/fft.h"
#include "dsp/filter.h"
#include "dsp/qpsk.h"
#include "dsp/window.h"

using namespace dsp;

namespace {

constexpr double kTwoPi = 6.283185307179586476925286766559;

// 可复现的伪随机数，避免测试结果随运行变化
class Lcg {
public:
    explicit Lcg(std::uint64_t seed) : s_(seed) {}
    double next() {
        s_ = s_ * 6364136223846793005ULL + 1442695040888963407ULL;
        const std::uint64_t hi = s_ >> 11;
        return static_cast<double>(hi) / static_cast<double>(1ULL << 53) * 2.0 - 1.0;
    }

private:
    std::uint64_t s_;
};

Signal randomSignal(std::size_t n, std::uint64_t seed = 987654321ULL) {
    Lcg rng(seed);
    Signal x(n);
    for (std::size_t i = 0; i < n; ++i) x[i] = Complex(rng.next(), rng.next());
    return x;
}

double maxAbsDiff(const Signal& a, const Signal& b) {
    REQUIRE(a.size() == b.size());
    double m = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i) m = std::max(m, std::abs(a[i] - b[i]));
    return m;
}

double energy(const Signal& x) {
    double e = 0.0;
    for (const auto& v : x) e += std::norm(v);
    return e;
}

std::size_t argmaxAbs(const Signal& x) {
    std::size_t best = 0;
    for (std::size_t i = 1; i < x.size(); ++i) {
        if (std::abs(x[i]) > std::abs(x[best])) best = i;
    }
    return best;
}

}  // namespace

// 自定义 main：先切 UTF-8 控制台编码，否则中文用例名在 Windows 终端显示为乱码
int main(int argc, char** argv) {
    dsp::enableUtf8Console();
    doctest::Context ctx;
    ctx.applyCommandLine(argc, argv);
    return ctx.run();
}

TEST_CASE("isPowerOfTwo 判定") {
    CHECK(isPowerOfTwo(1));
    CHECK(isPowerOfTwo(2));
    CHECK(isPowerOfTwo(1024));
    CHECK_FALSE(isPowerOfTwo(0));
    CHECK_FALSE(isPowerOfTwo(3));
    CHECK_FALSE(isPowerOfTwo(1000));
}

TEST_CASE("bitReversalOrder 已知小样例") {
    const auto n8 = bitReversalOrder(8);
    const std::vector<std::size_t> expect8 = {0, 4, 2, 6, 1, 5, 3, 7};
    CHECK(n8 == expect8);

    const auto n4 = bitReversalOrder(4);
    const std::vector<std::size_t> expect4 = {0, 2, 1, 3};
    CHECK(n4 == expect4);

    const auto n2 = bitReversalOrder(2);
    const std::vector<std::size_t> expect2 = {0, 1};
    CHECK(n2 == expect2);
}

TEST_CASE("位反转表是一个双射") {
    for (std::size_t n : {2u, 4u, 8u, 16u, 64u, 256u}) {
        const auto order = bitReversalOrder(n);
        std::vector<int> seen(n, 0);
        for (auto v : order) {
            REQUIRE(v < n);
            ++seen[v];
        }
        for (std::size_t i = 0; i < n; ++i) CHECK(seen[i] == 1);
    }
}

TEST_CASE("非 2 的幂长度应当抛异常") {
    const Signal x(6, Complex(1.0, 0.0));
    CHECK_THROWS_AS(fft(x), std::invalid_argument);
    CHECK_THROWS_AS(fftRecursiveDIF(x), std::invalid_argument);
    CHECK_THROWS_AS(ifft(x), std::invalid_argument);
}

TEST_CASE("FFT 与朴素 DFT 结果一致") {
    for (std::size_t n : {2u, 4u, 8u, 16u, 128u, 1024u}) {
        const Signal x = randomSignal(n);
        const double err = maxAbsDiff(dft(x), fft(x));
        CHECK(err < 1e-9);
    }
}

TEST_CASE("迭代 DIT 与递归 DIF 两条独立路径一致") {
    for (std::size_t n : {2u, 8u, 64u, 512u, 2048u}) {
        const Signal x = randomSignal(n, 424242ULL);
        const double err = maxAbsDiff(fft(x), fftRecursiveDIF(x));
        CHECK(err < 1e-9);
    }
}

TEST_CASE("Parseval 定理：时域能量等于频域能量除以 N") {
    for (std::size_t n : {16u, 256u, 1024u}) {
        const Signal x = randomSignal(n, 555ULL);
        const Signal X = fft(x);
        const double et = energy(x);
        const double ef = energy(X) / static_cast<double>(n);
        CHECK(std::abs(et - ef) < 1e-9 * std::max(1.0, et));
    }
}

TEST_CASE("FFT 与 IFFT 往返可还原") {
    for (std::size_t n : {2u, 8u, 256u, 1024u}) {
        const Signal x = randomSignal(n, 31337ULL);
        const double err = maxAbsDiff(x, ifft(fft(x)));
        CHECK(err < 1e-9);
    }
}

TEST_CASE("IDFT 与 IFFT 结果一致") {
    const std::size_t n = 128;
    const Signal x = randomSignal(n);
    const Signal X = dft(x);
    CHECK(maxAbsDiff(idft(X), ifft(X)) < 1e-9);
}

TEST_CASE("单位冲激的频谱恒为 1") {
    const std::size_t n = 64;
    Signal x(n, Complex(0.0, 0.0));
    x[0] = Complex(1.0, 0.0);
    const Signal X = fft(x);
    for (const auto& v : X) CHECK(std::abs(v - Complex(1.0, 0.0)) < 1e-12);
}

TEST_CASE("直流信号的频谱只有第 0 点非零") {
    const std::size_t n = 64;
    const double amp = 2.0;
    Signal x(n, Complex(amp, 0.0));
    const Signal X = fft(x);

    CHECK(std::abs(X[0] - Complex(amp * n, 0.0)) < 1e-9);
    for (std::size_t k = 1; k < n; ++k) CHECK(std::abs(X[k]) < 1e-9);
}

TEST_CASE("单频正弦的频谱峰值落在对应频点") {
    const std::size_t n = 256;
    for (std::size_t k0 : {1u, 19u, 64u, 127u}) {
        Signal x(n);
        for (std::size_t i = 0; i < n; ++i) {
            x[i] = Complex(std::cos(kTwoPi * static_cast<double>(k0) * i / n), 0.0);
        }
        const Signal X = fft(x);
        // 实余弦的频谱在 k0 与 N-k0 是一对共轭对称的峰，幅度同为 N/2。
        // 两者在浮点意义下可能相差最后几位，因此峰值落在任一个上都算正确。
        const std::size_t peak = argmaxAbs(X);
        CHECK((peak == k0 || peak == n - k0));
        CHECK(std::abs(std::abs(X[k0]) - n / 2.0) < 1e-9);
        CHECK(std::abs(std::abs(X[n - k0]) - n / 2.0) < 1e-9);
    }
}

TEST_CASE("实信号的频谱满足共轭对称") {
    const std::size_t n = 256;
    Signal x(n);
    for (std::size_t i = 0; i < n; ++i) {
        x[i] = Complex(std::sin(0.1 * i) + 0.3 * std::cos(0.37 * i), 0.0);
    }
    const Signal X = fft(x);
    for (std::size_t k = 1; k < n / 2; ++k) {
        CHECK(std::abs(X[k] - std::conj(X[n - k])) < 1e-9);
    }
}

TEST_CASE("线性性质：FFT 是线性变换") {
    const std::size_t n = 128;
    const Signal a = randomSignal(n, 1ULL);
    const Signal b = randomSignal(n, 2ULL);
    const double alpha = 1.7, beta = -0.4;

    Signal mix(n);
    for (std::size_t i = 0; i < n; ++i) mix[i] = alpha * a[i] + beta * b[i];

    const Signal lhs = fft(mix);
    const Signal A = fft(a), B = fft(b);
    Signal rhs(n);
    for (std::size_t i = 0; i < n; ++i) rhs[i] = alpha * A[i] + beta * B[i];

    CHECK(maxAbsDiff(lhs, rhs) < 1e-9);
}

TEST_CASE("时移对应频域乘旋转因子") {
    const std::size_t n = 64;
    const Signal x = randomSignal(n, 99ULL);
    const std::size_t shift = 7;

    Signal shifted(n);
    for (std::size_t i = 0; i < n; ++i) shifted[(i + shift) % n] = x[i];

    const Signal X = fft(x), Y = fft(shifted);
    for (std::size_t k = 0; k < n; ++k) {
        const Complex w = std::exp(Complex(0.0, -kTwoPi * static_cast<double>(k) *
                                                    static_cast<double>(shift) / n));
        CHECK(std::abs(Y[k] - X[k] * w) < 1e-9);
    }
}

// ==================== 阶段二：窗函数与滤波器 ====================

TEST_CASE("窗函数：矩形窗恒为 1，其余窗关于中心对称") {
    for (std::size_t n : {2u, 5u, 16u, 65u, 256u}) {
        const auto r = window(WindowType::Rectangular, n);
        REQUIRE(r.size() == n);
        for (double v : r) CHECK(v == 1.0);

        for (WindowType t : {WindowType::Hann, WindowType::Hamming, WindowType::Blackman}) {
            const auto w = window(t, n);
            REQUIRE(w.size() == n);
            for (std::size_t i = 0; i < n; ++i) {
                CHECK(std::abs(w[i] - w[n - 1 - i]) < 1e-12);
            }
        }
    }
}

TEST_CASE("窗函数端点值与中心值符合定义") {
    const std::size_t n = 64;
    CHECK(std::abs(window(WindowType::Hann, n)[0] - 0.0) < 1e-12);
    CHECK(std::abs(window(WindowType::Hamming, n)[0] - 0.08) < 1e-12);
    CHECK(std::abs(window(WindowType::Blackman, n)[0] - 0.0) < 1e-12);

    // 奇数长度时中心正好落在采样点上，三种窗的中心值都等于 1
    const std::size_t odd = 65;
    CHECK(std::abs(window(WindowType::Hann, odd)[32] - 1.0) < 1e-12);
    CHECK(std::abs(window(WindowType::Hamming, odd)[32] - 1.0) < 1e-12);
    CHECK(std::abs(window(WindowType::Blackman, odd)[32] - 1.0) < 1e-12);
}

TEST_CASE("FIR 低通：直流增益为 1、系数对称、群延迟正确") {
    const std::size_t taps = 65;
    for (WindowType t : {WindowType::Rectangular, WindowType::Hann, WindowType::Hamming,
                         WindowType::Blackman}) {
        const auto h = designFirLowpass(0.1, taps, t);
        REQUIRE(h.size() == taps);

        double sum = 0.0;
        for (double v : h) sum += v;
        CHECK(std::abs(sum - 1.0) < 1e-12);   // 直流增益归一化

        for (std::size_t i = 0; i < taps; ++i) {
            CHECK(std::abs(h[i] - h[taps - 1 - i]) < 1e-12);   // 线性相位 => 系数偶对称
        }
        CHECK(std::abs(groupDelay(taps) - 32.0) < 1e-12);
    }
}

TEST_CASE("FIR 频响：直流处为 1，奈奎斯特处被压得很低") {
    const auto h = designFirLowpass(0.1, 65, WindowType::Hamming);
    const auto H = freqResponse(h, 257);
    CHECK(std::abs(H.front() - Complex(1.0, 0.0)) < 1e-9);   // ω = 0
    CHECK(std::abs(H.back()) < 1e-3);                        // ω = π
}

TEST_CASE("FIR 卷积：冲激输入得到系数本身") {
    const auto h = designFirLowpass(0.2, 31, WindowType::Hann);
    Signal imp(1, Complex(1.0, 0.0));
    const Signal y = firFilter(h, imp);
    REQUIRE(y.size() == h.size());
    for (std::size_t i = 0; i < h.size(); ++i) {
        CHECK(std::abs(y[i] - Complex(h[i], 0.0)) < 1e-12);
    }
}

TEST_CASE("FIR 低通对直流信号无失真（稳态增益为 1）") {
    const auto h = designFirLowpass(0.1, 65, WindowType::Hamming);
    const std::size_t n = 200;
    Signal x(n, Complex(1.0, 0.0));
    const Signal y = firFilter(h, x);
    for (std::size_t i = h.size(); i < n; ++i) {
        CHECK(std::abs(y[i] - Complex(1.0, 0.0)) < 1e-9);
    }
}

TEST_CASE("IIR 双二阶：直流增益为 1") {
    const Biquad bq = biquadLowpass(0.1, 0.70710678118654752);
    const auto H = biquadResponse(bq, 129);
    CHECK(std::abs(H.front() - Complex(1.0, 0.0)) < 1e-9);

    const std::size_t n = 400;
    Signal x(n, Complex(1.0, 0.0));
    const Signal y = biquadFilter(bq, x);
    for (std::size_t i = 200; i < n; ++i) {
        CHECK(std::abs(y[i] - Complex(1.0, 0.0)) < 1e-6);
    }
}

TEST_CASE("IIR 冲激响应最终衰减到 0（说明滤波器稳定）") {
    const Biquad bq = biquadLowpass(0.1, 0.70710678118654752);
    const std::size_t n = 256;
    Signal imp(n, Complex(0.0, 0.0));
    imp[0] = Complex(1.0, 0.0);
    const Signal ir = biquadFilter(bq, imp);

    double tail = 0.0;
    for (std::size_t i = 64; i < n; ++i) tail = std::max(tail, std::abs(ir[i]));
    CHECK(tail < 1e-9);
}

// ==================== 阶段三：QPSK 调制解调 ====================

TEST_CASE("QPSK 映射/判决往返一致，星座点为能量归一") {
    for (int b0 = 0; b0 <= 1; ++b0) {
        for (int b1 = 0; b1 <= 1; ++b1) {
            const Complex s = qpskMap(b0, b1);
            CHECK(std::abs(std::abs(s) - 1.0) < 1e-12);   // 单位能量

            int r0 = -1, r1 = -1;
            qpskDemap(s, r0, r1);
            CHECK(r0 == b0);
            CHECK(r1 == b1);
        }
    }
}

TEST_CASE("QPSK 格雷映射：差 1 比特的码点恰是相邻星座点") {
    struct Code { int b0; int b1; };
    const Code codes[4] = {{0, 0}, {0, 1}, {1, 1}, {1, 0}};

    for (int i = 0; i < 4; ++i) {
        for (int j = i + 1; j < 4; ++j) {
            const int hamming = (codes[i].b0 != codes[j].b0 ? 1 : 0) +
                                (codes[i].b1 != codes[j].b1 ? 1 : 0);
            const double dist = std::abs(qpskMap(codes[i].b0, codes[i].b1) -
                                         qpskMap(codes[j].b0, codes[j].b1));
            if (hamming == 1) {
                CHECK(std::abs(dist - std::sqrt(2.0)) < 1e-12);   // 相邻
            } else {
                CHECK(std::abs(dist - 2.0) < 1e-12);              // 对角（差 2 比特）
            }
        }
    }
}

TEST_CASE("RRC 滤波器：单位能量、对称、峰值在中心") {
    for (double beta : {0.22, 0.35, 0.5}) {
        for (std::size_t sps : {4u, 8u}) {
            const std::size_t span = 6;
            const auto h = rrcFilter(beta, span, sps);
            const std::size_t len = 2 * span * sps + 1;
            REQUIRE(h.size() == len);

            double e = 0.0;
            for (double v : h) e += v * v;
            CHECK(std::abs(e - 1.0) < 1e-12);

            for (std::size_t i = 0; i < len; ++i) CHECK(std::abs(h[i] - h[len - 1 - i]) < 1e-12);

            const std::size_t center = span * sps;
            for (std::size_t i = 0; i < len; ++i) CHECK(h[center] >= h[i] - 1e-12);
        }
    }
}

TEST_CASE("上采样把符号放到每隔 sps 的位置") {
    const Signal sym = {Complex(1.0, 0.0), Complex(2.0, 0.0), Complex(3.0, 0.0)};
    const std::size_t sps = 3;
    const Signal up = upsample(sym, sps);
    REQUIRE(up.size() == 9);
    for (std::size_t i = 0; i < 9; ++i) {
        if (i % sps == 0) {
            CHECK(std::abs(up[i] - sym[i / sps]) < 1e-12);
        } else {
            CHECK(std::abs(up[i]) < 1e-12);
        }
    }
}

TEST_CASE("加高斯白噪声：零方差不动，方差与设定相符") {
    Signal x(1000, Complex(0.0, 0.0));
    addAwgn(x, 0.0, 1ULL);
    for (const auto& v : x) CHECK(std::abs(v) < 1e-15);

    const double var = 0.5;
    Signal y(100000, Complex(0.0, 0.0));
    addAwgn(y, var, 20260920ULL);

    double mean = 0.0;
    for (const auto& v : y) mean += v.real();
    mean /= static_cast<double>(y.size());
    double m2 = 0.0;
    for (const auto& v : y) m2 += (v.real() - mean) * (v.real() - mean);
    const double sampleVar = m2 / static_cast<double>(y.size() - 1);

    CHECK(std::abs(sampleVar - var) < 0.05 * var);
}

TEST_CASE("无噪链路：误码率为 0，判决幅度接近单位能量") {
    const QpskLinkResult r = runQpskLink(300.0, 40000, 0.35, 8, 6, 20260920ULL);

    std::size_t errors = 0;
    for (std::size_t i = 0; i < r.txBits.size(); ++i) {
        if (r.txBits[i] != r.rxBits[i]) ++errors;
    }
    CHECK(errors == 0);

    // 幅度约束能同时验证 RRC 的能量归一与成形/匹配/抽样对齐是正确的
    for (const auto& s : r.decisionSymbols) CHECK(std::abs(std::abs(s) - 1.0) < 0.1);
}

TEST_CASE("实测误码率贴近理论曲线") {
    const double theory = qpskTheoreticalBer(4.0);
    const double measured = simulateQpskBer(4.0, 200000, 0.35, 8, 6, 20260920ULL);
    CHECK(measured > 0.6 * theory);
    CHECK(measured < 1.5 * theory);
}

TEST_CASE("误码率随 Eb/N0 增大而下降") {
    const double b0 = simulateQpskBer(0.0, 50000, 0.35, 8, 6, 1ULL);
    const double b4 = simulateQpskBer(4.0, 50000, 0.35, 8, 6, 2ULL);
    const double b8 = simulateQpskBer(8.0, 50000, 0.35, 8, 6, 3ULL);
    CHECK(b0 > b4);
    CHECK(b4 > b8);
}

TEST_CASE("理论误码率随信噪比单调下降") {
    double prev = 1.0;
    for (double db = 0.0; db <= 10.0; db += 1.0) {
        const double ber = qpskTheoreticalBer(db);
        CHECK(ber < prev);
        prev = ber;
    }
    // 与已知量级核对：4 dB 时 QPSK 理论 BER ≈ 1.19e-2
    CHECK(std::abs(qpskTheoreticalBer(4.0) - 0.0119) < 0.001);
}

