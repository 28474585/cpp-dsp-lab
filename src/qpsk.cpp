#include "dsp/qpsk.h"

#include <cmath>
#include <stdexcept>
#include <string>

#include "dsp/filter.h"   // firFilter：成形与匹配滤波都复用它

namespace dsp {

namespace {

constexpr double kPi = 3.14159265358979323846;
constexpr double kSqrt2 = 1.4142135623730950488016887242097;

// 可复现的线性同余发生器，与阶段一的写法保持一致（不依赖 rand()）。
class Lcg {
public:
    explicit Lcg(std::uint64_t seed) : s_(seed ? seed : 0x9E3779B97F4A7C15ULL) {}

    std::uint64_t nextU64() {
        s_ = s_ * 6364136223846793005ULL + 1442695040888963407ULL;
        return s_;
    }

    // [0, 1)
    double uniform() {
        const std::uint64_t hi = nextU64() >> 11;   // 取高 53 位，低位质量差
        return static_cast<double>(hi) / static_cast<double>(1ULL << 53);
    }

    // 取状态最高位当随机比特（LCG 的高位统计性质最好）
    int bit() { return static_cast<int>(nextU64() >> 63); }

private:
    std::uint64_t s_;
};

// Box-Muller：由两个均匀样本产生一个标准正态样本。
double gaussian(Lcg& rng) {
    double u1 = rng.uniform();
    if (u1 < 1e-300) u1 = 1e-300;   // 防 log(0)
    const double u2 = rng.uniform();
    return std::sqrt(-2.0 * std::log(u1)) * std::cos(2.0 * kPi * u2);
}

}  // namespace

Complex qpskMap(int b0, int b1) {
    const double i = (b1 == 0) ? 1.0 : -1.0;
    const double q = (b0 == 0) ? 1.0 : -1.0;
    return Complex(i / kSqrt2, q / kSqrt2);
}

void qpskDemap(const Complex& symbol, int& b0, int& b1) {
    b1 = (symbol.real() >= 0.0) ? 0 : 1;
    b0 = (symbol.imag() >= 0.0) ? 0 : 1;
}

std::vector<double> rrcFilter(double beta, std::size_t span, std::size_t sps) {
    if (beta <= 0.0 || beta > 1.0) {
        throw std::invalid_argument("beta 必须在 (0, 1] 内，当前为 " + std::to_string(beta));
    }
    if (span == 0) throw std::invalid_argument("span 必须大于 0");
    if (sps == 0) throw std::invalid_argument("sps 必须大于 0");

    const std::size_t len = 2 * span * sps + 1;
    const double center = static_cast<double>(span * sps);

    std::vector<double> h(len);
    for (std::size_t i = 0; i < len; ++i) {
        // 以符号周期为单位的时间变量
        const double t = (static_cast<double>(i) - center) / static_cast<double>(sps);

        double v;
        if (std::abs(t) < 1e-8) {
            // t -> 0 的极限值
            v = 1.0 - beta + 4.0 * beta / kPi;
        } else if (std::abs(std::abs(4.0 * beta * t) - 1.0) < 1e-6) {
            // t = ±1/(4β) 处解析式为 0/0，取洛必达得到的极限。
            // RRC 是偶函数，用 +1/(4β) 处的值即可（注意两侧都要判，否则负侧会 0/0 出 NaN）。
            const double tt = 1.0 / (4.0 * beta);
            const double thm = kPi * (1.0 - beta) * tt;
            const double thp = kPi * (1.0 + beta) * tt;
            const double nPrime = kPi * (1.0 - beta) * std::cos(thm) + 4.0 * beta * std::cos(thp) -
                                  kPi * (1.0 + beta) * std::sin(thp);
            v = nPrime / (-2.0 * kPi);
        } else {
            const double b4t = 4.0 * beta * t;
            const double num = std::sin(kPi * t * (1.0 - beta)) + b4t * std::cos(kPi * t * (1.0 + beta));
            const double den = kPi * t * (1.0 - b4t * b4t);
            v = num / den;
        }
        h[i] = v;
    }

    // 单位能量归一化：这样 (h⊛h) 在零点恰好为 1，匹配滤波后的符号幅度就是 1，
    // 噪声定标才能干净地写成 1/(4·Eb/N0)。
    double e = 0.0;
    for (double v : h) e += v * v;
    const double s = std::sqrt(e);
    for (double& v : h) v /= s;
    return h;
}

Signal upsample(const Signal& symbols, std::size_t sps) {
    if (sps == 0) throw std::invalid_argument("sps 必须大于 0");
    Signal out(symbols.size() * sps, Complex(0.0, 0.0));
    for (std::size_t k = 0; k < symbols.size(); ++k) {
        out[k * sps] = symbols[k];
    }
    return out;
}

void addAwgn(Signal& x, double variancePerComponent, std::uint64_t seed) {
    if (variancePerComponent <= 0.0) return;
    const double sigma = std::sqrt(variancePerComponent);
    Lcg rng(seed);
    for (auto& v : x) {
        v += Complex(sigma * gaussian(rng), sigma * gaussian(rng));
    }
}

double qpskTheoreticalBer(double ebN0Db) {
    // QPSK 与 BPSK 在"每比特信噪比"口径下理论误比特率相同：
    //   BER = Q(√(2·Eb/N0)) = ½·erfc(√(Eb/N0))
    const double ebN0 = std::pow(10.0, ebN0Db / 10.0);
    return 0.5 * std::erfc(std::sqrt(ebN0));
}

QpskLinkResult runQpskLink(double ebN0Db, std::size_t numBits, double beta, std::size_t sps,
                           std::size_t span, std::uint64_t seed) {
    if (numBits < 2) throw std::invalid_argument("numBits 至少为 2");

    const std::size_t nBits = numBits - (numBits % 2);   // 只用偶数个比特（每符号 2 比特）
    const std::size_t nSym = nBits / 2;

    Lcg rng(seed);

    QpskLinkResult r;

    // 1) 产生随机比特并做格雷映射
    r.txBits.resize(nBits);
    for (auto& b : r.txBits) b = rng.bit();

    r.txSymbols.resize(nSym);
    for (std::size_t k = 0; k < nSym; ++k) {
        r.txSymbols[k] = qpskMap(r.txBits[2 * k], r.txBits[2 * k + 1]);
    }

    // 2) 上采样 + 发射端 RRC 成形
    const std::vector<double> h = rrcFilter(beta, span, sps);
    const Signal up = upsample(r.txSymbols, sps);
    r.txShaped = firFilter(h, up);

    // 3) 加高斯白噪声。
    //    符号单位能量 Es = 1，每符号 2 比特 => Es = 2·Eb。
    //    RRC 单位能量，匹配滤波后输出噪声每分量方差 = 每采样点方差。
    //    要求判决点每分量方差 = 1/(4·Eb/N0)，故每采样点方差取同一值。
    const double ebN0 = std::pow(10.0, ebN0Db / 10.0);
    const double varPerComponent = 1.0 / (4.0 * ebN0);

    Signal rx = r.txShaped;
    addAwgn(rx, varPerComponent, seed ^ 0xA5A5A5A5A5A5A5A5ULL);

    // 4) 接收端匹配滤波（还是同一个 RRC）
    r.matchedOutput = firFilter(h, rx);

    // 5) 在符号时刻抽样并判决。
    //    两次卷积的总延迟 = (h.size()-1)；符号 k 的峰值出现在 k*sps + delay。
    const std::size_t delay = h.size() - 1;
    r.decisionSymbols.resize(nSym);
    r.rxBits.resize(nBits);
    for (std::size_t k = 0; k < nSym; ++k) {
        const std::size_t idx = k * sps + delay;
        const Complex s =
            (idx < r.matchedOutput.size()) ? r.matchedOutput[idx] : Complex(0.0, 0.0);
        r.decisionSymbols[k] = s;

        int b0 = 0, b1 = 0;
        qpskDemap(s, b0, b1);
        r.rxBits[2 * k] = b0;
        r.rxBits[2 * k + 1] = b1;
    }
    return r;
}

double simulateQpskBer(double ebN0Db, std::size_t numBits, double beta, std::size_t sps,
                       std::size_t span, std::uint64_t seed) {
    const QpskLinkResult r = runQpskLink(ebN0Db, numBits, beta, sps, span, seed);
    if (r.txBits.empty()) return 0.0;

    std::size_t errors = 0;
    for (std::size_t i = 0; i < r.txBits.size(); ++i) {
        if (r.txBits[i] != r.rxBits[i]) ++errors;
    }
    return static_cast<double>(errors) / static_cast<double>(r.txBits.size());
}

}  // namespace dsp
