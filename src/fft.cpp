#include "dsp/fft.h"

#include <cmath>
#include <stdexcept>
#include <string>

namespace dsp {

namespace {

constexpr double kTwoPi = 6.283185307179586476925286766559;

// 旋转因子 W_N^k = exp(-j*2*pi*k/N)
Complex twiddle(std::size_t k, std::size_t n, bool inverse) {
    const double sign = inverse ? 1.0 : -1.0;
    const double angle = sign * kTwoPi * static_cast<double>(k) / static_cast<double>(n);
    return Complex(std::cos(angle), std::sin(angle));
}

void requirePowerOfTwo(std::size_t n) {
    if (!isPowerOfTwo(n)) {
        throw std::invalid_argument("FFT 长度必须是 2 的幂，当前为 " + std::to_string(n));
    }
}

// 频域抽取（DIF）核心：先把序列对折成两个长度减半的子问题，再递归。
// 输出是位反转序，由 fftRecursiveDIF 负责整理。
Signal dif(const Signal& x) {
    const std::size_t n = x.size();
    if (n == 1) return x;

    const std::size_t half = n / 2;
    Signal even(half), odd(half);
    for (std::size_t k = 0; k < half; ++k) {
        const Complex a = x[k];
        const Complex b = x[k + half];
        even[k] = a + b;
        odd[k] = (a - b) * twiddle(k, n, false);
    }
    const Signal e = dif(even);
    const Signal o = dif(odd);

    Signal out(n);
    for (std::size_t k = 0; k < half; ++k) {
        out[k] = e[k];
        out[k + half] = o[k];
    }
    return out;
}

}  // namespace

bool isPowerOfTwo(std::size_t n) { return n != 0 && (n & (n - 1)) == 0; }

Signal dft(const Signal& x) {
    const std::size_t n = x.size();
    Signal out(n);
    for (std::size_t k = 0; k < n; ++k) {
        Complex acc(0.0, 0.0);
        for (std::size_t t = 0; t < n; ++t) {
            acc += x[t] * twiddle(k * t % n, n, false);
        }
        out[k] = acc;
    }
    return out;
}

Signal idft(const Signal& X) {
    const std::size_t n = X.size();
    Signal out(n);
    const double inv_n = 1.0 / static_cast<double>(n);
    for (std::size_t k = 0; k < n; ++k) {
        Complex acc(0.0, 0.0);
        for (std::size_t j = 0; j < n; ++j) {
            acc += X[j] * twiddle(k * j % n, n, true);
        }
        out[k] = acc * inv_n;
    }
    return out;
}

std::vector<std::size_t> bitReversalOrder(std::size_t n) {
    // 求 log2(n)，等于蝶形的级数
    std::size_t bits = 0;
    while ((static_cast<std::size_t>(1) << bits) < n) ++bits;

    std::vector<std::size_t> order(n);
    for (std::size_t i = 0; i < n; ++i) {
        std::size_t rev = 0;
        for (std::size_t b = 0; b < bits; ++b) {
            if (i & (static_cast<std::size_t>(1) << b)) {
                rev |= static_cast<std::size_t>(1) << (bits - 1 - b);
            }
        }
        order[i] = rev;
    }
    return order;
}

Signal fft(const Signal& x) {
    const std::size_t n = x.size();
    requirePowerOfTwo(n);

    // 第一步：按位反转重排。递归的“按奇偶分治”展开到最底层之后，元素的
    // 排列顺序恰好就是下标的位反转序，因此可以一次重排 + 逐级蝶形代替递归。
    const auto order = bitReversalOrder(n);
    Signal a(n);
    for (std::size_t i = 0; i < n; ++i) a[i] = x[order[i]];

    // 第二步：逐级蝶形。第 s 级的蝶形跨度为 2^s，组内两个输入相距 half，
    // 共用同一个旋转因子，一次算出上下两个输出（这就是“蝶形”的形状）。
    for (std::size_t span = 2; span <= n; span <<= 1) {
        const std::size_t half = span >> 1;
        for (std::size_t start = 0; start < n; start += span) {
            for (std::size_t k = 0; k < half; ++k) {
                const Complex w = twiddle(k, span, false);
                const Complex even = a[start + k];
                const Complex odd = a[start + k + half] * w;
                a[start + k] = even + odd;
                a[start + k + half] = even - odd;
            }
        }
    }
    return a;
}

Signal fftRecursiveDIF(const Signal& x) {
    const std::size_t n = x.size();
    requirePowerOfTwo(n);
    if (n == 1) return x;

    const Signal scrambled = dif(x);
    // 整理位反转序 -> 自然序
    const auto order = bitReversalOrder(n);
    Signal out(n);
    for (std::size_t i = 0; i < n; ++i) out[i] = scrambled[order[i]];
    return out;
}

Signal ifft(const Signal& X) {
    const std::size_t n = X.size();
    requirePowerOfTwo(n);

    // IDFT(X) = conj(DFT(conj(X))) / N，于是可以复用正向 FFT。
    Signal conj_in(n);
    for (std::size_t i = 0; i < n; ++i) conj_in[i] = std::conj(X[i]);

    Signal y = fft(conj_in);

    const double inv_n = 1.0 / static_cast<double>(n);
    for (auto& v : y) v = std::conj(v) * inv_n;
    return y;
}

}  // namespace dsp
