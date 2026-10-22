#pragma once

#include <complex>
#include <cstddef>
#include <vector>

namespace dsp {

using Complex = std::complex<double>;
using Signal = std::vector<Complex>;

// 朴素离散傅里叶变换，O(N^2)。
// 它唯一的用途是充当 FFT 的基准真值：FFT 的结果必须和它在极小的误差内一致。
Signal dft(const Signal& x);

// 朴素逆变换，O(N^2)。
Signal idft(const Signal& X);

// 判断是否为 2 的幂（基 2 算法的前提）。
bool isPowerOfTwo(std::size_t n);

// 返回 n 的二进制位反转下标表，长度 n。n 必须是 2 的幂。
std::vector<std::size_t> bitReversalOrder(std::size_t n);

// 迭代式基 2 Cooley-Tukey FFT：位反转重排 + 蝶形运算，原地计算，O(N log N)。
// n 必须是 2 的幂，否则抛 std::invalid_argument。
Signal fft(const Signal& x);

// 频率抽取法（DIF）的递归实现，与 fft() 互为独立实现，用于交叉验证。
// 说明：DIT（时域抽取）先做位反转重排再逐级蝶形；DIF（频域抽取）先逐级蝶形，
// 最后得到的输出是位反转序。两条路径结果应当完全一致。
Signal fftRecursiveDIF(const Signal& x);

// 逆变换：利用共轭关系 conj(FFT(conj(X)))/N 实现，O(N log N)。
Signal ifft(const Signal& X);

}  // namespace dsp
