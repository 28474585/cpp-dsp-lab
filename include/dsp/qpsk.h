#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

#include "dsp/fft.h"     // Complex / Signal

namespace dsp {

// ======================= 星座映射 =======================

// QPSK 格雷映射：两个比特 -> 一个单位能量星座点。
// 约定（b0 b1，格雷循环）：
//   00 -> (+1,+1)/√2    01 -> (-1,+1)/√2
//   11 -> (-1,-1)/√2    10 -> (+1,-1)/√2
// 即 I 的符号由 b1 决定、Q 的符号由 b0 决定。这样相邻星座点只差一个比特，
// 因此在相近信噪比下"一个符号判错"通常只带来"一个比特错"。
Complex qpskMap(int b0, int b1);

// 最近邻判决（等价于看实部、虚部各自的符号），写回两个比特。
void qpskDemap(const Complex& symbol, int& b0, int& b1);

// ======================= 成形滤波 =======================

// 根升余弦（RRC）滤波器系数，已做单位能量归一化（sum h² = 1）。
//   beta : 滚降系数，0 < β ≤ 1。β 越小频谱越接近理想矩形、越省带宽，
//          但时域拖尾越长、对定时误差越敏感；工程上常取 0.22 ~ 0.5
//   span : 单侧跨越的符号数（决定截断长度，一般 4 ~ 10）
//   sps  : 每符号采样点数（过采样倍数）
// 长度 = 2*span*sps + 1。发射端与接收端各用一个 RRC，级联后得到升余弦响应，
// 在整数倍符号时刻满足无码间串扰（ISI）。
std::vector<double> rrcFilter(double beta, std::size_t span, std::size_t sps);

// 每符号插入 sps-1 个零（上采样），为成形滤波做准备。
Signal upsample(const Signal& symbols, std::size_t sps);

// ======================= 噪声 =======================

// 向信号叠加复高斯白噪声：实部与虚部各自服从 N(0, variancePerComponent)。
// 使用固定种子的线性同余发生器 + Box-Muller 变换，保证同一 seed 结果完全可复现。
// variancePerComponent ≤ 0 时不做任何事（便于构造无噪链路）。
void addAwgn(Signal& x, double variancePerComponent, std::uint64_t seed);

// ======================= 误码率 =======================

// 理论 QPSK 误比特率：BER = Q(√(2·Eb/N0)) = ½·erfc(√(Eb/N0))，输入 Eb/N0（dB）。
double qpskTheoreticalBer(double ebN0Db);

// 一次完整链路的中间结果，供绘图（星座图 / 眼图）与误码统计。
struct QpskLinkResult {
    std::vector<int> txBits;                 // 发送比特
    std::vector<int> rxBits;                 // 判决比特
    Signal txShaped;                         // 成形后（加噪前）的发射信号
    Signal matchedOutput;                    // 匹配滤波后的接收信号（含噪）
    std::vector<Complex> decisionSymbols;    // 符号时刻抽样得到的判决点
    std::vector<Complex> txSymbols;          // 发送的星座点（无噪对照）
};

// 跑一次链路并保留全部中间结果。
//   ebN0Db : 每比特信噪比（dB）。给一个很大的值（如 300）即近似无噪。
// 链路顺序：比特 -> 格雷映射 -> 上采样 -> RRC 成形 -> 加高斯白噪声
//          -> RRC 匹配滤波 -> 符号时刻抽样 -> 判决 -> 得到接收比特。
QpskLinkResult runQpskLink(double ebN0Db, std::size_t numBits, double beta, std::size_t sps,
                           std::size_t span, std::uint64_t seed);

// 同上，但只返回误比特率。
double simulateQpskBer(double ebN0Db, std::size_t numBits, double beta, std::size_t sps,
                       std::size_t span, std::uint64_t seed);

}  // namespace dsp
