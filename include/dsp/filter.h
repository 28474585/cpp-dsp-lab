#pragma once

#include <complex>
#include <cstddef>
#include <vector>

#include "dsp/fft.h"     // Complex / Signal
#include "dsp/window.h"

namespace dsp {

// ============================ FIR 滤波器 ============================
//
// FIR（有限冲激响应）：输出只依赖当前和过去的输入，没有反馈。
// 好处是天然稳定、可以做到严格线性相位（群延迟恒定，波形不失真）；
// 代价是达到同样的选择性需要更高的阶数。

// 窗函数法设计线性相位低通 FIR。
//   cutoff  : 归一化截止频率，以采样率为单位，取值 (0, 0.5)，0.5 即奈奎斯特频率
//   numTaps : 抽头数。奇数得到整数群延迟、严格线性相位（第 I 类线性相位）
//   win     : 截断用的窗函数，决定旁瓣高低与过渡带宽度
// 返回系数已归一化到直流增益 sum(h) = 1，即通带增益为 1。
std::vector<double> designFirLowpass(double cutoff, std::size_t numTaps, WindowType win);

// FIR 直接卷积滤波（线性卷积，不做循环环绕）。
// 输出长度 = x.size() + h.size() - 1，群延迟为 groupDelay(h.size()) 个采样点。
Signal firFilter(const std::vector<double>& h, const Signal& x);

// 频率响应 H(e^{jω})，在 ω ∈ [0, π] 上均匀取 points 个点（对应数字频率 0 到奈奎斯特）。
// 取 abs() 即幅度谱，可用来观察通带、阻带和窗函数带来的旁瓣差异。
std::vector<Complex> freqResponse(const std::vector<double>& h, std::size_t points);

// 线性相位 FIR 的群延迟（采样点数），恒为 (numTaps-1)/2。
double groupDelay(std::size_t numTaps);

// ============================ IIR 滤波器 ============================
//
// IIR（无限冲激响应）：带反馈，用很低的阶数就能做出陡峭的过渡带，
// 计算量和内存都远小于同指标的 FIR。代价是相位非线性、且必须关心稳定性
//（极点要在单位圆内）。

// 二阶（双二阶 / biquad）节：
//   y[n] = b0·x[n] + b1·x[n-1] + b2·x[n-2] − a1·y[n-1] − a2·y[n-2]
// 系数均已对 a0 归一化。
struct Biquad {
    double b0 = 1.0;
    double b1 = 0.0;
    double b2 = 0.0;
    double a1 = 0.0;
    double a2 = 0.0;
};

// 二阶低通 IIR 设计（RBJ 双二阶公式，即音频领域广为使用的 cookbook 形式）。
//   cutoff : 归一化截止频率，以采样率为单位，取值 (0, 0.5)
//   q      : 品质因数。巴特沃斯（最平坦）响应取 1/√2 ≈ 0.7071；
//            Q 越大通带越接近截止频率处越"翘"，产生谐振峰
// 直流增益恒为 1。
Biquad biquadLowpass(double cutoff, double q);

// IIR 直接 I 型实现：内部维护 x/y 的延时线。输出长度与输入相同。
Signal biquadFilter(const Biquad& bq, const Signal& x);

// 双二阶的频率响应 H(e^{jω})，ω ∈ [0, π] 上均匀取 points 个点。
std::vector<Complex> biquadResponse(const Biquad& bq, std::size_t points);

}  // namespace dsp
