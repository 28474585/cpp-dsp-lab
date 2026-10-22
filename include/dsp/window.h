#pragma once

#include <cstddef>
#include <vector>

namespace dsp {

// 窗函数类型。
//
// 一个窗序列有两种用途：
//   1. 设计 FIR 滤波器时，把无限长的理想冲激响应（sinc）截断成有限长，
//      窗的形状决定了截断带来的旁瓣高低；
//   2. 做频谱分析时，对有限长度的采样加窗，抑制因"非整数周期截断"而生的频谱泄漏。
// 两者本质是同一件事：都是在跟"有限长度"这个现实妥协。
enum class WindowType {
    Rectangular,  // 矩形窗：等价于不加窗。主瓣最窄（频率分辨率最好），但旁瓣最高（泄漏最严重）
    Hann,         // 汉宁窗：旁瓣衰减快（-31 dB 起，且快速下降），通用频谱分析用它
    Hamming,      // 汉明窗：第一旁瓣最低（-43 dB），滤波器设计常用
    Blackman,     // 布莱克曼窗：旁瓣最低（-58 dB 起），代价是主瓣最宽
};

// 生成 n 点窗序列。n == 0 返回空，n == 1 返回 {1.0}。
// 采用对称定义（分母 N-1），保证 w[i] == w[n-1-i]，这是线性相位的前提。
std::vector<double> window(WindowType type, std::size_t n);

// 窗的中文名，用于输出。
const char* windowName(WindowType type);

}  // namespace dsp
