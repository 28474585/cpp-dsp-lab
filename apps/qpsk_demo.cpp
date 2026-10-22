// 阶段三演示：QPSK 调制解调全链路 + 误码率实测。
//
// 输出：
//   - 终端：无噪自检、BER 实测 vs 理论对照表
//   - output/qpsk_ber.csv          ：BER-曲线数据（实测与理论）
//   - output/qpsk_constellation.csv：判决点星座图
//   - output/qpsk_eye.csv          ：匹配滤波输出的眼图（I 路）
//
// 用法：qpsk_demo [每点比特数]，默认 200000。比特数越多，低误码率处统计越稳。

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "dsp/console.h"
#include "dsp/csv.h"
#include "dsp/qpsk.h"

using namespace dsp;

namespace {

constexpr double kBeta = 0.35;      // 滚降系数
constexpr std::size_t kSps = 8;     // 每符号采样点数
constexpr std::size_t kSpan = 6;    // RRC 单侧跨度（符号）

// 无噪自检：BER 必须为 0，且判决点应精确落在发送星座点上。
void noiselessCheck() {
    const QpskLinkResult r = runQpskLink(300.0, 40000, kBeta, kSps, kSpan, 20260920ULL);

    std::size_t errors = 0;
    for (std::size_t i = 0; i < r.txBits.size(); ++i) {
        if (r.txBits[i] != r.rxBits[i]) ++errors;
    }

    double maxErr = 0.0;
    for (std::size_t k = 0; k < r.txSymbols.size(); ++k) {
        maxErr = std::max(maxErr, std::abs(r.decisionSymbols[k] - r.txSymbols[k]));
    }

    std::cout << "===== 无噪自检 =====\n";
    std::cout << "误比特数: " << errors << "（必须为 0）\n";
    std::cout << "判决点与发送星座点的最大偏差: " << maxErr
              << "（应接近 0，说明成形/匹配/抽样对齐与能量归一都正确）\n\n";
}

void berSweep(std::size_t numBits) {
    std::cout << "===== 误码率：实测 vs 理论 =====\n";
    std::cout << "每点仿真 " << numBits << " 比特  (beta=" << kBeta << ", sps=" << kSps
              << ", span=" << kSpan << ")\n";
    std::cout << "  Eb/N0(dB)      实测 BER       理论 BER     比值\n";

    std::vector<CsvRow> rows;
    for (double db = 0.0; db <= 8.0001; db += 1.0) {
        const std::uint64_t seed = 20260920ULL + static_cast<std::uint64_t>(db * 10.0);
        const double measured = simulateQpskBer(db, numBits, kBeta, kSps, kSpan, seed);
        const double theory = qpskTheoreticalBer(db);

        // 用科学计数法展示 BER（数量级跨度很大）
        auto sci = [](double v) {
            if (v <= 0.0) return std::string("0");
            std::ostringstream ss;
            ss << std::scientific << std::setprecision(3) << v;
            return ss.str();
        };

        const double ratio = (theory > 0.0) ? measured / theory : 0.0;
        std::cout << std::setw(10) << fmt(db, 1) << std::setw(16) << sci(measured)
                  << std::setw(16) << sci(theory) << std::setw(10) << fmt(ratio, 3) << "\n";

        rows.push_back({fmt(db, 1), sci(measured), sci(theory)});
    }
    writeCsv("output/qpsk_ber.csv", {"eb_n0_db", "measured_ber", "theoretical_ber"}, rows);
    std::cout << "\n实测点应紧贴理论曲线；低 Eb/N0 处比值≈1，高 Eb/N0 处统计误差变大。\n";
    std::cout << "误码率曲线已写入 output/qpsk_ber.csv\n\n";
}

void constellationDemo() {
    // 在 9 dB 下画星座图：四个点应聚成四团，噪声越大越散。
    const QpskLinkResult r = runQpskLink(9.0, 20000, kBeta, kSps, kSpan, 777ULL);

    std::vector<CsvRow> rows;
    const std::size_t limit = std::min<std::size_t>(r.decisionSymbols.size(), 1500);
    for (std::size_t k = 0; k < limit; ++k) {
        rows.push_back({fmt(r.decisionSymbols[k].real(), 6), fmt(r.decisionSymbols[k].imag(), 6)});
    }
    writeCsv("output/qpsk_constellation.csv", {"i", "q"}, rows);
    std::cout << "星座图已写入 output/qpsk_constellation.csv（9 dB，四簇围绕 ±0.707±0.707j）\n\n";
}

void eyeDemo() {
    // 眼图：把匹配滤波输出的 I 路按 2 个符号周期为一段叠画。
    // 眼睛张开的程度即"判决裕量"。
    const QpskLinkResult r = runQpskLink(12.0, 8000, kBeta, kSps, kSpan, 4242ULL);

    const std::vector<double> h = rrcFilter(kBeta, kSpan, kSps);
    const std::size_t delay = h.size() - 1;
    const std::size_t seg = 2 * kSps;      // 每段 2 个符号周期
    const std::size_t traces = 200;        // 叠 200 条

    std::vector<CsvRow> rows;
    for (std::size_t tr = 0; tr < traces; ++tr) {
        const std::size_t start = tr * kSps + delay;
        for (std::size_t j = 0; j < seg; ++j) {
            const std::size_t idx = start + j;
            const double v = (idx < r.matchedOutput.size()) ? r.matchedOutput[idx].real() : 0.0;
            rows.push_back({fmtInt(static_cast<long long>(tr)), fmtInt(static_cast<long long>(j)),
                            fmt(v, 6)});
        }
    }
    writeCsv("output/qpsk_eye.csv", {"trace", "offset", "value"}, rows);
    std::cout << "眼图已写入 output/qpsk_eye.csv（12 dB，200 条轨迹叠加）\n";
}

}  // namespace

int main(int argc, char** argv) {
    enableUtf8Console();
    try {
        std::size_t numBits = 200000;
        if (argc > 1) {
            try {
                numBits = static_cast<std::size_t>(std::stoul(argv[1]));
            } catch (...) {
                std::cerr << "参数解析失败，使用默认 " << numBits << " 比特\n";
            }
        }

        std::cout << "QPSK 调制解调与误码率演示（阶段三）\n\n";
        noiselessCheck();
        berSweep(numBits);
        constellationDemo();
        eyeDemo();
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "错误: " << e.what() << "\n";
        return 1;
    }
}
