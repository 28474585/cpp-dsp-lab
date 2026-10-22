#!/usr/bin/env python3
"""把 output/ 下的 CSV 画成图。

这些 CSV 是 C++ 程序跑出来的（fft_demo / filter_demo / qpsk_demo）。
本脚本只负责可视化，不做任何计算——项目本身不依赖 Python，是可选的辅助工具。

用法：
    python tools/plot.py                 # 画全部（缺哪个 CSV 就跳过哪个）
    python tools/plot.py --only ber      # 只画 BER 曲线
    python tools/plot.py --only filter   # 只画滤波器组图

依赖：numpy、matplotlib
"""

import argparse
import csv
import os
import sys

# 让中文提示在 Windows 控制台正常显示（与 C++ 侧 CP_UTF8 的做法对齐）
try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

import matplotlib

matplotlib.use("Agg")  # 无显示器环境也能出图
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "output")
FIG = os.path.join(OUT, "figures")

# 让中文标题正常显示（Windows 上优先用微软雅黑，退回黑体）
for _fam in ("Microsoft YaHei", "SimHei", "DejaVu Sans"):
    try:
        matplotlib.font_manager.findfont(_fam, fallback_to_default=False)
        plt.rcParams["font.sans-serif"] = [_fam]
        break
    except Exception:
        continue
plt.rcParams["axes.unicode_minus"] = False


def load(name):
    """读一个 CSV，返回列名与按列组成的 numpy 数组。文件不存在返回 None。"""
    path = os.path.join(OUT, name)
    if not os.path.exists(path):
        print(f"  [跳过] 找不到 {name}")
        return None
    with open(path, newline="", encoding="utf-8") as f:
        reader = csv.reader(f)
        header = next(reader)
        rows = [r for r in reader if r]
    if not rows:
        print(f"  [跳过] {name} 没有数据")
        return None
    cols = np.array(rows)
    data = {}
    for i, key in enumerate(header):
        try:
            data[key] = cols[:, i].astype(float)
        except ValueError:
            data[key] = cols[:, i]
    return data


def save(fig, name):
    os.makedirs(FIG, exist_ok=True)
    path = os.path.join(FIG, name)
    fig.savefig(path, dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(f"  已保存 {os.path.relpath(path, ROOT)}")


def plot_spectrum():
    d = load("fft_spectrum.csv")
    if d is None:
        return
    fig, ax = plt.subplots(figsize=(9, 4))
    ax.plot(d["bin"], d["magnitude"], lw=1.0)
    ax.set_title("阶段一 · 示例信号幅度谱（峰值应在 bin 10 与 33）")
    ax.set_xlabel("频点 bin")
    ax.set_ylabel("幅度")
    ax.grid(alpha=0.3)
    save(fig, "01-fft-spectrum.png")


def plot_bench():
    d = load("fft_bench.csv")
    if d is None:
        return
    fig, ax = plt.subplots(figsize=(9, 4))
    n = d["n"]
    ax.loglog(n, d["dft_ms"], "o-", label="朴素 DFT  O(N²)")
    ax.loglog(n, d["fft_ms"], "s-", label="FFT  O(N log N)")
    ax.set_title("阶段一 · DFT 与 FFT 耗时随长度的变化（双对数）")
    ax.set_xlabel("N")
    ax.set_ylabel("耗时 (ms)")
    ax.grid(alpha=0.3, which="both")
    ax.legend()
    save(fig, "02-fft-bench.png")


def plot_freq_response():
    d = load("filter_freq_response.csv")
    if d is None:
        return
    fig, ax = plt.subplots(figsize=(10, 5))
    f = d["freq_norm"]
    for key, label in [
        ("rect_db", "矩形窗 FIR"),
        ("hann_db", "汉宁窗 FIR"),
        ("hamming_db", "汉明窗 FIR"),
        ("blackman_db", "布莱克曼窗 FIR"),
        ("iir_db", "IIR biquad（2 阶）"),
    ]:
        if key in d:
            ax.plot(f, d[key], lw=1.2, label=label)
    ax.set_ylim(-120, 10)
    ax.set_xlim(0, 0.5)
    ax.axhline(-60, color="gray", ls=":", lw=0.8, label="-60 dB 参考")
    ax.set_title("阶段二 · 低通滤波器幅频响应（截止 0.1）")
    ax.set_xlabel("归一化频率（1.0 = 采样率）")
    ax.set_ylabel("幅度 (dB)")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)
    save(fig, "03-filter-freq-response.png")


def plot_time_domain():
    d = load("filter_time_domain.csv")
    if d is None:
        return
    fig, ax = plt.subplots(figsize=(10, 4))
    ax.plot(d["n"], d["input"], lw=0.8, alpha=0.6, label="滤波前（低频 + 高频）")
    ax.plot(d["n"], d["output"], lw=1.4, label="滤波后（只剩低频，已对齐延迟）")
    ax.set_title("阶段二 · 低通滤波前后时域对比")
    ax.set_xlabel("采样点 n")
    ax.set_ylabel("幅度")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=9)
    save(fig, "04-filter-time-domain.png")


def plot_leakage():
    d = load("spectral_leakage.csv")
    if d is None:
        return
    fig, ax = plt.subplots(figsize=(10, 4.5))
    ax.plot(d["bin"], d["rect_db"], lw=1.0, label="矩形窗")
    ax.plot(d["bin"], d["hamming_db"], lw=1.0, label="汉明窗")
    ax.set_ylim(-80, 5)
    ax.set_title("阶段二 · 频谱泄漏：非整数周期正弦（10.3 周期）加不同窗")
    ax.set_xlabel("频点 bin")
    ax.set_ylabel("相对主瓣电平 (dB)")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=9)
    save(fig, "05-spectral-leakage.png")


def plot_impulse():
    d = load("filter_impulse_response.csv")
    if d is None:
        return
    fig, ax = plt.subplots(figsize=(9, 4))
    ax.stem(d["n"][:80], d["iir_impulse"][:80], basefmt=" ")
    ax.set_title("阶段二 · IIR 冲激响应（衰减到 0，极点在单位圆内，稳定）")
    ax.set_xlabel("采样点 n")
    ax.set_ylabel("幅度")
    ax.grid(alpha=0.3)
    save(fig, "06-iir-impulse-response.png")


def plot_constellation():
    d = load("qpsk_constellation.csv")
    if d is None:
        return
    fig, ax = plt.subplots(figsize=(5.5, 5.5))
    ax.plot(d["i"], d["q"], ".", ms=3, alpha=0.4)
    ideal = 1.0 / np.sqrt(2.0)
    for sx in (-ideal, ideal):
        for sy in (-ideal, ideal):
            ax.plot([sx], [sy], "r+", ms=14, mew=2)
    ax.axhline(0, color="gray", lw=0.5)
    ax.axvline(0, color="gray", lw=0.5)
    ax.set_title("阶段三 · QPSK 星座图（9 dB，红叉为理想点）")
    ax.set_xlabel("同相 I")
    ax.set_ylabel("正交 Q")
    ax.set_aspect("equal")
    ax.grid(alpha=0.3)
    save(fig, "07-qpsk-constellation.png")


def plot_eye():
    d = load("qpsk_eye.csv")
    if d is None:
        return
    fig, ax = plt.subplots(figsize=(7, 5))
    trace = d["trace"]
    for t in np.unique(trace):
        m = trace == t
        ax.plot(d["offset"][m], d["value"][m], lw=0.5, alpha=0.25, color="steelblue")
    ax.set_title("阶段三 · QPSK 眼图（12 dB，200 条轨迹叠加）")
    ax.set_xlabel("段内采样偏移（1 段 = 2 个符号周期）")
    ax.set_ylabel("幅度 (I 路)")
    ax.grid(alpha=0.3)
    save(fig, "08-qpsk-eye.png")


def plot_ber():
    d = load("qpsk_ber.csv")
    if d is None:
        return
    fig, ax = plt.subplots(figsize=(8.5, 6))
    ax.semilogy(d["eb_n0_db"], d["theoretical_ber"], "k-", lw=1.5,
                label="理论 QPSK：½·erfc(√(Eb/N0))")
    ax.semilogy(d["eb_n0_db"], d["measured_ber"], "ro", ms=6, mfc="none",
                label="实测（蒙特卡洛仿真）")
    ax.set_ylim(1e-5, 0.2)
    ax.set_title("阶段三 · QPSK 误码率实测 vs 理论")
    ax.set_xlabel("Eb/N0 (dB)")
    ax.set_ylabel("误比特率 BER")
    ax.grid(alpha=0.3, which="both")
    ax.legend()
    save(fig, "09-qpsk-ber.png")
    print("  （这条曲线就是全链路的'总检验'：实测点贴合理论，说明调制、成形、加噪、"
          "匹配滤波、判决每一步都对。）")


def main():
    ap = argparse.ArgumentParser(description="把 cpp-dsp-lab 的输出 CSV 画成图")
    ap.add_argument("--only", default="all",
                    choices=["all", "fft", "filter", "ber"],
                    help="只画某一类图")
    args = ap.parse_args()

    print(f"读取目录: {os.path.relpath(OUT, ROOT)}")
    if args.only in ("all", "fft"):
        print("[阶段一] FFT")
        plot_spectrum()
        plot_bench()
    if args.only in ("all", "filter"):
        print("[阶段二] 滤波器与窗函数")
        plot_freq_response()
        plot_time_domain()
        plot_leakage()
        plot_impulse()
    if args.only in ("all", "ber"):
        print("[阶段三] QPSK")
        plot_constellation()
        plot_eye()
        plot_ber()
    print(f"\n图已输出到 {os.path.relpath(FIG, ROOT)}/")
    return 0


if __name__ == "__main__":
    sys.exit(main())
