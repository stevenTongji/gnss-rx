#!/usr/bin/env python3
"""冒烟测试：验证环境装对了，并且"捕获"这一步真的能工作。

运行（在仓库根目录）：
    uv run python scripts/00_smoke_test.py

它会做四件事：
    1. 打印 Python / numpy / matplotlib 版本       —— 确认环境
    2. 自检 C/A 码生成器                            —— 确认基础算法正确
    3. 合成一段含 5 颗星的 GPS L1 中频信号          —— 答案已知
    4. 对全部 32 个 PRN 做捕获，和真值比对          —— 确认捕获算法正确

全部通过才会打印 ✅，并输出一张图到 results/00_smoke_test.png。
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")          # 无图形界面环境也能出图
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import gnssrx.plotting  # noqa: F401,E402  —— 注册中文字体配置
from gnssrx.acquisition import acquire_all, detect            # noqa: E402
from gnssrx.ca_code import self_test as ca_self_test          # noqa: E402
from gnssrx.sim import DEFAULT_FS, DEFAULT_IF, synthesize     # noqa: E402

TRUE_PRNS = [7, 11, 19, 24, 30]
FS, F_IF = DEFAULT_FS, DEFAULT_IF
N_MS = 50


def main() -> int:
    print("=" * 62)
    print("① 环境检查")
    print(f"   Python     {sys.version.split()[0]}")
    print(f"   numpy      {np.__version__}")
    print(f"   matplotlib {matplotlib.__version__}")
    print(f"   采样率 {FS/1e6:.3f} MHz | 中频 {F_IF/1e6:.3f} MHz | 时长 {N_MS} ms")

    print("\n② C/A 码生成器自检")
    ca_self_test()

    print("\n③ 合成中频信号（真值已知）")
    sig, truth = synthesize(TRUE_PRNS, fs=FS, f_if=F_IF, n_ms=N_MS, cn0_dbhz=45)
    for t in truth:
        print(f"   PRN {t['prn']:>2}  码延迟 {t['code_delay_samples']:>5} 采样  "
              f"多普勒 {t['doppler_hz']:>+8.1f} Hz")

    print("\n④ 对 32 个 PRN 做并行码相位捕获")
    results = acquire_all(sig, FS, F_IF, ms=1,
                          doppler_half_range=6000, doppler_step=500)
    found = detect(results, threshold=2.5)
    found_prns = sorted(r["prn"] for r in found)

    print(f"   检出 {len(found_prns)} 颗：{found_prns}")
    print(f"   真值 {len(TRUE_PRNS)} 颗：{sorted(TRUE_PRNS)}")

    print("\n   排名前 8（含若干噪声，观察峰/次峰比如何区分）：")
    print(f"   {'PRN':>4} {'峰/次峰':>8} {'峰/均值':>9} {'多普勒Hz':>10} {'码相位(采样)':>12}")
    for r in results[:8]:
        flag = "  ← 真值" if r["prn"] in TRUE_PRNS else ""
        print(f"   {r['prn']:>4} {r['peak_ratio']:>8.2f} {r['peak_over_mean']:>9.2f} "
              f"{r['doppler_hz']:>+10.0f} {r['code_phase_samples']:>12}{flag}")

    # ---- 判分 ----
    ok = True
    if found_prns == sorted(TRUE_PRNS):
        print("\n   ✅ 捕获结果与真值完全一致")
    else:
        ok = False
        print(f"\n   ❌ 不一致：漏检 {set(TRUE_PRNS) - set(found_prns)}，"
              f"虚警 {set(found_prns) - set(TRUE_PRNS)}")

    # 码相位误差应在一个码片（4 个采样点）以内
    for t in truth:
        r = next(x for x in results if x["prn"] == t["prn"])
        err = abs(r["code_phase_samples"] - t["code_delay_samples"])
        if err > FS / 1.023e6:
            ok = False
            print(f"   ❌ PRN {t['prn']} 码相位偏差 {err} 采样，超过 1 码片")
    if ok:
        print("   ✅ 所有卫星码相位误差 < 1 码片")

    # ---- 出图 ----
    out = ROOT / "results" / "00_smoke_test.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    prns = [r["prn"] for r in results]
    ratios = [r["peak_ratio"] for r in results]
    colors = ["#3ddc97" if p in TRUE_PRNS else "#4a5560" for p in prns]

    fig, ax = plt.subplots(figsize=(11, 4.2), dpi=130)
    ax.bar(prns, ratios, color=colors)
    ax.axhline(2.5, ls="--", lw=1.2, color="#ffb454", label="检测门限 2.5")
    ax.set_xlabel("PRN")
    ax.set_ylabel("峰 / 次峰 比")
    ax.set_title("并行码相位 FFT 捕获结果（绿色 = 真实存在的卫星）")
    ax.set_xticks(prns)
    ax.tick_params(axis="x", labelsize=8)
    ax.legend(loc="upper right", fontsize=9)
    ax.grid(axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(out)
    print(f"\n   图已保存：{out}")

    print("=" * 62)
    print("✅ 全部通过，环境可用。" if ok else "❌ 有检查未通过，请看上面输出。")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
