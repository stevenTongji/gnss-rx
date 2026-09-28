#!/usr/bin/env python3
"""真实信号验证：在 Nottingham GPS L1 公开数据集上跑通捕获 + 跟踪。

数据集说明
----------
文件    gps.samples.1bit.I.fs5456.if4092.bin（55.8 MB，1-bit 量化实数采样）
来源    JiaoXianjun/GNSS-GPS-SDR 仓库转存，原始采集见 jks.com/gps/gps.html
参数    采样率 fs = 5.456 MHz，模拟中频 4.092 MHz，1 bit/采样点

关于数字中频（这里有个坑）：
    4.092 MHz > fs/2 = 2.728 MHz，属于欠采样（带通采样），信号会混叠。
    因 4.092 / 5.456 = 0.75 恰好为有理数，混叠后的数字中频为
        f_digital = |4.092 - 5.456| = 1.364 MHz
    用 1.364 MHz 和 4.092 MHz 都能得到相同的相关峰（cos 是偶函数，
    且 4.092/fs = 0.75 ⇒ cos(1.5πn) = cos(0.5πn)），但多普勒符号相反。
    本脚本统一用混叠后的 1.364 MHz，并在结论里标注符号需要取反。

关于比特序：
    1-bit 数据每字节打包 8 个采样点。该数据集是 **LSB 优先**。
    排错的话相关峰会掉到 1/4 左右 —— 所以脚本两种都试，取峰/次峰比高的那种。

运行
----
    uv run python scripts/02_real_data_test.py

数据文件不入库（*.bin 已被 .gitignore 排除）。下载方式见 data/README.md。
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import gnssrx.plotting  # noqa: F401
from gnssrx.acquisition import acquire_all, detect, refine_doppler   # noqa: E402
from gnssrx.io_if import read_1bit_i                                  # noqa: E402
from gnssrx.tracking import TrackingConfig, track_all                 # noqa: E402

DATA = ROOT / "data" / "raw" / "nottingham_gps_l1_1bit.bin"
FS = 5.456e6
F_IF_ALIASED = 1.364e6       # 混叠后的数字中频
PROCESS_MS = 800             # 处理时长
SETTLE_MS = 200              # 统计指标跳过开头的牵引段


def pick_bit_order(samples: int = int(FS * 0.01)) -> bool:
    """用一小段数据判断字节内的比特顺序：峰/次峰比高的那个就是对的。"""
    best = (0.0, True)
    for msb in (True, False):
        d = read_1bit_i(DATA, n_samples=samples, msb_first=msb)
        res = acquire_all(d, FS, F_IF_ALIASED, ms=5,
                          doppler_half_range=6000, doppler_step=500)
        top = max(r["peak_ratio"] for r in res)
        if top > best[0]:
            best = (top, msb)
    return best[1]


def main() -> int:
    if not DATA.exists():
        print(f"找不到数据文件：{DATA}")
        print("下载方式见 data/README.md")
        return 1

    print("=" * 72)
    print("真实 GPS L1 信号验证 · Nottingham 公开数据集")
    print(f"  fs = {FS/1e6:.3f} MHz | 数字中频 = {F_IF_ALIASED/1e6:.3f} MHz | 1-bit 量化")

    print("\n① 判定字节内比特顺序")
    msb_first = pick_bit_order()
    print(f"   采用 {'MSB' if msb_first else 'LSB'} 优先")

    print(f"\n② 读取 {PROCESS_MS} ms 真实数据")
    data = read_1bit_i(DATA, n_samples=int(FS * PROCESS_MS / 1e3), msb_first=msb_first)
    print(f"   {data.size} 采样点，均值 {data.mean():+.3f}")

    print("\n③ 捕获（32 个 PRN × 全部多普勒频点，5 ms 非相干累积）")
    acqs = acquire_all(data, FS, F_IF_ALIASED, ms=5,
                       doppler_half_range=6000, doppler_step=250)
    found = sorted(detect(acqs, threshold=2.5), key=lambda a: -a["peak_ratio"])
    print(f"   检出 {len(found)} 颗：{[a['prn'] for a in found]}")

    print("\n④ 精频估计 + 闭环跟踪")
    for a in found:
        a["doppler_refined_hz"] = refine_doppler(
            data, FS, F_IF_ALIASED, a["prn"],
            code_phase_samples=a["code_phase_samples"],
            coarse_doppler_hz=a["doppler_hz"])
    results = track_all(data, FS, F_IF_ALIASED, found, n_ms=PROCESS_MS,
                        cfg=TrackingConfig())
    tail = slice(SETTLE_MS, None)

    rows = []
    print(f"\n   {'PRN':>4} {'峰/次峰':>7} {'多普勒*':>9} {'C/N0':>7} {'抖动':>7} {'|I_P|保持':>9} {'判定':>6}")
    for res, a in zip(results, found):
        cn0 = res.cn0_dbhz(pdi_s=1e-3, start=SETTLE_MS)
        jit = float(np.std(res.phase_error_deg()[tail]))
        dop = float(np.mean(res.doppler_hz[tail]))
        keep = float(np.mean(np.abs(res.ip[-50:])) / max(np.mean(np.abs(res.ip[:50])), 1e-9))
        locked = cn0 >= 38.0 and jit <= 20.0
        rows.append((res, a["peak_ratio"], dop, cn0, jit, keep, locked))
        print(f"   {res.prn:>4} {a['peak_ratio']:>7.2f} {dop:>+9.1f} {cn0:>7.1f} "
              f"{jit:>6.1f}° {keep*100:>8.0f}% {'锁定' if locked else '弱':>6}")
    print("   * 欠采样导致多普勒符号相反，真实值需取负号")

    n_locked = sum(1 for r in rows if r[6])
    print(f"\n⑤ 判定：{n_locked} / {len(rows)} 颗稳定锁定")
    ok = n_locked >= 4
    print("   ✅ 真实信号上跟踪正常" if ok else "   ❌ 锁定数量不足")

    # ---- 出图 ----
    out = ROOT / "results" / "02_real_data.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), dpi=120)

    ax = axes[0, 0]
    prns = [r[0].prn for r in rows]
    cn0s = [r[3] for r in rows]
    ax.bar([str(p) for p in prns], cn0s,
           color=["#3ddc97" if r[6] else "#8a6d3b" for r in rows])
    ax.axhline(38, ls="--", lw=1, color="#ffb454")
    ax.set_title("各卫星 C/N₀ 估计（真实信号）")
    ax.set_ylabel("dB-Hz"); ax.set_xlabel("PRN"); ax.grid(axis="y", alpha=0.25)

    ax = axes[0, 1]
    for res, *_ in [(r[0],) for r in rows[:5]]:
        ax.plot(np.arange(res.n_blocks), np.abs(res.ip), lw=0.6,
                label=f"PRN {res.prn}")
    ax.set_title("即时支路幅度 |I_P|（前 5 颗）")
    ax.set_xlabel("时间 (ms)"); ax.set_ylabel("|I_P|")
    ax.legend(fontsize=8); ax.grid(alpha=0.25)

    ax = axes[1, 0]
    best = rows[0][0]
    ax.scatter(best.ip[SETTLE_MS:], best.qp[SETTLE_MS:], s=3, alpha=0.3, color="#4da3ff")
    ax.axhline(0, color="#5c6b7a", lw=0.6); ax.axvline(0, color="#5c6b7a", lw=0.6)
    ax.set_title(f"PRN {best.prn} 星座图（真实信号）")
    ax.set_xlabel("I_P"); ax.set_ylabel("Q_P"); ax.grid(alpha=0.25)

    ax = axes[1, 1]
    ax.plot(np.arange(best.n_blocks), best.doppler_hz, lw=0.8, color="#4da3ff")
    ax.set_title(f"PRN {best.prn} 多普勒牵引（欠采样，符号相反）")
    ax.set_xlabel("时间 (ms)"); ax.set_ylabel("多普勒 (Hz)"); ax.grid(alpha=0.25)

    fig.suptitle("真实 GPS L1 信号验证 · Nottingham 数据集 · fs 5.456 MHz · 1-bit",
                 fontsize=13)
    fig.tight_layout()
    fig.savefig(out)
    print(f"\n   图已保存：{out}")
    print("=" * 72)
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
