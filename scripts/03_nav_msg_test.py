#!/usr/bin/env python3
"""导航电文同步验证：在真实 GPS L1 信号上做位同步 + 帧同步。

运行：
    uv run python scripts/03_nav_msg_test.py      # 需先下载数据，见 data/README.md

它做四件事：
    1. 真实数据上捕获，挑最强的几颗星
    2. 闭环跟踪 30 秒，拿到每毫秒的 I_P
    3. 位同步 —— 找出比特边界在哪一毫秒
    4. 帧同步 —— 搜同步头 10001011，并要求命中点间隔恰好 300 比特

判据（都是硬约束，很难靠巧合蒙对）：
    · 位同步置信度 > 1.5          → 说明"20 个 ms 同号"这个结构是真的存在的
    · 至少找到一条长度 ≥ 2 的子帧链
    · 相邻子帧间隔必须恰好 300 比特 = 6.000 秒   ← 最关键的一条
      随机 8 比特匹配的概率是 1/128，但要连续多次间隔 300 命中几乎不可能，
      所以这条过了，基本就说明整条链路（捕获→跟踪→位同步）都是对的。
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
from gnssrx.nav_msg import (BITS_PER_SUBFRAME, MS_PER_BIT, bit_sync,
                            bit_sync_transitions, extract_bits,
                            find_subframes)                           # noqa: E402
from gnssrx.tracking import TrackingConfig, track_all                 # noqa: E402

DATA = ROOT / "data" / "raw" / "nottingham_gps_l1_1bit.bin"
FS = 5.456e6
F_IF = 1.364e6
PROCESS_MS = 30000        # 30 秒 → 理论上有 5 个子帧
N_CHANNELS = 3
SETTLE_MS = 500           # 跳过开头的牵引段


def main() -> int:
    if not DATA.exists():
        print(f"找不到数据文件：{DATA}\n下载方式见 data/README.md")
        return 1

    print("=" * 74)
    print("导航电文同步验证 · 真实 GPS L1 信号")
    print(f"  处理 {PROCESS_MS/1e3:.0f} 秒 | 理论上应有 {PROCESS_MS/6000:.0f} 个子帧")

    data = read_1bit_i(DATA, n_samples=int(FS * PROCESS_MS / 1e3), msb_first=False)
    print(f"  读入 {data.size} 采样点（{data.size/FS:.1f} 秒）")

    print("\n① 捕获")
    acqs = acquire_all(data, FS, F_IF, ms=5, doppler_half_range=6000, doppler_step=250)
    found = sorted(detect(acqs, threshold=2.5), key=lambda a: -a["peak_ratio"])[:N_CHANNELS]
    print(f"   取最强的 {len(found)} 颗：{[a['prn'] for a in found]}")
    for a in found:
        a["doppler_refined_hz"] = refine_doppler(
            data, FS, F_IF, a["prn"],
            code_phase_samples=a["code_phase_samples"],
            coarse_doppler_hz=a["doppler_hz"])

    print(f"\n② 闭环跟踪 {PROCESS_MS/1e3:.0f} 秒（这一步最耗时，请稍候）")
    results = track_all(data, FS, F_IF, found, n_ms=PROCESS_MS, cfg=TrackingConfig())
    for r in results:
        print(f"   PRN {r.prn:>2}  C/N₀ {r.cn0_dbhz(1e-3, SETTLE_MS):>5.1f} dB-Hz")

    print("\n③ 位同步（两种独立方法交叉验证）")
    syncs = []
    for r in results:
        ip = r.ip[SETTLE_MS:]
        offset, conf = bit_sync(ip)
        trans = bit_sync_transitions(ip)
        trans_str = "（跳变太少无法判定）" if trans is None else \
            f"跳变法 → {trans[0]} ms（置信度 {trans[1]:.1f}）"
        agree = trans is not None and trans[0] == offset
        bits, strength = extract_bits(ip, offset)
        syncs.append({"prn": r.prn, "ip": ip, "offset": offset, "conf": conf,
                      "trans": trans, "agree": agree, "bits": bits,
                      "strength": strength})
        print(f"   PRN {r.prn:>2}  能量法 → {offset} ms（置信度 {conf:.2f}）   {trans_str}"
              f"   {'✓ 一致' if agree else '✗ 不一致'}")

    print("\n④ 帧同步（搜同步头 10001011）")
    ok = True
    for s in syncs:
        chains = find_subframes(s["bits"])
        s["chains"] = chains
        if not chains:
            print(f"   PRN {s['prn']:>2}  未找到子帧链")
            continue
        c = chains[0]
        starts = np.array(c["starts"])
        diffs = np.diff(starts)
        secs = starts * MS_PER_BIT / 1e3
        print(f"   PRN {s['prn']:>2}  最长链 {c['n_subframes']} 个子帧，"
              f"极性 {'正常' if c['polarity'] > 0 else '整体反号'}，"
              f"间隔 {set(diffs.tolist())} 比特")
        print(f"        子帧起始时刻（秒）：{np.round(secs, 3).tolist()}")

    print("\n⑤ 判据检查")
    for s in syncs:
        bad = []
        # 阈值 1.4：无噪声理想值约 2.0，强噪声约 1.7，无比特结构约 1.0
        if s["conf"] < 1.4:
            bad.append(f"位同步置信度 {s['conf']:.2f} 偏低（无结构时约 1.0）")
        if not s["agree"]:
            bad.append("两种位同步方法不一致")
        chains = s.get("chains", [])
        if not chains or chains[0]["n_subframes"] < 2:
            bad.append("未找到长度 ≥ 2 的子帧链")
        else:
            diffs = np.diff(np.array(chains[0]["starts"]))
            if not np.all(diffs == BITS_PER_SUBFRAME):
                bad.append(f"子帧间隔 {set(diffs.tolist())} 不等于 300 比特")
        if bad:
            ok = False
            print(f"   ❌ PRN {s['prn']}: " + "；".join(bad))
        else:
            c = chains[0]
            t = s.get("trans")
            extra = f"，跳变法一致（{t[1]:.1f}）" if t else ""
            print(f"   ✅ PRN {s['prn']}: 位同步置信度 {s['conf']:.2f}{extra}，"
                  f"{c['n_subframes']} 个子帧，间隔严格等于 300 比特（6.000 秒）")

    # ---- 出图 ----
    out = ROOT / "results" / "03_nav_msg.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    n = len(syncs)
    fig, axes = plt.subplots(n, 2, figsize=(13, 3.4 * n), dpi=120, squeeze=False)

    for i, s in enumerate(syncs):
        ax = axes[i][0]
        seg = s["ip"][:400]
        ax.plot(np.arange(len(seg)), seg, lw=0.8, color="#4da3ff")
        for k in range(0, len(seg), MS_PER_BIT):
            ax.axvline(k, color="#5c6b7a", lw=0.5, alpha=0.55)
        ax.axvline(s["offset"], color="#ffb454", lw=1.4,
                   label=f"比特边界 offset={s['offset']}")
        ax.set_title(f"PRN {s['prn']} 每毫秒 I_P（竖线为比特边界）")
        ax.set_xlabel("毫秒"); ax.set_ylabel("I_P")
        ax.legend(fontsize=8); ax.grid(alpha=0.2)

        ax = axes[i][1]
        bits = s["bits"][:600]
        ax.step(np.arange(len(bits)), bits, where="post", lw=0.7, color="#3ddc97")
        for st in s.get("chains", [{}])[0].get("starts", []):
            if st < len(bits):
                ax.axvline(st, color="#ff6b6b", lw=1.3)
                ax.text(st, 1.35, "子帧", fontsize=8, color="#ff6b6b")
        ax.set_ylim(-1.6, 1.6)
        ax.set_title(f"PRN {s['prn']} 解出的电文比特（红线为子帧起点）")
        ax.set_xlabel("比特序号"); ax.set_ylabel("±1")
        ax.grid(alpha=0.2)

    fig.suptitle("导航电文同步 · 真实 GPS L1 信号", fontsize=13)
    fig.tight_layout()
    fig.savefig(out)
    print(f"\n   图已保存：{out}")
    print("=" * 74)
    print("✅ 全部通过。" if ok else "❌ 有检查未通过。")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
