#!/usr/bin/env python3
"""审计 SiGe 数据里所有被跟踪通道的星历解码结果。

目的：搞清楚 PRN 17/8/28/1 等"失败"通道到底是
  (a) C/A 互相关重影（解码出的轨道其实是另一颗真星，IODE 等也该一致 → 但可能
      因弱信号而错位）
  (b) 真实弱星但个别比特出错导致 IODE 自洽校验失败（可通过放宽校验 / 纠错恢复）
  (c) 纯噪声

直接复用 sige_ephemeris_blocks.pkl（已由 05 用 28 通道、门限 1.0 生成），
无需重新跑耗时的跟踪。
"""
from __future__ import annotations

import pickle
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from gnssrx.ephemeris import decode_subframes, strip_parity, read_tow, sfid_from_tow  # noqa: E402

CACHE = ROOT / "data" / "processed" / "sige_ephemeris_blocks.pkl"


def main() -> int:
    if not CACHE.exists():
        print(f"找不到缓存 {CACHE}，请先跑 scripts/05_sige_ephemeris_test.py")
        return 1
    with open(CACHE, "rb") as f:
        blocks_by_prn = pickle.load(f)

    print(f"缓存里共 {len(blocks_by_prn)} 个通道：{sorted(blocks_by_prn)}")

    print("\n=== 逐通道解码审计 ===")
    decoded = {}
    for prn in sorted(blocks_by_prn):
        blocks = blocks_by_prn[prn]
        eph, info = decode_subframes(prn, blocks)
        if eph is None:
            # 哪些子帧缺失
            sfids = sorted({p["sfid"] for p in info})
            print(f"  PRN {prn:>2}  三类子帧不全  已捕获子帧号={sfids}  "
                  f"块数={len(blocks)}")
            continue
        probs = eph.self_check()
        decoded[prn] = eph
        status = "✅ 通过" if not probs else "❌ 失败"
        print(f"  PRN {prn:>2}  {status}  sqrtA={eph.sqrtA:.1f} e={eph.e:.4f} "
              f"i0={eph.i0*180/np.pi:.2f}°  IODE2={eph.iode2} IODE3={eph.iode3} "
              f"IODC={eph.iodc}")
        if probs:
            print(f"          问题: {probs}")

    clean = {p: e for p, e in decoded.items() if not e.self_check()}
    ghosts = {p: e for p, e in decoded.items() if e.self_check()}
    print(f"\n通过自检: {sorted(clean)}")
    print(f"失败但解出三类子帧: {sorted(ghosts)}")

    # 用轨道身份指纹判断失败通道是否其实是某颗真星的重影
    print("\n=== 轨道身份指纹对比（检测重影 / 真实弱星）===")
    print("身份 = (round(sqrtA,1), iode2, iode3)")

    def keyof(e):
        return (round(e.sqrtA, 1), e.iode2, e.iode3)

    clean_keys = {keyof(e): p for p, e in clean.items()}
    print("  干净星指纹:")
    for p, e in sorted(clean.items()):
        print(f"    PRN {p:>2}: {keyof(e)}  A={e.A/1e3:.1f}km e={e.e:.5f}")

    for p, e in sorted(ghosts.items()):
        k = keyof(e)
        match = clean_keys.get(k)
        verdict = f"↔ 与 PRN {match} 指纹相同 = 重影/同一星" if match else "独立轨道 = 可能是真实弱星"
        print(f"    PRN {p:>2}: {k}  A={e.A/1e3:.1f}km e={e.e:.5f}  → {verdict}")

    # 对"失败但物理合理"的通道尝试放宽校验（仅检查轨道物理性，不查 IODE 自洽）
    print("\n=== 放宽校验：仅看轨道参数是否物理合理（疑似真实弱星）===")
    for p, e in sorted(ghosts.items()):
        phys_ok = (4700.0 < e.sqrtA < 5600.0) and (0 < e.e <= 0.2) and \
                  (abs(abs(e.i0) - 0.96) < 0.15) and (e.A > 0)
        # 是否与其他真星指纹碰撞
        is_ghost = keyof(e) in clean_keys
        tag = "物理合理" if phys_ok else "物理异常"
        tag += " / 重影" if is_ghost else " / 独立"
        print(f"    PRN {p:>2}: {tag}  sqrtA={e.sqrtA:.1f} e={e.e:.4f} "
              f"i0={e.i0*180/np.pi:.2f}°")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
