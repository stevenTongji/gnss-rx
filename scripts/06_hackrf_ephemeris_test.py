#!/usr/bin/env python3
"""星历解码验证（HackRF 真实数据，8-bit I/Q → 实数 IF）。

数据集：ION GNSS SDR 元数据标准的 HackRF 样本
    fs = 10 MHz, IF = 420 kHz, 8-bit 双补码 I/Q 交替（每样本 2 字节）
    收录卫星（ION 标注）：PRN 2,5,6,7,9,13,29,30 共 8 颗（信号强）

本脚本先把 I/Q 文件转换成「实数 IF」int8（复用现有实数链路），再做
捕获 → 跟踪 → 帧同步 → 星历解码。HackRF 卫星多且强，目标是让 ≥4 颗
解出完整 SF1/2/3，为 PVT 铺路。

运行：
    uv run python scripts/06_hackrf_ephemeris_test.py          # 复用缓存
    uv run python scripts/06_hackrf_ephemeris_test.py --regen  # 重新转换/捕获/跟踪
"""
from __future__ import annotations

import pickle
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from gnssrx.acquisition import acquire_all, dedupe_detections, detect, \
    refine_doppler                                                        # noqa: E402
from gnssrx.ephemeris import decode_subframes, sfid_from_tow               # noqa: E402
from gnssrx.io_if import iq_to_real_int8, read_real_int8                  # noqa: E402
from gnssrx.nav_msg import BITS_PER_SUBFRAME, bit_sync, extract_bits, \
    find_subframes                                                       # noqa: E402
from gnssrx.tracking import TrackingConfig, track_all                     # noqa: E402

DATA_IQ = ROOT / "data" / "raw" / "hackrf_gps_l1_8bit_iq.int8"
DATA_REAL = ROOT / "data" / "raw" / "hackrf_gps_l1_8bit_real.dat"
FS = 10.0e6
F_IF = 420.0e3
PROCESS_MS = 30000          # 30 秒：覆盖完整 5 子帧主帧，足够 PVT
N_CHANNELS = 14             # HackRF 强信号多，14 通道即可，避免跟踪过慢
THRESHOLD = 1.3             # 门限适中，重影靠 self_check 剔除
SETTLE_MS = 500
CACHE = ROOT / "data" / "processed" / "hackrf_ephemeris_blocks.pkl"


def _ensure_real_if() -> None:
    """把 I/Q 文件转换成实数 IF int8（带磁盘缓存）。"""
    if DATA_REAL.exists():
        print(f"   实数 IF 已存在：{DATA_REAL.name}，跳过转换")
        return
    print("   转换 I/Q → 实数 IF（int8，fs=10MHz, IF=420kHz）…")
    n = iq_to_real_int8(DATA_IQ, DATA_REAL, FS, F_IF)
    print(f"   转换完成：{n} 个实数 IF 采样点 → {DATA_REAL.name}")


def get_blocks(regen: bool = False) -> dict[int, list[np.ndarray]]:
    """捕获 + 跟踪 + 帧同步，返回 {prn: [300 比特子帧块(±1) ...]}。带磁盘缓存。"""
    if CACHE.exists() and not regen:
        with open(CACHE, "rb") as f:
            return pickle.load(f)                                     # type: ignore

    _ensure_real_if()
    DATA = DATA_REAL

    print("① 捕获全部 32 颗（门限 %g），取最强 %d 颗" % (THRESHOLD, N_CHANNELS))
    data = read_real_int8(DATA, n_samples=int(FS * PROCESS_MS / 1e3))
    acqs = acquire_all(data, FS, F_IF, ms=5, doppler_half_range=6000, doppler_step=250)
    found = sorted(detect(acqs, threshold=THRESHOLD), key=lambda a: -a["peak_ratio"])
    print("   候选 PRN 与 peak_ratio:")
    for a in found[:N_CHANNELS]:
        print(f"     PRN {a['prn']:>2}  peak_ratio={a['peak_ratio']:.2f}  "
              f"doppler={a['doppler_hz']:+.0f} Hz  code_phase={a['code_phase_samples']}")
    found = dedupe_detections(found, FS)                 # 剔除同相位的 C/A 互相关重影
    found = found[:N_CHANNELS]
    print(f"   去重后 PRN: {[a['prn'] for a in found]}")
    for a in found:
        a["doppler_refined_hz"] = refine_doppler(
            data, FS, F_IF, a["prn"],
            code_phase_samples=a["code_phase_samples"],
            coarse_doppler_hz=a["doppler_hz"])

    print(f"② 闭环跟踪 {PROCESS_MS/1e3:.0f} 秒")
    results = track_all(data, FS, F_IF, found, n_ms=PROCESS_MS, cfg=TrackingConfig())
    for r in results:
        print(f"   PRN {r.prn:>2}  C/N₀ {r.cn0_dbhz(1e-3, SETTLE_MS):>5.1f} dB-Hz")

    print("③ 位同步 + 帧同步，收集子帧块")
    blocks_by_prn: dict[int, list[np.ndarray]] = {}
    cn0_by_prn: dict[int, float] = {}
    for r in results:
        ip = r.ip[SETTLE_MS:]
        offset, _ = bit_sync(ip)
        bits, _ = extract_bits(ip, offset)
        chains = find_subframes(bits)
        if not chains:
            print(f"   PRN {r.prn:>2}  无子帧链，跳过")
            continue
        chain = chains[0]
        blks = []
        for start in chain["starts"]:
            if start + BITS_PER_SUBFRAME > len(bits):
                continue
            blks.append(bits[start:start + BITS_PER_SUBFRAME].astype(np.int8))
        blocks_by_prn[r.prn] = blks
        cn0_by_prn[r.prn] = r.cn0_dbhz(1e-3, SETTLE_MS)
        print(f"   PRN {r.prn:>2}  收集到 {len(blks)} 个子帧块")

    CACHE.parent.mkdir(parents=True, exist_ok=True)
    with open(CACHE, "wb") as f:
        pickle.dump(blocks_by_prn, f)
    print(f"   已缓存到 {CACHE}")
    return blocks_by_prn


def main() -> int:
    if not DATA_IQ.exists():
        print(f"找不到数据文件：{DATA_IQ}\n下载方式见 data/README.md")
        return 1

    regen = "--regen" in sys.argv
    blocks_by_prn = get_blocks(regen=regen)

    print("\n④ 星历解码（全部跟踪通道）")
    all_decoded = {}
    for prn, blocks in blocks_by_prn.items():
        eph, info = decode_subframes(prn, blocks)
        if eph:
            all_decoded[prn] = (eph, info)
            sfids = sorted({s["sfid"] for s in info})
            print(f"   PRN {prn:>2}  解出子帧 {sfids}  TLM同步头全对="
                  f"{all(s['preamble_ok'] for s in info)}")

    clean = {}
    ghosts = {}
    for prn, (eph, info) in all_decoded.items():
        if not eph.self_check():
            clean[prn] = eph
        else:
            ghosts[prn] = eph

    seen = {}
    final = {}
    for prn, eph in clean.items():
        key = (round(eph.sqrtA, 1), eph.iode2, eph.iode3)
        if key in seen:
            print(f"   （PRN {prn} 与 PRN {seen[key]} 轨道身份相同，合并）")
            continue
        seen[key] = prn
        final[prn] = eph

    print(f"\n   解码出 {len(all_decoded)} 颗含三类子帧；"
          f"通过自检 {len(clean)} 颗；去重后 {len(final)} 颗真星。")

    print("\n⑤ 校验 + 星历摘要（通过自检的真星）")
    for prn, eph in final.items():
        rad = eph.radius_km(tk=0.0)
        rad_trans = eph.radius_km(tk=(eph._raw["sf1"]["toc"] - eph.toes))
        print(f"\n✅ PRN {prn:>2}")
        print(f"   周内时 WN={eph.week}  IODE2={eph.iode2}  IODE3={eph.iode3}  "
              f"IODC={eph.iodc}")
        print(f"   半长轴 A={eph.A/1e3:.1f} km  e={eph.e:.5f}  "
              f"i0={eph.i0*180/np.pi:.2f}°  sqrtA={eph.sqrtA:.4f}")
        print(f"   时钟: af0={eph.af0*1e6:+.3f} µs  af1={eph.af1*1e6:+.3f} µs/s"
              f"  Tgd={eph.tgd*1e9:+.2f} ns  toc={eph.toc:.0f}s")
        print(f"   卫星位置半径: @toe={rad:.1f} km   @toc={rad_trans:.1f} km"
              f"   （GPS 标称 ≈ 26560 km）")

    if ghosts:
        print("\n   被自检剔除的可疑通道（C/A 互相关重影 / 弱信号比特错误）：")
        for prn, eph in ghosts.items():
            probs = eph.self_check()
            print(f"     PRN {prn:>2}  sqrtA={eph.sqrtA:.1f} e={eph.e:.3f}  "
                  f"问题: {'; '.join(p.split('：')[0] for p in probs)}")

    print("\n" + "=" * 78)
    n = len(final)
    if n >= 4:
        print(f"✅ {n} 颗卫星星历自检通过 —— 满足 PVT 最低 4 星要求。")
        rc = 0
    elif n > 0:
        print(f"⚠️  仅 {n} 颗卫星星历自检通过，PVT 还差 {4 - n} 颗；"
              f"可再放低门限或换更长数据。")
        rc = 0
    else:
        print("❌ 没有卫星的星历自检通过，请检查链路。")
        rc = 1
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
