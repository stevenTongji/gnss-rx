#!/usr/bin/env python3
"""时基自检：贯穿「捕获 → 跟踪 → 位同步 → 帧同步 → t_user」的合成验证。

目的
----
PVT 的伪距由 `t_user = m_blk/1000 + τ/fs` 得到（m_blk = 整数毫秒，τ = 码相位子毫秒）。
本脚本用已知真值的合成信号，检验这个「整数毫秒 + 亚毫秒」口径是否正确。

原理（关键）
-----------
合成 K 颗卫星，让它们携带【完全相同】的导航电文比特流（同一 TOW 序列），
且【没有传播时延】（信号同时刻到达）。各星的码相位、多普勒各自不同。
于是对同一个子帧：
    各星解出的「接收时刻」 t_user = m_blk/1000 + τ/fs 必须完全相同（相差一个系统常数）。
若各星 t_user 不一致（散布 ~ 亚毫秒以上），说明 m_blk 或 τ 的口径有 bug。

这与真实数据的区别：真实数据里各星有传播时延（t_user 之差 = 距离差/c），
所以真实数据无法直接判定「口径」对不对；合成数据（无时延）可以。

运行：
    uv run python scripts/15_timing_selfcheck.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from gnssrx.acquisition import acquire_all, dedupe_detections, detect, \
    refine_doppler                                                        # noqa: E402
from gnssrx.ca_code import CODE_LENGTH, CODE_RATE_HZ, L1_HZ, generate_ca  # noqa: E402
from gnssrx.ephemeris import read_tow, strip_parity                       # noqa: E402
from gnssrx.nav_msg import BITS_PER_SUBFRAME, PREAMBLE, bit_sync, \
    extract_bits, find_subframes                                          # noqa: E402
from gnssrx.tracking import TrackingConfig, track_all                     # noqa: E402

FS = 4.092e6
F_IF = 1.023e6
N_MS = 18000                       # 18 s = 900 比特 = 3 个子帧（够 find_subframes 成链）
PRNS = [7, 11, 19, 24]
CN0 = 45.0
SETTLE_MS = 0
PDI_MS = 1


def build_nav_bits(n_bits: int, tow0_raw: int = 100000, seed: int = 7) -> np.ndarray:
    """连续子帧比特流：前导 + HOW(TOW) 正确设置，其余数据位随机
    （随机位是为了让位同步/帧同步有足够的跳变，能被可靠检出）。"""
    rng = np.random.default_rng(seed)
    bits = np.zeros(n_bits, dtype=np.int8)
    pos, j = 0, 0
    while pos + BITS_PER_SUBFRAME <= n_bits:
        d240 = rng.integers(0, 2, size=240).astype(np.int8)
        d240[0:8] = PREAMBLE
        raw = tow0_raw + j                      # 每个子帧 +1 → TOW +6s
        for k in range(17):
            d240[24 + k] = (raw >> (16 - k)) & 1
        w = np.zeros(300, dtype=np.int8)
        for wd in range(10):                    # 数据位入位，奇偶位=0（本测试不需要校验）
            w[wd * 30:wd * 30 + 24] = d240[wd * 24:(wd + 1) * 24]
        bits[pos:pos + 300] = w
        pos += 300
        j += 1
    return bits


def synthesize(prns, nav_bits, seed=20260929):
    """合成多星 IF 信号。导航电文调制在码上（比特边界严格对齐到码历元，与真实 GPS 一致）。"""
    rng = np.random.default_rng(seed)
    spc = int(round(FS / 1e3))
    n = spc * N_MS
    t = np.arange(n) / FS
    spc_chip = FS / CODE_RATE_HZ

    sig = np.zeros(n)
    truth = []
    for prn in prns:
        code = generate_ca(prn).astype(np.float64)
        delay = float(rng.uniform(0.0, spc))
        dop = float(rng.uniform(-3000.0, 3000.0))
        ph = float(rng.uniform(0.0, 2 * np.pi))
        cr = CODE_RATE_HZ * (1.0 + dop / L1_HZ)
        # 连续码相位（chip）：chip_cont = n*cr/fs − delay/spc_chip
        chip_cont = np.arange(n) * (cr / FS) - delay / spc_chip
        chip = np.floor(chip_cont).astype(np.int64) % CODE_LENGTH
        # 导航比特边界严格对齐到「码历元」（每 20 个码历元一个比特），这才符合真实 GPS
        nav_idx = np.floor(chip_cont / (20.0 * CODE_LENGTH)).astype(np.int64)
        nav_idx = np.clip(nav_idx, 0, len(nav_bits) - 1)
        nav = nav_bits[nav_idx].astype(np.float64) * 2.0 - 1.0
        sig += code[chip] * nav * np.cos(2 * np.pi * (F_IF + dop) * t + ph)
        truth.append({"prn": prn, "delay_samples": delay, "doppler_hz": dop,
                      "code_rate_hz": cr})
    sigma = np.sqrt((0.5 / 10.0 ** (CN0 / 10.0)) * FS / 2.0)
    sig += rng.normal(0.0, sigma, size=n)
    return sig, truth


def true_subframe_time(tr: dict, j: int) -> float:
    """合成信号里，第 j 个子帧起始的【真实到达时刻】（秒）。

    子帧 j 从码历元 20460*300*j 起（每个比特 = 20 个码历元 = 20*1023 chip）。
    chip_cont = n*cr/fs − delay/spc_chip ⇒ n = (target_chip + delay/spc_chip)*fs/cr。
    """
    spc_chip = FS / CODE_RATE_HZ
    target_chip = 20.0 * CODE_LENGTH * BITS_PER_SUBFRAME * j
    n_true = (target_chip + tr["delay_samples"] / spc_chip) * FS / tr["code_rate_hz"]
    return n_true / FS


def unwrap_code_phase(tau, fs, pdi_ms):
    spc = fs / 1.0e3
    spc_chip = fs / 1.023e6
    n = len(tau)
    u = np.zeros(n)
    u[0] = tau[0] / spc_chip
    for m in range(1, n):
        d = float(tau[m] - tau[m - 1])
        d -= round(d / spc) * spc
        u[m] = u[m - 1] + d / spc_chip
    return u


def main() -> int:
    print("=" * 74)
    print("时基自检：合成多星（同电文、无传播时延）→ 各星 t_user 必须一致")
    n_bits = int(N_MS / 20) + 2
    nav = build_nav_bits(n_bits)
    print(f"  合成 {len(PRNS)} 颗星（{PRNS}），{N_MS/1e3:.0f} s，C/N₀ {CN0} dB-Hz，"
          f"电文 {n_bits} 比特（含 {n_bits//300} 个子帧）")
    sig, truth = synthesize(PRNS, nav)

    print("\n① 捕获")
    acqs = acquire_all(sig, FS, F_IF, ms=5, doppler_half_range=6000, doppler_step=250)
    found = sorted(detect(acqs, threshold=2.0), key=lambda a: -a["peak_ratio"])
    found = dedupe_detections(found, FS)
    print(f"   检出 {[a['prn'] for a in found]}")

    print(f"\n② 跟踪 {N_MS} ms")
    for a in found:
        a["doppler_refined_hz"] = refine_doppler(
            sig, FS, F_IF, a["prn"], code_phase_samples=a["code_phase_samples"],
            coarse_doppler_hz=a["doppler_hz"])
    results = track_all(sig, FS, F_IF, found, n_ms=N_MS, cfg=TrackingConfig())

    print("\n③ 位同步 + 帧同步 → t_user = m_blk/1000 + τ/fs")
    rows = []
    for r in results:
        ip = r.ip[SETTLE_MS:]
        offset, conf = bit_sync(ip)
        bits, _ = extract_bits(ip, offset)
        chains = find_subframes(bits)
        if not chains:
            print(f"   PRN {r.prn:>2}  无子帧链")
            continue
        starts = chains[0]["starts"]
        u = unwrap_code_phase(r.code_phase_samples, FS, PDI_MS)
        # 对链上每个子帧算 t_user
        tus = []
        tows = []
        mbs = []
        taus = []
        for s in starts:
            if s + BITS_PER_SUBFRAME > len(bits):
                continue
            m_blk = SETTLE_MS + offset + s * 20 * PDI_MS
            if m_blk >= len(u):
                m_blk = len(u) - 1
            tus.append(m_blk / 1000.0 - u[m_blk] / 1.023e6)
            tows.append(read_tow(strip_parity(bits[s:s + BITS_PER_SUBFRAME])[0]))
            mbs.append(m_blk)
            taus.append(u[m_blk])
        rows.append((r.prn, offset, tus, tows, mbs, taus, r.cn0_dbhz(1e-3, 0)))
        print(f"   PRN {r.prn:>2}  offset={offset}ms conf={conf:.2f}  C/N₀={r.cn0_dbhz(1e-3,0):.1f}  "
              f"子帧 TOW={tows}")

    print("\n④ 判定：把推得的 t_user 与合成里已知的真实子帧到达时刻对比")
    print("   （delta = t_user − t_true 必须是公共常数；同时打印 m_blk / τ / n_true 便于定位）")
    spc = int(round(FS / 1e3))
    tr_by_prn = {t["prn"]: t for t in truth}
    deltas = []
    for prn, off, tus, tows, mbs, taus, _ in rows:
        tr = tr_by_prn[prn]
        for tuser, tow, mb, tau in zip(tus, tows, mbs, taus):
            j = (tow - 600000) // 6
            n_true = true_subframe_time(tr, j) * FS
            deltas.append((prn, j, tuser - n_true / FS))
            frac_true = (n_true - mb * spc) / spc
            print(f"   PRN {prn:>2} 子帧{j}: m_blk={mb:6d} τ={tau:8.2f}chip({tau/spc:+.4f}ms) "
                  f"n_true={n_true:10.1f}({frac_true:+.4f}ms)  delta={(tuser-n_true/FS)*1e3:+.4f}ms")
    ds = np.array([d for _, _, d in deltas])
    spread_ms = (ds.max() - ds.min()) * 1e3
    ok = spread_ms < 0.03
    print(f"\n   delta 散布 = {spread_ms:.4f} ms  {'✅ 口径正确' if ok else '❌ 口径有 bug'}")

    print("\n" + "=" * 74)
    print("✅ 时基口径正确：t_user 精确恢复真实到达时刻（仅差公共常数）。" if ok
          else "❌ 时基口径有 bug：t_user 与真实到达时刻之差随卫星/子帧变化。")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
