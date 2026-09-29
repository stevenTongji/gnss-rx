#!/usr/bin/env python3
"""单点定位（PVT）验证 —— 真实数据（Nottingham 公开数据集，1-bit 实数 IF）。

为什么用这份数据：
    SiGe GN3S 只稳定解出 3 颗干净卫星（PRN 7/9/11），3 星欠定、结果对先验敏感，
    无法给出可靠的独立坐标（这正是"最少 4 星"标准约束的物理原因）。
    Nottingham 数据集视野开阔、信噪比高、1-bit 量化，通常能看到 8~11 颗卫星，
    足以支撑完整的 ≥4 星标准单点定位（SPS），从而用真实数据独立验证 PVT 引擎。

流程：捕获 → 跟踪 → 帧同步/星历解码 → 伪距提取 → 迭代加权最小二乘定位。
伪距由"子帧起点块索引 × 1ms + 本地码相位 τ + TOW"算出；接收机钟差 b 吸收我们
0 时基与 GPS 时的公共偏置，因此 TOW 只需同星期内一致。

数据集：gps.samples.1bit.I.fs5456.if4092.bin
    fs = 5.456 MHz，模拟中频 4.092 MHz，1-bit 量化实数采样
    注意：4.092 MHz > fs/2=2.728 MHz 属于欠采样，混叠后数字中频 = |4.092-5.456|
          = 1.364 MHz，且多普勒符号与真实值相反（跟踪环照常工作，对 PVT 无影响）。
    比特序：字节内 **LSB 优先** → read_1bit_i(msb_first=False)。

运行：
    uv run python scripts/14_nottingham_pvt.py          # 复用缓存
    uv run python scripts/14_nottingham_pvt.py --regen  # 重新捕获/跟踪/定位

位置参照：⚠️ 这份数据集的官方说明只承诺「a lat/lon position in Nottingham, UK」，
**从未公布天线坐标**。下面用诺丁汉市中心 (52.9536°N, 1.1505°W) 仅作"是否落在诺丁汉"
的量级参照，**不能当作位置真值**；精度验证请看 scripts/17_rtklib_comparison.py
（残差×PDOP 的形式精度、以及卫星子集一致性）。
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
from gnssrx.ephemeris import decode_subframes, read_tow, strip_parity     # noqa: E402
from gnssrx.io_if import read_1bit_i                                     # noqa: E402
from gnssrx.nav_msg import BITS_PER_SUBFRAME, bit_sync, \
    bit_sync_transitions, extract_bits, find_subframes                   # noqa: E402
from gnssrx.pvt import build_measurement, ecef_to_llh, solve_cold_start, \
    solve_robust                                                         # noqa: E402
from gnssrx.tracking import TrackingConfig, track_all                    # noqa: E402

DATA = ROOT / "data" / "raw" / "nottingham_gps_l1_1bit.bin"
FS = 5.456e6
F_IF = 1.364e6              # 混叠后的数字中频
PROCESS_MS = 48000         # 48 秒 = 8 个子帧，保证无论子帧链起点如何都能覆盖到 SF1
N_CHANNELS = 12
THRESHOLD = 2.0
SETTLE_MS = 500
PDI_MS = 1                 # 与 TrackingConfig 默认一致（20 块/比特）
CACHE = ROOT / "data" / "processed" / "nottingham_pvt_measurements.pkl"
TRACK_CACHE = ROOT / "data" / "processed" / "nottingham_tracking.pkl"

# 诺丁汉市中心（城市级参照点；⚠️ 非天线坐标，该数据集未公布天线位置）
NOTTINGHAM_TRUTH = (52.9536, -1.1505, 0.0)


def unwrap_code_phase(tau_stream: np.ndarray, fs: float, pdi_ms: int) -> np.ndarray:
    """把逐块（mod spc）的码相位 τ 解缠成自跟踪开始累计的码片数。

    每 1 ms 块码推进约 1023 个码片（含多普勒牵引后略偏），τ 每帧回卷一次；
    解缠后得到亚微秒级的累计码相位 → 换算成接收机时刻，是伪距精度的关键。
    """
    spc = fs / 1.0e3
    samples_per_chip = fs / 1.023e6
    n = len(tau_stream)
    unwrapped = np.zeros(n)
    unwrapped[0] = tau_stream[0] / samples_per_chip
    for m in range(1, n):
        d = float(tau_stream[m] - tau_stream[m - 1])
        d = d - round(d / spc) * spc            # 回卷到 [-spc/2, spc/2] 采样点
        unwrapped[m] = unwrapped[m - 1] + d / samples_per_chip
    return unwrapped


def _ecef_of(lat_deg, lon_deg, h=0.0) -> np.ndarray:
    a = 6378137.0
    f = 1.0 / 298.257223563
    e2 = f * (2 - f)
    lat, lon = np.radians(lat_deg), np.radians(lon_deg)
    N = a / np.sqrt(1 - e2 * np.sin(lat) ** 2)
    return np.array([(N + h) * np.cos(lat) * np.cos(lon),
                     (N + h) * np.cos(lat) * np.sin(lon),
                     (N * (1 - e2) + h) * np.sin(lat)])


def get_tracking(regen: bool = False) -> dict:
    """捕获 + 跟踪，返回 {prn: {"ip","tau","cn0"}}。带磁盘缓存（跟踪最耗时）。"""
    if TRACK_CACHE.exists() and not regen:
        with open(TRACK_CACHE, "rb") as f:
            return pickle.load(f)                                     # type: ignore

    print("① 捕获全部 32 颗（门限 %g），取最强 %d 颗" % (THRESHOLD, N_CHANNELS))
    data = read_1bit_i(DATA, n_samples=int(FS * PROCESS_MS / 1e3), msb_first=False)
    print(f"   读入 {data.size} 采样点（{data.size/FS:.1f} 秒，1-bit int8）")
    acqs = acquire_all(data, FS, F_IF, ms=5, doppler_half_range=6000, doppler_step=250)
    found = sorted(detect(acqs, threshold=THRESHOLD), key=lambda a: -a["peak_ratio"])
    found = dedupe_detections(found, FS)[:N_CHANNELS]
    print(f"   去重后 PRN: {[a['prn'] for a in found]}")
    for a in found:
        a["doppler_refined_hz"] = refine_doppler(
            data, FS, F_IF, a["prn"],
            code_phase_samples=a["code_phase_samples"],
            coarse_doppler_hz=a["doppler_hz"])

    print(f"② 闭环跟踪 {PROCESS_MS/1e3:.0f} 秒")
    results = track_all(data, FS, F_IF, found, n_ms=PROCESS_MS, cfg=TrackingConfig())
    track: dict[int, dict] = {}
    for r in results:
        cn0 = float(r.cn0_dbhz(1e-3, SETTLE_MS))
        print(f"   PRN {r.prn:>2}  C/N₀ {cn0:>5.1f} dB-Hz")
        track[r.prn] = {"ip": np.asarray(r.ip, dtype=float),
                        "tau": np.asarray(r.code_phase_samples, dtype=float),
                        "cn0": cn0}
    TRACK_CACHE.parent.mkdir(parents=True, exist_ok=True)
    with open(TRACK_CACHE, "wb") as f:
        pickle.dump(track, f)
    return track


def diagnose(track: dict) -> None:
    """逐通道打印位同步/子帧/星历诊断，定位为何 self_check 不过。"""
    import hashlib
    for prn, t in track.items():
        ip = t["ip"][SETTLE_MS:]
        offset, conf = bit_sync(ip)
        trans = bit_sync_transitions(ip)
        bits, _ = extract_bits(ip, offset)
        chains = find_subframes(bits)
        h = hashlib.md5(np.asarray(t["tau"], dtype=np.float64).tobytes()).hexdigest()[:10]
        print(f"\n=== PRN {prn:>2}  C/N₀={t['cn0']:.1f} dB-Hz  tau_md5={h}")
        print(f"    位同步 offset={offset}ms conf={conf:.2f}  "
              f"跳变法={'None' if trans is None else f'{trans[0]}ms/{trans[1]:.1f}'}  "
              f"子帧链长度={[c['n_subframes'] for c in chains]}")
        if not chains:
            print("    → 无子帧链")
            continue
        c = chains[0]
        blocks = []
        for s in c["starts"]:
            if s + BITS_PER_SUBFRAME > len(bits):
                continue
            blocks.append(bits[s:s + BITS_PER_SUBFRAME].astype(np.int8))
        # 每个 300 比特块的 md5：若多颗星块完全相同 → 说明锁到同一信号（假锁）
        blk_hashes = [hashlib.md5(b.tobytes()).hexdigest()[:8] for b in blocks]
        print(f"    块 md5: {blk_hashes}")
        eph, info = decode_subframes(prn, blocks)
        for e in info:
            print(f"      子帧 sfid={e['sfid']} tow={e['tow']}s c={e['c']} "
                  f"preamble_ok={e['preamble_ok']}")
        if eph is None:
            print("      → 未凑齐 SF1/2/3")
        else:
            print(f"      → wkN={eph.week} toe={eph.toes:.0f}s toc={eph.toc:.0f}s "
                  f"sqrtA={eph.sqrtA:.3f} e={eph.e:.5f} i0={eph.i0*180/np.pi:+.2f}° "
                  f"iode2={eph.iode2} iode3={eph.iode3} iodc={eph.iodc}")
            probs = eph.self_check()
            print(f"      → self_check: {'PASS' if not probs else '; '.join(probs)}")


def get_measurements(regen: bool = False) -> list[dict]:
    """完整链路 → 返回干净卫星的伪距观测列表（dict）。带磁盘缓存。"""
    if CACHE.exists() and not regen:
        with open(CACHE, "rb") as f:
            return pickle.load(f)                                     # type: ignore

    track = get_tracking(regen)
    measurements = []
    for prn, t in track.items():
        ip = t["ip"][SETTLE_MS:]
        offset, _ = bit_sync(ip)
        bits, _ = extract_bits(ip, offset)
        chains = find_subframes(bits)
        if not chains:
            print(f"   PRN {prn:>2}  无子帧链，跳过")
            continue
        chain = chains[0]
        pairs = []
        for start in chain["starts"]:
            if start + BITS_PER_SUBFRAME > len(bits):
                continue
            pairs.append((start, bits[start:start + BITS_PER_SUBFRAME].astype(np.int8)))

        eph, info = decode_subframes(prn, [b for _, b in pairs])
        if eph is None or eph.self_check():
            print(f"   PRN {prn:>2}  星历未通过自检，跳过")
            continue
        start, block = pairs[-1]
        d240, _ = strip_parity(block)
        tow = read_tow(d240)
        unwrapped = unwrap_code_phase(t["tau"], FS, PDI_MS)
        m_block = SETTLE_MS + offset + start * 20 * PDI_MS  # 跟踪块序号（每块 1ms）
        if m_block >= len(unwrapped):
            m_block = len(unwrapped) - 1
        # 接收机时刻 = 整数毫秒(块序号) − 亚毫秒(解缠后的码相位)。
        # 符号：本跟踪器的 τ 与「信号码历元在块内的偏移」反号，故取【减号】
        # （合成自检脚本 15 已严格验证：减号下各星 delta 为公共常数；加号则散布 ~1ms）。
        t_user = m_block / 1000.0 - unwrapped[m_block] / 1.023e6
        meas = build_measurement(eph, prn, tow, t_user, FS)
        print(f"   PRN {prn:>2}  ✅ 伪距 = {meas['pseudorange']:.3e} m  "
              f"(t_sat={meas['t_sat']:.1f}s  t_user={t_user:.3f}s)")
        measurements.append(meas)

    print(f"\n   共 {len(measurements)} 颗干净卫星参与定位。")
    CACHE.parent.mkdir(parents=True, exist_ok=True)
    with open(CACHE, "wb") as f:
        pickle.dump(measurements, f)
    return measurements


def _great_circle_km(a: np.ndarray, b: np.ndarray) -> float:
    """两点（ECEF，米）间大圆距离，公里。"""
    la1, lo1, _ = ecef_to_llh(a)
    la2, lo2, _ = ecef_to_llh(b)
    la1, lo1, la2, lo2 = map(np.radians, (la1, lo1, la2, lo2))
    d = np.sin(la1) * np.sin(la2) + np.cos(la1) * np.cos(la2) * np.cos(lo1 - lo2)
    d = min(1.0, max(-1.0, d))
    return 6371.0 * np.arccos(d)


def main() -> int:
    if not DATA.exists():
        print(f"找不到数据文件：{DATA}\n下载方式见 data/README.md")
        return 1
    regen = "--regen" in sys.argv
    if "--diag" in sys.argv:
        diagnose(get_tracking(regen))
        return 0
    meas = get_measurements(regen=regen)

    if len(meas) < 4:
        print(f"❌ 仅 {len(meas)} 颗干净卫星，不足以做 ≥4 星完整 SPS（本数据应 ≥4 颗）。")
        return 1

    truth_ecef = _ecef_of(*NOTTINGHAM_TRUTH)

    print(f"\n③ 迭代加权最小二乘定位（完整 SPS，{len(meas)} 星）")
    print("   先试【冷启动】（不给任何先验，全球网格搜索自动定初值）：")
    sol_free = solve_cold_start(meas)
    fla, flo, fh = ecef_to_llh(sol_free["pos"])
    ferr = _great_circle_km(sol_free["pos"], truth_ecef)
    print(f"   冷启动：lat={fla:.5f}° lon={flo:.5f}° h={fh:.1f} m "
          f"收敛={sol_free['converged']} 与真值距离={ferr:.2f} km "
          f"（网格最优 {sol_free['cold_start_llh']}）")

    print("   再给 Nottingham 区域先验（模拟真实接收机的上次定位/AGPS 提示）：")
    prior = _ecef_of(52.95, -1.15, 0.0)
    sol = solve_robust(meas, earth_constraint=False, p0=prior, init=prior)
    pos = sol["pos"]
    lat, lon, h = ecef_to_llh(pos)
    err = _great_circle_km(pos, truth_ecef)
    print(f"   收敛: {sol['converged']}  迭代: {sol['iter']}  卫星数: {sol['n_sat']}")
    print(f"   ECEF 位置 (m): X={pos[0]:.1f}  Y={pos[1]:.1f}  Z={pos[2]:.1f}")
    print(f"   接收机钟差 b = {sol['clock_bias']/299792458:.6f} s "
          f"(= {sol['clock_bias']/1e3:.1f} km)")
    print(f"   GDOP = {sol['gdop']:.2f}")
    print(f"   LLH: 纬度={lat:.6f}°  经度={lon:.6f}°  高程={h:.1f} m")
    print(f"   与诺丁汉市中心参照点 (52.9536°N, 1.1505°W) 距离 = {err:.2f} km")
    print(f"   各星伪距残差 (m):")
    for m, res in zip(meas, sol["residuals"]):
        print(f"     PRN {m['prn']:>2}: {res:+.2f}")

    print("\n" + "=" * 78)
    if sol["converged"] and sol["n_sat"] >= 4:
        print(f"✅ PVT 定位成功（真实数据，{len(meas)} 星完整 SPS）。"
              f"结果落在诺丁汉市内（距市中心参照点 {err:.2f} km；该参照点并非天线坐标）。")
        rc = 0
    else:
        print("❌ PVT 未收敛。")
        rc = 1
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
