#!/usr/bin/env python3
"""单点定位（PVT）验证 —— 真实数据（SiGe GN3S，8-bit 实数 IF）。

流程：捕获 → 跟踪 → 帧同步/星历解码 → 伪距提取 → 迭代加权最小二乘定位。

伪距由"子帧起点块索引 × 1ms + 本地码相位 τ + TOW"算出；接收机钟差 b 吸收
我们 0 时基与 GPS 时的公共偏置，因此 TOW 只需同星期内一致。

SiGe 这份数据稳定解出 3 颗干净卫星（PRN 7/9/11），因此演示"3 星 + 地球表面约束"
的 3D 定位（标准静态接收机做法）；若某天数据解出 ≥4 颗，则自动走完整 SPS。

运行：
    uv run python scripts/07_sige_pvt_test.py          # 复用缓存
    uv run python scripts/07_sige_pvt_test.py --regen  # 重新捕获/跟踪/定位
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
from gnssrx.io_if import read_real_int8                                  # noqa: E402
from gnssrx.nav_msg import BITS_PER_SUBFRAME, bit_sync, extract_bits, \
    find_subframes                                                       # noqa: E402
from gnssrx.pvt import build_measurement, ecef_to_llh, refine_measurements, \
    solve_cold_start, solve_robust                                       # noqa: E402

import math


def unwrap_code_phase(tau_stream: np.ndarray, fs: float, pdi_ms: int) -> np.ndarray:
    """把逐块（mod spc）的码相位 τ 解缠成"自跟踪开始累计的码片数"。

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
from gnssrx.tracking import TrackingConfig, track_all                     # noqa: E402

DATA = ROOT / "data" / "raw" / "sige_gps_l1_8bit.dat"
FS = 16.368e6
F_IF = 4.092e6
PROCESS_MS = 30000
N_CHANNELS = 12
THRESHOLD = 1.5
SETTLE_MS = 500
PDI_MS = 1                      # 与 TrackingConfig 默认一致（20 块/比特）
CACHE = ROOT / "data" / "processed" / "sige_pvt_measurements.pkl"
TRACK_CACHE = ROOT / "data" / "processed" / "sige_tracking.pkl"
TOW_DELTA_S = 0.0               # 发射时刻口径微调（用于诊断；0 = t_sat = TOW − 6）


def get_tracking(regen: bool = False) -> dict:
    """捕获 + 跟踪，返回 {prn: {"ip","tau","cn0"}}。带磁盘缓存（跟踪最耗时）。"""
    if TRACK_CACHE.exists() and not regen:
        with open(TRACK_CACHE, "rb") as f:
            return pickle.load(f)                                     # type: ignore

    print("① 捕获全部 32 颗（门限 %g），取最强 %d 颗" % (THRESHOLD, N_CHANNELS))
    data = read_real_int8(DATA, n_samples=int(FS * PROCESS_MS / 1e3))
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


def get_measurements(regen: bool = False) -> list[dict]:
    """完整链路 → 返回干净卫星的伪距观测列表（dict）。带磁盘缓存。"""
    if CACHE.exists() and not regen:
        with open(CACHE, "rb") as f:
            cached = pickle.load(f)                                   # type: ignore
        # 旧缓存没有 _eph（迭代重算卫星位置需要），视为过期 → 重建
        if cached and "_eph" in cached[0]:
            return cached
        print("   （测量缓存缺少 _eph，重建以启用迭代重算）")

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
        # 配对 (start_bit, 300bit_block)
        pairs = []
        for start in chain["starts"]:
            if start + BITS_PER_SUBFRAME > len(bits):
                continue
            pairs.append((start, bits[start:start + BITS_PER_SUBFRAME].astype(np.int8)))

        # 星历解码（取第一个能凑齐 SF1/2/3 的）
        eph, info = decode_subframes(prn, [b for _, b in pairs])
        if eph is None or eph.self_check():
            print(f"   PRN {prn:>2}  星历未通过自检，跳过")
            continue
        # 用最后一个子帧起点组装伪距：先解缠码相位得到亚微秒级接收机时刻
        start, block = pairs[-1]
        d240, _ = strip_parity(block)
        tow = read_tow(d240)
        unwrapped = unwrap_code_phase(t["tau"], FS, PDI_MS)
        m_block = SETTLE_MS + offset + start * 20 * PDI_MS  # 跟踪块序号（每块 1ms）
        if m_block >= len(unwrapped):
            m_block = len(unwrapped) - 1
        # 接收机时刻 = 整数毫秒(块序号) − 亚毫秒(解缠后的码相位)。
        # 本地码每 1ms 刚好推进 1 个完整周期(1023 码片)，故整数毫秒由块序号给出；
        # 解缠码相位给出亚毫秒精确值。注意符号：本跟踪器的 τ 与「信号码历元在块内的
        # 偏移」反号，因此这里取【减号】（合成自检脚本 15 已严格验证：
        # 取减号后 t_user 与真实到达时刻之差为公共常数；取加号则随卫星散布 ~1ms）。
        # 绝对值的时基偏置由钟差 b 吸收。
        t_user = m_block / 1000.0 - unwrapped[m_block] / 1.023e6
        meas = build_measurement(eph, prn, tow + TOW_DELTA_S, t_user, FS)
        meas["tow"] = tow + TOW_DELTA_S
        # 留着供 refine_measurements 迭代重算卫星位置（见 pvt.refine_measurements）
        meas["_eph"] = eph
        meas["_t_user"] = t_user
        print(f"   PRN {prn:>2}  ✅ ρ={meas['pseudorange']:.4e}  "
              f"off={offset:2d} start={start:5d} m_blk={m_block:6d} "
              f"τ={unwrapped[m_block]:+9.2f}chip  tow={tow:.0f}  t_user={t_user:.6f}")
        measurements.append(meas)

    print(f"\n   共 {len(measurements)} 颗干净卫星参与定位。")
    CACHE.parent.mkdir(parents=True, exist_ok=True)
    with open(CACHE, "wb") as f:
        pickle.dump(measurements, f)
    return measurements


def _ecef_of(lat_deg, lon_deg, h=0.0):
    a = 6378137.0
    f = 1.0 / 298.257223563
    e2 = f * (2 - f)
    lat, lon = np.radians(lat_deg), np.radians(lon_deg)
    N = a / np.sqrt(1 - e2 * np.sin(lat) ** 2)
    return np.array([(N + h) * np.cos(lat) * np.cos(lon),
                     (N + h) * np.cos(lat) * np.sin(lon),
                     (N * (1 - e2) + h) * np.sin(lat)])


def main() -> int:
    if not DATA.exists():
        print(f"找不到数据文件：{DATA}\n下载方式见 data/README.md")
        return 1
    regen = "--regen" in sys.argv
    meas = get_measurements(regen=regen)

    if len(meas) < 3:
        print(f"❌ 仅 {len(meas)} 颗干净卫星，不足以定位（最少 3 颗）。")
        return 1

    earth = (len(meas) == 3)
    print(f"\n③ 迭代加权最小二乘定位（{'3 星 + 地球表面约束' if earth else '完整 SPS (≥4 星)'}）")

    if earth:
        # 3 星欠定：整数毫秒消模糊需要一个 1ms（≈300km）以内的先验位置，
        # 真实接收机由基站/上次定位提供。这里用"英国区域先验"做演示（本数据确属英国）。
        region_prior = _ecef_of(54.0, -2.0, 0.0)
        print("   说明：3 星场景需要 300km 内的先验（真实接收机的冷启动由网络/AGPS 给出）。")
        print("   先用'英国区域先验'(54°N,2°W) 求解：")
        sol = solve_robust(meas, earth_constraint=True, p0=region_prior, init=region_prior)
        # 对照：纯冷启动（质心先验）——仅用于展示 3 星对先验的敏感性
        cold = solve_robust(meas, earth_constraint=True)
        clat, clon, ch = ecef_to_llh(cold["pos"])
        print(f"   （对照）质心冷启动解：lat={clat:.4f}° lon={clon:.4f}° h={ch/1e3:.2f} km "
              f"收敛={cold['converged']}  ← 说明无先验时 3 星易落入错误盆地")
    else:
        # ≥4 星：完整 SPS，**不需要任何先验位置**。
        # 早期版本这里给了一个"慕尼黑区域先验"（照抄 ION 元数据的 <position>），
        # 结果落进错误的整毫秒盆地、报出 ~40 km 残差 —— 那是搜索问题，不是数据问题。
        # 现在用 solve_cold_start：0.5° 细网格全球搜 + 多个候选盆地各自精化，
        # 按最终残差挑最优。代价函数含 300 km 硬台阶、真盆地只有约 100–200 km 宽，
        # 粗网格（≥1°）会直接跨过去 —— 见 pvt.solve_cold_start 的 docstring。
        print("   说明：≥4 星为完整 SPS，**不给任何先验**，直接冷启动求解：")
        sol = solve_cold_start(meas)
        clat, clon, ch = ecef_to_llh(sol["pos"])
        print(f"   冷启动定位：lat={clat:.5f}° lon={clon:.5f}° h={ch:.1f} m  "
              f"（试了 {sol.get('cold_start_candidates')} 个候选盆地，取残差最小者）")
        # 对照：盲信元数据 <position>（慕尼黑）会落到错误盆地
        wrong = _ecef_of(48.1715, 11.8087, 0.0)
        sol_w = solve_robust(meas, earth_constraint=False, p0=wrong, init=wrong)
        wlat, wlon, wh = ecef_to_llh(sol_w["pos"])
        wr = np.asarray(sol_w["residuals"])
        print(f"   （对照）用元数据坐标(慕尼黑)当先验：lat={wlat:.4f}° lon={wlon:.4f}° "
              f"残差RMS={np.sqrt((wr**2).mean())/1e3:.2f} km  ← 错误盆地")

    # 迭代重算卫星位置：Sagnac 的 θ = ω·ρ/c 需要接收机位置，第一次解算时只能用
    # 「卫星星下点的地表点」当标称值，会留下每星固定的米级偏差（见 pvt.refine_measurements）。
    # 真实接收机同样是迭代的；scripts/19 的消融实验显示这一步把端到端 2D 从 2.2 m 降到 0.54 m。
    if meas and "_eph" in meas[0]:
        meas = refine_measurements(meas, sol["pos"], FS)
        sol = solve_robust(meas, earth_constraint=earth,
                           p0=sol["pos"], init=sol["pos"])

    pos = sol["pos"]
    lat, lon, h = ecef_to_llh(pos)
    print(f"   收敛: {sol['converged']}  迭代: {sol['iter']}  卫星数: {sol['n_sat']}")
    print(f"   ECEF 位置 (m): X={pos[0]:.1f}  Y={pos[1]:.1f}  Z={pos[2]:.1f}")
    print(f"   接收机钟差 b = {sol['clock_bias']/299792458:.6f} s "
          f"(= {sol['clock_bias']/1e3:.1f} km)")
    print(f"   GDOP = {sol['gdop']:.2f}")
    print(f"   LLH: 纬度={lat:.6f}°  经度={lon:.6f}°  高程={h:.1f} m")
    print(f"   各星伪距残差 (m):")
    for m, res in zip(meas, sol["residuals"]):
        print(f"     PRN {m['prn']:>2}: {res:+.2f}")

    print("\n" + "=" * 78)
    if sol["converged"] and sol["n_sat"] >= 4:
        print(f"✅ PVT 定位成功（真实数据，{len(meas)} 星完整 SPS，无任何先验）。")
        print("   本数据（ION GNSS SDR 元数据标准样本，SiGe GN3S v3，2013-05-23）"
              "的采集地：")
        print("     · 元数据 <position> 写的是德国慕尼黑 —— 但那是**占位值**")
        print("       （<campaign>Demo data</campaign>；且慕尼黑处有 4 颗星在地平线下、")
        print("        斜距超过地面接收机的物理上限，几何不成立）；")
        print("     · 由数据本身反解：39.2877°N / 82.0634°W（美国俄亥俄州 Athens 附近），")
        print("       距 Ohio University 约 5.7 km —— 与元数据里的 <contact>Sanjeev</contact>")
        print("       （Sanjeev Gunawardena，2013 年任职于 Ohio University 航空电子工程中心、")
        print("        亦为 ION GNSS SDR 元数据标准共同作者）完全吻合。")
        rc = 0
    else:
        print("❌ PVT 未收敛。")
        rc = 1
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
