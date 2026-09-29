#!/usr/bin/env python3
"""PVT 结果检查器 —— 复用 07 脚本缓存的伪距观测，快速重解并做详细诊断。

因为 30s 跟踪很慢，本脚本只读取 data/processed/sige_pvt_measurements.pkl，
不重新跟踪，便于在定位异常时快速迭代（换初猜、关/开地球约束、打印卫星几何等）。

运行：
    uv run python scripts/08_pvt_inspect.py
"""
from __future__ import annotations

import pickle
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from gnssrx.pvt import build_measurement, ecef_to_llh, solve, solve_robust  # noqa: E402

CACHE = ROOT / "data" / "processed" / "sige_pvt_measurements.pkl"


def great_circle_km(a: np.ndarray, b: np.ndarray) -> float:
    ra, rb = np.linalg.norm(a), np.linalg.norm(b)
    cosang = float(np.dot(a, b) / (ra * rb))
    cosang = max(-1.0, min(1.0, cosang))
    return 6371.0 * np.arccos(cosang)


def main() -> int:
    if not CACHE.exists():
        print(f"找不到缓存 {CACHE}，请先运行 07_sige_pvt_test.py --regen")
        return 1
    with open(CACHE, "rb") as f:
        meas = pickle.load(f)

    print(f"缓存中卫星数: {len(meas)}")
    print("-" * 78)
    for m in meas:
        sat = np.asarray(m["sat_ecef"])
        print(f"PRN {m['prn']:>2}  ρ={m['pseudorange']:.4e} m  "
              f"sat_clk={m['sat_clock_m']:+.1f} m  t_sat={m['t_sat']:.1f}s")
        print(f"        sat ECEF = [{sat[0]:.0f}, {sat[1]:.0f}, {sat[2]:.0f}]  "
              f"|r|={np.linalg.norm(sat)/1e3:.1f} km")

    # 用 3 星地球约束求解（先消除整数毫秒模糊度）
    earth = (len(meas) == 3)
    print("-" * 78)
    print(f"求解模式: {'3 星 + 地球表面约束 (+整数毫秒模糊度消除)' if earth else '完整 SPS (≥4 星)'}")
    sol = solve_robust(meas, earth_constraint=earth)
    pos = sol["pos"]
    lat, lon, h = ecef_to_llh(pos)
    print(f"收敛={sol['converged']}  迭代={sol['iter']}  n_sat={sol['n_sat']}")
    print(f"ECEF: X={pos[0]:.1f} Y={pos[1]:.1f} Z={pos[2]:.1f}")
    print(f"LLH:  纬度={lat:.6f}° 经度={lon:.6f}° 高程={h:.1f} m")
    print(f"接收机钟差 b = {sol['clock_bias']:.3e} m  →  b/c = {sol['clock_bias']/299792458:.3f} s")
    print(f"GDOP={sol['gdop']:.2f}")
    print("残差 (m):")
    for m, r in zip(meas, sol["residuals"]):
        print(f"   PRN {m['prn']:>2}: {r:+.2f}")

    # 参考点：Nottingham, UK (~52.95°N, 1.15°W) 来自 SiGe 静态采集
    ref = np.array([
        -3.858e6 + 0,  # 占位，下面用经纬度算
    ])
    import math
    lat0, lon0 = math.radians(52.95), math.radians(-1.15)
    a, f = 6378137.0, 1.0/298.257223563
    e2 = f*(2-f)
    n0 = a/np.sqrt(1-e2*np.sin(lat0)**2)
    ref_pt = np.array([
        (n0)*np.cos(lat0)*np.cos(lon0),
        (n0)*np.cos(lat0)*np.sin(lon0),
        (n0*(1-e2))*np.sin(lat0),
    ])
    print("-" * 78)
    dkm = great_circle_km(pos, ref_pt)
    print(f"与参考点(Nottingham ~52.95°N,1.15°W)球面距离 ≈ {dkm:.1f} km")

    # 尝试一个远离初猜，验证解是否稳定
    if len(meas) >= 4:
        for trial in [np.array([0.0, 0.0, 6371e3]),
                     np.array([1.4e7, 0.0, 0.0]),
                     np.array([-1.2e7, 5e6, -8e6])]:
            s2 = solve(meas, init=trial)
            p2 = s2["pos"]
            d = np.linalg.norm(p2 - pos)
            print(f"   初猜 {trial} → 解偏移 {d:.1f} m  (收敛={s2['converged']})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
