#!/usr/bin/env python3
"""3 星几何定界：给定 3 颗卫星的 ECEF 与（消模糊后的）伪距，求 3 个球面的
两个交点，判断 149.7 km 偏差的来源。

三种可能：
  (a) 3 球交点是两个，接收机真实位置是其中之一；我们解的恰是"另一个" → 3 星固有歧义。
  (b) 两个交点都离参考点很远 → 参考点（Nottingham 假设）本身不对，或伪距有偏差。
  (c) 伪距本身量错了（整数毫秒歧义没消干净 / 码相位误差）→ 交点整体上移。

直接复用 07 缓存的伪距观测（已含消模糊）。
"""
from __future__ import annotations

import pickle
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from gnssrx.pvt import resolve_ms_ambiguity, ecef_to_llh, R_EARTH  # noqa: E402
from gnssrx.ephemeris import decode_subframes                          # noqa: E402

CACHE = ROOT / "data" / "processed" / "sige_pvt_measurements.pkl"

# 参考点（假设采集地）：Nottingham
REF_LAT, REF_LON = 52.95, -1.15


def llh_to_ecef(lat_deg, lon_deg, h=0.0):
    a = 6378137.0
    f = 1.0 / 298.257223563
    e2 = f * (2 - f)
    lat = np.radians(lat_deg)
    lon = np.radians(lon_deg)
    N = a / np.sqrt(1 - e2 * np.sin(lat) ** 2)
    x = (N + h) * np.cos(lat) * np.cos(lon)
    y = (N + h) * np.cos(lat) * np.sin(lon)
    z = (N * (1 - e2) + h) * np.sin(lat)
    return np.array([x, y, z])


def intersect_three_spheres(sats, ranges):
    """3 个球心 sats[i]、半径 ranges[i] 的两个交点（解析）。返回 (p1, p2)。"""
    P0, P1, P2 = (np.asarray(s, float) for s in sats)
    r0, r1, r2 = (float(r) for r in ranges)
    # 以 P0 为原点，P1 在 +x，P2 在 xy 平面
    e_x = (P1 - P0) / np.linalg.norm(P1 - P0)
    i = np.dot(e_x, P2 - P0)
    e_y = (P2 - P0 - i * e_x)
    e_y = e_y / np.linalg.norm(e_y)
    e_z = np.cross(e_x, e_y)
    d = np.linalg.norm(P1 - P0)
    j = np.dot(e_y, P2 - P0)
    x = (r0 ** 2 - r1 ** 2 + d ** 2) / (2 * d)
    y = (r0 ** 2 - r2 ** 2 + i ** 2 + j ** 2) / (2 * j) - (i / j) * x
    z2 = r0 ** 2 - x ** 2 - y ** 2
    if z2 < -1e-6:
        return None, None
    z = np.sqrt(max(z2, 0.0))
    center = P0 + x * e_x + y * e_y
    p1 = center + z * e_z
    p2 = center - z * e_z
    return p1, p2


def great_circle_km(a, b):
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    R = 6371.0
    lat1, lon1, _ = ecef_to_llh(a)
    lat2, lon2, _ = ecef_to_llh(b)
    p1 = np.radians([lat1, lon1])
    p2 = np.radians([lat2, lon2])
    dlat = p2[0] - p1[0]
    dlon = p2[1] - p1[1]
    h = np.sin(dlat / 2) ** 2 + np.cos(p1[0]) * np.cos(p2[0]) * np.sin(dlon / 2) ** 2
    return 2 * R * np.arcsin(np.sqrt(h))


def main() -> int:
    with open(CACHE, "rb") as f:
        meas = pickle.load(f)
    print(f"缓存观测: {[(m['prn'], round(m['pseudorange'])) for m in meas]}")

    ref = llh_to_ecef(REF_LAT, REF_LON)

    for label, prior in [("Nottingham 先验", ref),
                         ("质心投影先验", None)]:
        if prior is None:
            cvec = np.mean([np.asarray(m["sat_ecef"]) for m in meas], axis=0)
            prior = cvec / np.linalg.norm(cvec) * R_EARTH
        resolved = resolve_ms_ambiguity(meas, prior)
        sats = [np.asarray(m["sat_ecef"]) for m in resolved]
        ranges = [m["pseudorange"] for m in resolved]
        print(f"\n=== {label} ===")
        for m, rr in zip(resolved, ranges):
            print(f"  PRN {m['prn']:>2}: N={m['_ambig_n']}  ρ={rr/1e3:.1f} km  "
                  f"|sat|={np.linalg.norm(m['sat_ecef'])/1e3:.1f} km")
        p1, p2 = intersect_three_spheres(sats, ranges)
        if p1 is None:
            print("  三球面不相交（伪距内部不一致）")
            continue
        for name, p in [("交点A", p1), ("交点B", p2)]:
            lat, lon, h = ecef_to_llh(p)
            print(f"    {name}: lat={lat:.4f}° lon={lon:.4f}° h={h/1e3:.2f} km  "
                  f"距参考 {great_circle_km(p, ref):.1f} km")

    # 直接用 LS 解（earth constraint）看收敛到哪个交点
    from gnssrx.pvt import solve
    sol = solve(resolved, init=ref, earth_constraint=True)
    print(f"\nLS(earth) 解: lat={ecef_to_llh(sol['pos'])[0]:.4f}° "
          f"lon={ecef_to_llh(sol['pos'])[1]:.4f}° "
          f"h={ecef_to_llh(sol['pos'])[2]/1e3:.2f} km")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
