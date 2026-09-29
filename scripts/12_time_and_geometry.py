#!/usr/bin/env python3
"""核验时间基与卫星几何：解码星历里的 WN/TOW，确认卫星位置算在正确的 GPS 时刻；
并计算三颗星的星下点，判断接收机究竟被几何"逼"到了英国还是美国。"""
from __future__ import annotations

import pickle
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from gnssrx.ephemeris import decode_subframes, read_tow, strip_parity  # noqa: E402
from gnssrx.pvt import ecef_to_llh, R_EARTH, solve, resolve_ms_ambiguity  # noqa: E402

EPH_CACHE = ROOT / "data" / "processed" / "sige_ephemeris_blocks.pkl"
MEAS_CACHE = ROOT / "data" / "processed" / "sige_pvt_measurements.pkl"

GPS_EPOCH = datetime(1980, 1, 6, 0, 0, 0, tzinfo=timezone.utc)


def gps_week_seconds_to_utc(wn: int, tow: int) -> datetime:
    return GPS_EPOCH + __import__("datetime").timedelta(weeks=wn, seconds=tow)


def subpoint(eph, t_sat):
    """卫星星下点（纬度、经度，单位度）。"""
    pos = np.asarray(eph.position(t_sat - eph.toes), float)
    # ECEF → 地理（简化，忽略地球扁率对星下点影响极小）
    x, y, z = pos
    lon = np.degrees(np.arctan2(y, x))
    p = np.hypot(x, y)
    lat = np.degrees(np.arctan2(z, p))
    return lat, lon


def main() -> int:
    with open(EPH_CACHE, "rb") as f:
        blocks_by_prn = pickle.load(f)
    # 取一颗干净星解码 WN / TOW / 星历时刻
    for prn in (7, 9, 11):
        if prn not in blocks_by_prn:
            continue
        eph, info = decode_subframes(prn, blocks_by_prn[prn])
        if eph is None or eph.self_check():
            continue
        # 找一块有 TOW 的
        blk = blocks_by_prn[prn][0]
        d240, _ = strip_parity(blk)
        tow = read_tow(d240)
        print(f"PRN {prn}: WN={eph.week}  TOW={tow}s  toc={eph.toc}s  toes={eph.toes}s")
        dt = gps_week_seconds_to_utc(eph.week, tow)
        print(f"   对应 UTC 时刻: {dt.isoformat()}  (采集自 2013-05-23 的 ION SiGe 样本)")
        # 发射时刻 t_sat = TOW-6
        t_sat = tow - 6.0
        lat, lon = subpoint(eph, t_sat)
        print(f"   卫星星下点 @t_sat: lat={lat:.2f}° lon={lon:.2f}°  高度={np.linalg.norm(eph.position(t_sat-eph.toes))/1e3:.1f} km")
        break

    # 三颗星星下点 + 接收机几何
    print("\n=== 三颗卫星星下点（判断几何所在半球）===")
    sats, ranges = [], []
    with open(MEAS_CACHE, "rb") as f:
        meas = pickle.load(f)
    ref_uk = np.array([ecef_of(52.95, -1.15)])
    # 先把伪距消模糊（用英国先验）
    resolved = resolve_ms_ambiguity(meas, ref_uk[0])
    for m in resolved:
        print(f"  PRN {m['prn']:>2}: 星下点 lat={subpoint(eph_prn(m['prn'], blocks_by_prn), m['t_sat'])[0]:.2f}° "
              f"lon={subpoint(eph_prn(m['prn'], blocks_by_prn), m['t_sat'])[1]:.2f}°  "
              f"ρ={m['pseudorange']/1e3:.1f} km")
    # 两个不同初值的 3 星+地球解
    print("\n=== 两个不同初值的 3 星+地球解（看是否出现两个交点）===")
    for name, init in [("英国(Nottingham)", ref_uk[0]),
                       ("美国(Boulder,CO)", ecef_of(40.0, -105.27)),
                       ("澳洲(Perth)", ecef_of(-31.95, 115.86))]:
        sol = solve(resolved, init=init, earth_constraint=True)
        lat, lon, h = ecef_to_llh(sol["pos"])
        print(f"  初值 {name:>18s}: lat={lat:.4f}° lon={lon:.4f}° h={h/1e3:.2f} km "
              f"收敛={sol['converged']} 钟差={sol['clock_bias']/299792458*1e3:.1f} km")
    return 0


def ecef_of(lat_deg, lon_deg, h=0.0):
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


_eph_cache = {}


def eph_prn(prn, blocks_by_prn):
    if prn not in _eph_cache:
        e, _ = decode_subframes(prn, blocks_by_prn[prn])
        _eph_cache[prn] = e
    return _eph_cache[prn]


if __name__ == "__main__":
    raise SystemExit(main())
