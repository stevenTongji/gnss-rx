#!/usr/bin/env python3
"""用 IGS 广播星历（RINEX 导航文件）独立校验我们解出的星历。

为什么需要它
------------
我们上一轮把 SiGe 数据集的元数据坐标（德国慕尼黑）当成"精确真值"，但用**解码出的
星历**反算发现：7 颗被跟踪的卫星里有 4 颗在该坐标处位于地平线以下（−6°/−15°/−19°/−27°），
而斜距高达 26,700–28,700 km（地面接收机到 GPS 卫星的物理上限只有 ≈25,776 km）。
这有两种可能：

    A. 我们的星历解码有静默位错（`self_check` 只查 IODE 一致性与 sqrtA/e/i0 的量级，
       查不出 M0/Ω0 这类"位置相关"参数的位错）；
    B. 元数据里的 `<position>` 是占位值（该条目的 `<campaign>` 写的是 "Demo data"），
       数据其实不是在那里采的。

本脚本用**权威第三方源**判决：从 BKG 的 IGS 归档下载当天（DOY 由 GPS 周内秒算出）
的广播星历 RINEX，逐参数比对。若完全一致 → 我们的解码没问题（结论 B）；若有差异
→ 是我们的 bug（结论 A），可据此定位。

用法：
    uv run python scripts/18_igs_ephemeris_check.py            # 校验 SiGe
    uv run python scripts/18_igs_ephemeris_check.py --prn 30 29 31 1 21 5 13 23 16 --week 1609 --tow 466752
"""
from __future__ import annotations

import datetime as dt
import gzip
import math
import pickle
import subprocess
import sys
import urllib.request
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from gnssrx.ephemeris import decode_subframes, read_tow, strip_parity  # noqa: E402
from gnssrx.nav_msg import BITS_PER_SUBFRAME, bit_sync, extract_bits, \
    find_subframes                                                     # noqa: E402
from gnssrx.pvt import C_LIGHT, ecef_to_llh, llh_to_ecef                # noqa: E402

CACHE = ROOT / "data" / "raw" / "igs"
MIRRORS = [
    "https://igs.bkg.bund.de/root_ftp/IGS/BRDC/{y}/{doy3}/brdc{doy3}0.{yy}n.Z",
    "https://igs.bkg.bund.de/root_ftp/IGS/BRDC/{y}/brdc{doy3}0.{yy}n.Z",
]
GPS_EPOCH = dt.datetime(1980, 1, 6, tzinfo=dt.timezone.utc)


# ------------------------------------------------------------------ 下载/解析


def fetch_brdc(week: int, sow: float) -> Path | None:
    """按 GPS 周内秒算出日期，下载当天的 IGS 广播星历（RINEX 2 nav）。"""
    t = GPS_EPOCH + dt.timedelta(weeks=week, seconds=sow)
    doy3, yy = f"{int(t.strftime('%j')):03d}", f"{t:%y}"
    CACHE.mkdir(parents=True, exist_ok=True)
    out = CACHE / f"brdc{doy3}0.{yy}n"
    if out.exists() and out.stat().st_size > 1000:
        return out
    print(f"  下载 IGS 广播星历：{t:%Y-%m-%d} (DOY {doy3}) …")
    for m in MIRRORS:
        url = m.format(y=t.year, doy3=doy3, yy=yy)
        try:
            raw = urllib.request.urlopen(url, timeout=60).read()
        except Exception as exc:                                   # noqa: BLE001
            print(f"    ✗ {url}  ({exc.__class__.__name__})")
            continue
        if len(raw) < 1000:
            continue
        if raw[:2] == b"\x1f\x8b":                # gzip
            text = gzip.decompress(raw)
        elif raw[:2] == b"\x1f\x9d":              # Unix compress (.Z)
            tmp = out.with_suffix(".Z")
            tmp.write_bytes(raw)
            r = subprocess.run(["zcat", str(tmp)], capture_output=True)
            text = r.stdout
        else:
            text = raw
        if len(text) > 1000:
            out.write_bytes(text)
            print(f"    ✓ {url}  ({len(text)/1e3:.0f} kB)")
            return out
    print("    ✗ 所有镜像都失败")
    return None


def _fld(s: str) -> float:
    """RINEX 数值字段：D 指数；空白字段返回 NaN（第 7 行末尾常为空白）。"""
    t = s.replace("D", "E").replace("d", "e").replace("\x00", " ").strip()
    return float(t) if t else float("nan")


def parse_rinex_nav(path: Path) -> list[dict]:
    """解析 RINEX 2 广播星历（每 8 行一条记录），返回参数 dict 列表。"""
    lines = path.read_text(errors="ignore").splitlines()
    try:
        i = next(k for k, ln in enumerate(lines) if "END OF HEADER" in ln) + 1
    except StopIteration:
        return []
    out = []
    while i < len(lines):
        ln = lines[i]
        if len(ln) < 22 or not ln[:2].strip().isdigit():
            i += 1
            continue
        try:
            rec = {
                "prn": int(ln[0:2]),
                "af0": _fld(ln[22:41]), "af1": _fld(ln[41:60]),
                "af2": _fld(ln[60:79]),
            }
            body = []
            for k in range(1, 8):
                bl = lines[i + k]
                body += [_fld(bl[3 + 19 * j: 3 + 19 * (j + 1)]) for j in range(4)]
            # 前 24 个字段是标准参数；第 7 行余下的是发播时刻/拟合区间（末尾可能空白）
            (rec["iode"], rec["crs"], rec["deln"], rec["m0"],
             rec["cuc"], rec["e"], rec["cus"], rec["sqrta"],
             rec["toe"], rec["cic"], rec["omg0"], rec["cis"],
             rec["i0"], rec["crc"], rec["omg"], rec["omgd"],
             rec["idot"], _c2, rec["week"], _l2,
             _acc, _hlt, _tgd, rec["iodc"]) = body[:24]
            out.append(rec)
            i += 8
        except Exception:                                          # noqa: BLE001
            i += 1
    return out


# ------------------------------------------------------------------ 我们的星历


def our_ephemerides(track_pkl: Path) -> dict[int, dict]:
    tr = pickle.load(open(track_pkl, "rb"))
    out = {}
    for prn, t in tr.items():
        ip = t["ip"][500:]
        off, _ = bit_sync(ip)
        bits, _ = extract_bits(ip, off)
        chains = find_subframes(bits)
        if not chains:
            continue
        pairs = [(s, bits[s:s + BITS_PER_SUBFRAME].astype(np.int8))
                 for s in chains[0]["starts"] if s + BITS_PER_SUBFRAME <= len(bits)]
        eph, _ = decode_subframes(prn, [b for _, b in pairs])
        if eph is None or eph.self_check():
            continue
        tow = read_tow(strip_parity(pairs[-1][1])[0])
        out[prn] = {"eph": eph, "tow": tow, "tsat": tow - 6.0}
    return out


# ------------------------------------------------------------------ 比对


# 每个参数换算成"物理量纲"后再比，避免被 RINEX 12 位有效数字的打印精度误导
FIELDS = [("sqrtA", "sqrta", 1.0, " m"), ("e", "e", 1.0, ""),
          ("i0", "i0", 180 / math.pi, " deg"), ("M0", "m0", 180 / math.pi, " deg"),
          ("OMG0", "omg0", 180 / math.pi, " deg"), ("OMGd", "omgd", 180 / math.pi, " deg/s"),
          ("omg", "omg", 180 / math.pi, " deg"),
          ("deln", "deln", 1e9, "e-9 rad/s"), ("idot", "idot", 180 / math.pi, " deg/s"),
          ("cuc", "cuc", 1e6, "e-6 rad"), ("cus", "cus", 1e6, "e-6 rad"),
          ("crc", "crc", 1.0, " m"), ("crs", "crs", 1.0, " m"),
          ("cic", "cic", 1e6, "e-6 rad"), ("cis", "cis", 1e6, "e-6 rad"),
          ("af0", "af0", 1e9, " ns"), ("af1", "af1", 1e12, "e-12 s/s"),
          ("af2", "af2", 1e18, "e-18 1/s2")]


def _sat_pos_from_params(v: dict, tk: float) -> np.ndarray:
    """用一组参数（RINEX 键名）算卫星 ECEF，公式同 RTKLIB eph2pos。"""
    return _eph2pos_rnx(v, tk)


def compare_to_igs(ours: dict, recs: list[dict]) -> None:
    """逐星：找出 IGS 中 toe 相同的记录，比较参数与由此算出的卫星位置。"""
    print(f"\n{'PRN':>4} {'IODE(我们/IGS)':>15} {'Δtoe(s)':>8} "
          f"{'|Δ卫星位置|':>12}  最大参数差（物理量纲）")
    nbad = 0
    for prn in sorted(ours):
        eph, tsat = ours[prn]["eph"], ours[prn]["tsat"]
        cands = [r for r in recs if r["prn"] == prn]
        if not cands:
            print(f"{prn:>4} {'—':>15} {'—':>8} {'—':>12}  （IGS 文件中无此 PRN）")
            continue
        r = min(cands, key=lambda x: abs(x["toe"] - eph.toes))
        dtoe = r["toe"] - eph.toes
        # 用"RINEX 键名"构造我们的参数集，才能喂给同一个位置函数做对比
        mine = {theirs_k: float(getattr(eph, mine_k))
                for mine_k, theirs_k, _, _ in FIELDS}
        mine["toe"] = float(eph.toes)
        pos_o = _sat_pos_from_params(mine, tsat - eph.toes)
        pos_i = _sat_pos_from_params(r, tsat - r["toe"])
        dpos = float(np.linalg.norm(pos_o - pos_i))
        worst, wval = "", 0.0
        for mine_k, theirs_k, scale, unit in FIELDS:
            d = abs(float(getattr(eph, mine_k)) - float(r[theirs_k])) * scale
            if d > wval:
                wval, worst = d, f"{mine_k}: {d:.3g}{unit}"
        flag = "" if dpos < 1.0 else "  ← 不一致"
        if dpos >= 1.0:
            nbad += 1
        print(f"{prn:>4} {str(eph.iode2)+'/'+str(int(r['iode'])):>15} "
              f"{dtoe:>8.0f} {dpos*1e3:>9.1f} mm  {worst}{flag}")
    print(f"\n位置不一致（>1 mm 视为异常）的卫星数：{nbad}/{len(ours)}")


def main() -> int:
    args = sys.argv[1:]
    if "--track" in args:
        track = Path(args[args.index("--track") + 1])
    else:
        track = ROOT / "data" / "processed" / "sige_tracking.pkl"
    ours = our_ephemerides(track)
    week = ours[next(iter(ours))]["eph"].week
    tow = ours[next(iter(ours))]["tow"]
    print(f"我们解出 {len(ours)} 颗星：PRN {sorted(ours)}；GPS 周 {week}，TOW {tow:.0f} s")

    nav = fetch_brdc(week, tow)
    if nav is None:
        return 1
    recs = parse_rinex_nav(nav)
    print(f"  IGS 广播星历：{len(recs)} 条记录，覆盖 PRN "
          f"{sorted({r['prn'] for r in recs})}")

    compare_to_igs(ours, recs)

    # 顺带：用 IGS 星历独立算一次"从某坐标能看到哪些星"
    print("\n附：用 IGS 星历从各候选坐标检查可见性（仰角 > 5° 计为可见）")
    sites = {"慕尼黑（元数据）": (48.17154012, 11.80868949, 576.86),
             "德州/路易斯安那（量测指向）": (30.0, -94.1, 60.0)}
    for name, (la, lo, h) in sites.items():
        P = llh_to_ecef(la, lo, h)
        els = []
        for prn in sorted(ours):
            cands = [r for r in recs if r["prn"] == prn]
            if not cands:
                continue
            r = min(cands, key=lambda x: abs(x["toe"] - ours[prn]["eph"].toes))
            tk = ours[prn]["tsat"] - r["toe"]
            s = _eph2pos_rnx(r, tk)
            els.append((prn, _elev(P, la, lo, s)))
        vis = [p for p, e in els if e > 5.0]
        print(f"  {name:<26} 可见 {len(vis)}/{len(els)}  "
              f"({', '.join(f'{p}:{e:.0f}°' for p, e in els)})")
    return 0


def _eph2pos_rnx(r: dict, tk: float) -> np.ndarray:
    """用 RINEX 参数算 ECEF（公式同 RTKLIB eph2pos），独立于我们自己的星历类。"""
    mu, omge, A = 3.9860050e14, 7.2921151467e-5, r["sqrta"] ** 2
    M = r["m0"] + (math.sqrt(mu / A ** 3) + r["deln"]) * tk
    E = M
    for _ in range(30):
        dE = (E - r["e"] * math.sin(E) - M) / (1 - r["e"] * math.cos(E))
        E -= dE
        if abs(dE) < 1e-14:
            break
    sE, cE = math.sin(E), math.cos(E)
    u = math.atan2(math.sqrt(1 - r["e"] ** 2) * sE, cE - r["e"]) + r["omg"]
    rad = A * (1 - r["e"] * cE)
    inc = r["i0"] + r["idot"] * tk
    s2, c2 = math.sin(2 * u), math.cos(2 * u)
    u += r["cus"] * s2 + r["cuc"] * c2
    rad += r["crs"] * s2 + r["crc"] * c2
    inc += r["cis"] * s2 + r["cic"] * c2
    x, y = rad * math.cos(u), rad * math.sin(u)
    ci = math.cos(inc)
    O = r["omg0"] + (r["omgd"] - omge) * tk - omge * r["toe"]
    sO, cO = math.sin(O), math.cos(O)
    return np.array([x * cO - y * ci * sO, x * sO + y * ci * cO, y * math.sin(inc)])


def _elev(P: np.ndarray, lat: float, lon: float, sat: np.ndarray) -> float:
    la, lo = math.radians(lat), math.radians(lon)
    up = np.array([math.cos(la) * math.cos(lo), math.cos(la) * math.sin(lo),
                   math.sin(la)])
    d = np.asarray(sat, float) - np.asarray(P, float)
    return math.degrees(math.asin(max(-1.0, min(1.0, float(d @ up)
                                                 / float(np.linalg.norm(d))))))


if __name__ == "__main__":
    raise SystemExit(main())
