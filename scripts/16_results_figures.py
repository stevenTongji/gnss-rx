#!/usr/bin/env python3
"""结果配图 —— 把 PVT 的关键结论画成图，供 README 引用。

出三张图（全部由已有缓存复现，不做重新捕获/跟踪）：

    docs/figures/pvt_result.png     定位结果：局部地图 + 各星伪距残差 + 卫星天空图
    docs/figures/pvt_coldstart.png  冷启动：全球代价地形（1° 网格）+ 采集点附近放大
    docs/figures/pvt_accuracy.png   精度阶梯：合成真值（厘米级）→ 真实数据（米/公里级）

风格与 00–03 的图一致：白底、中文字体（见 `gnssrx.plotting`）。

运行：
    uv run python scripts/16_results_figures.py
"""
from __future__ import annotations

import pickle
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from gnssrx import plotting                                              # noqa: E402,F401
from gnssrx.pvt import C_LIGHT, _grid_ecef, ecef_to_llh, llh_to_ecef, \
    solve_cold_start, solve_robust                                       # noqa: E402

PROC = ROOT / "data" / "processed"
FIGDIR = ROOT / "docs" / "figures"

# 公开采集点（仅作对照，不作为定位输入）
TRUTH = {
    "nottingham": ("Nottingham", 52.95, -1.15, 0.0),
    "sige": ("慕尼黑", 48.17154012, 11.80868949, 576.86),
}
# SiGe 的整数毫秒消模糊需要一个 ~150 km 内的先验（该数据集量测一致性有限，
# 冷启动会落到错误盆地）；真实接收机由上次定位/AGPS 提供该先验。
SIGE_PRIOR = (48.17154012, 11.80868949)

# 合成数据实测结果（取自 `scripts/13_synthetic_4sat_pvt.py` 的标准输出，可复现）
SYNTHETIC = [
    ("合成：6 星，无噪声", 0.010, "#2e9e5b"),
    ("合成：4 星（完整 SPS）", 0.010, "#2e9e5b"),
    ("合成：码相位噪声 σ=2 m", 6.312, "#c9a227"),
    ("合成：码相位噪声 σ=5 m", 15.767, "#d97b28"),
    ("合成：码相位噪声 σ=10 m", 31.525, "#d9534f"),
]


# ---------------------------------------------------------------- 小工具

def great_circle_km(p: np.ndarray, q: np.ndarray) -> float:
    la1, lo1, _ = ecef_to_llh(np.asarray(p, float))
    la2, lo2, _ = ecef_to_llh(np.asarray(q, float))
    la1, lo1, la2, lo2 = map(np.radians, (la1, lo1, la2, lo2))
    c = (np.sin(la1) * np.sin(la2)
         + np.cos(la1) * np.cos(la2) * np.cos(lo1 - lo2))
    return 6371.0 * np.arccos(min(1.0, max(-1.0, float(c))))


def enu_basis(lat_deg: float, lon_deg: float):
    la, lo = np.radians(lat_deg), np.radians(lon_deg)
    e = np.array([-np.sin(lo), np.cos(lo), 0.0])
    n = np.array([-np.sin(la) * np.cos(lo), -np.sin(la) * np.sin(lo), np.cos(la)])
    u = np.array([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)])
    return e, n, u


def enu_of(recv: np.ndarray, lat0: float, lon0: float, target: np.ndarray):
    """target 相对 (lat0, lon0) 的东/北/天分量（米）。"""
    e, n, u = enu_basis(lat0, lon0)
    d = np.asarray(target, float) - np.asarray(recv, float)
    return float(d @ e), float(d @ n), float(d @ u)


def load(name: str):
    with open(PROC / name, "rb") as f:
        return pickle.load(f)


def load_cn0() -> dict[int, float]:
    """从跟踪缓存取各星 C/N₀（dB-Hz）。"""
    out: dict[int, float] = {}
    path = PROC / "nottingham_tracking.pkl"
    if not path.exists():
        return out
    for prn, t in load("nottingham_tracking.pkl").items():
        v = np.asarray(t["cn0"], dtype=float).ravel()
        out[int(prn)] = float(v[-1])
    return out


# ---------------------------------------------------------------- 图 1：定位结果

def fig_result(nott, cn0) -> None:
    sol = solve_cold_start(nott)
    lat, lon, h = ecef_to_llh(sol["pos"])
    name, tlat, tlon, th = TRUTH["nottingham"]
    truth = llh_to_ecef(tlat, tlon, th)
    err_km = great_circle_km(sol["pos"], truth)
    res = np.array(sol["residuals"], float)
    prns = [m["prn"] for m in nott]

    fig = plt.figure(figsize=(15.5, 4.9), dpi=130)

    # ---- (a) 局部地图：东-北（km），真值在原点
    ax = fig.add_subplot(1, 3, 1)
    for r, alpha in ((0.5, 0.55), (1.0, 0.45), (2.0, 0.32)):
        th_ = np.linspace(0, 2 * np.pi, 200)
        ax.plot(r * np.cos(th_), r * np.sin(th_), color="#8a94a6",
                lw=0.7, alpha=alpha)
        ax.text(r * 0.707, r * 0.707 + 0.04, f"{r:g} km", fontsize=7.6,
                color="#8a94a6", ha="left", va="bottom", alpha=0.95)
    ax.plot(0, 0, marker="*", ms=20, color="#1f6feb",
            label=f"公开采集点 ({tlat:.2f}°, {abs(tlon):.2f}°W)")
    se, sn, _ = enu_of(sol["pos"], lat, lon, truth)       # 真值相对解算点
    ax.plot(-se / 1e3, -sn / 1e3, marker="o", ms=8, color="#d9534f",
            label="本接收机解算点（冷启动）")
    ax.annotate(f"{err_km:.2f} km", xy=(-se / 2e3, -sn / 2e3),
                xytext=(-se / 2e3 - 1.15, -sn / 2e3 + 0.62),
                fontsize=10, color="#d9534f",
                arrowprops=dict(arrowstyle="-", color="#d9534f", lw=0.8))
    ax.set_xlim(-2.3, 2.3); ax.set_ylim(-2.3, 2.3)
    ax.set_aspect("equal")
    ax.grid(alpha=0.22)
    ax.set_xlabel("东向偏移 (km)"); ax.set_ylabel("北向偏移 (km)")
    ax.set_title(f"(a) 定位结果（{name}，9 星完整 SPS）\n"
                 f"{lat:.5f}°N, {abs(lon):.5f}°W, h = {h:.1f} m", fontsize=10.5)
    ax.legend(fontsize=8.2, loc="upper left", framealpha=0.9)

    # ---- (b) 各星伪距残差
    ax = fig.add_subplot(1, 3, 2)
    colors = ["#2e9e5b" if abs(v) < 3 else "#d9534f" for v in res]
    ax.bar([str(p) for p in prns], res, color=colors, width=0.62)
    ax.axhline(0, color="#444", lw=0.8)
    ax.axhspan(-2, 2, color="#2e9e5b", alpha=0.10)
    for p, v in zip(prns, res):
        ax.text(str(p), v + (0.18 if v >= 0 else -0.42), f"{v:+.1f}",
                ha="center", va="bottom" if v >= 0 else "top",
                fontsize=8.2, color="#1f2937")
    ax.set_ylim(-2.9, 2.9)
    ax.set_xlabel("PRN"); ax.set_ylabel("伪距残差 (m)")
    ax.set_title(f"(b) 各星伪距残差，RMS = "
                 f"{np.sqrt((res ** 2).mean()):.2f} m\n"
                 f"阴影为 ±2 m 带", fontsize=10.5)
    ax.grid(alpha=0.22, axis="y")

    # ---- (c) 卫星天空图（方位-俯仰，颜色 = C/N₀）
    ax = fig.add_subplot(1, 3, 3, projection="polar")
    azs, els, cs = [], [], []
    for m in nott:
        sat = np.asarray(m["sat_ecef"], float)
        e, n, u = enu_of(sol["pos"], lat, lon, sat)
        rng = float(np.linalg.norm(sat - np.asarray(sol["pos"], float)))
        az = np.degrees(np.arctan2(e, n)) % 360.0
        el = np.degrees(np.arcsin(max(-1.0, min(1.0, u / rng))))
        azs.append(az); els.append(el)
        cs.append(cn0.get(m["prn"], np.nan))
    azs = np.array(azs); els = np.array(els); cs = np.array(cs, float)
    sc = ax.scatter(np.radians(azs), 90.0 - els, c=cs, s=95, cmap="viridis",
                    edgecolors="#333", linewidths=0.6, zorder=3)
    for m, a, el in zip(nott, azs, els):
        ax.text(np.radians(a), 90.0 - el - 4.6, f"{m['prn']}",
                ha="center", va="center", fontsize=7.6, color="#222")
    ax.set_theta_zero_location("N"); ax.set_theta_direction(-1)
    ax.set_rlim(0, 90)
    ax.set_yticks([0, 30, 60], labels=["90°", "60°", "30°"])
    ax.set_rlabel_position(112)
    ax.set_title("(c) 卫星天空图（方位 / 俯仰角）\n"
                 "数字为 PRN，颜色为 C/N₀", fontsize=10.5, pad=14)
    cb = fig.colorbar(sc, ax=ax, pad=0.10, fraction=0.046)
    cb.set_label("C/N₀ (dB-Hz)", fontsize=8.5)
    cb.ax.tick_params(labelsize=8)

    fig.suptitle("GPS L1 单点定位（真实数据 · 冷启动无先验 · Nottingham 9 星）",
                 fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    out = FIGDIR / "pvt_result.png"
    fig.savefig(out); plt.close(fig)
    print(f"  已保存 {out.name}  解算={lat:.5f},{lon:.5f} h={h:.1f}m  "
          f"误差={err_km:.2f} km  残差RMS={np.sqrt((res**2).mean()):.2f} m")


# ---------------------------------------------------------------- 图 2：冷启动

def fig_coldstart(nott) -> None:
    """代价函数结构（实测）：
        距采集点 d 处  cost ≈ 0.55·d（d ≲ 80 km，锥形），
        d ≳ 150 km 后因逐星取整错一个 1 ms（≈300 km）而饱和到 ~80 km 平台。
    因此极小值不"尖"但**唯一**：全球 1° 网格中仅 1 个格点 cost < 10 km。
    """
    sats = np.array([np.asarray(m["sat_ecef"], float) for m in nott])
    raw = np.array([m["pseudorange"] for m in nott], float)
    cms = C_LIGHT * 1e-3

    def cost_km(lats, lons):
        g = _grid_ecef(np.asarray(lats, float), np.asarray(lons, float))
        d = np.linalg.norm(g[:, :, None, :] - sats[None, None, :, :], axis=-1)
        n = np.round((d - raw) / cms)
        return (raw[None, None, :] + n * cms - d).std(axis=-1) / 1e3

    sol = solve_cold_start(nott)
    lat, lon, _ = ecef_to_llh(sol["pos"])
    _, tlat, tlon, _ = TRUTH["nottingham"]

    fig, axes = plt.subplots(1, 2, figsize=(14.8, 5.4), dpi=130)

    # ---- (a) 全球 1° 粗搜地形（唯一的低代价格点即采集点所在处）
    lats = np.arange(-60.0, 71.0, 1.0)
    lons = np.arange(-180.0, 180.0, 1.0)
    z = cost_km(lats, lons)
    jl, jo = np.unravel_index(int(np.argmin(z)), z.shape)
    ax = axes[0]
    im = ax.pcolormesh(lons, lats, np.maximum(z, 3.0), cmap="viridis_r",
                       norm=LogNorm(vmin=3.0, vmax=130.0), shading="auto")
    ax.plot(lons[jo], lats[jl], marker="s", ms=11, mfc="none", mec="#ff4d4d",
            mew=2.0, label=f"全局最小格点 {z[jl, jo]:.1f} km")
    ax.plot(tlon, tlat, marker="*", ms=17, color="#ffd166", mec="#222", mew=0.5,
            label=f"公开采集点 ({tlat:.2f}°, {abs(tlon):.2f}°W)")
    ax.annotate(f"中位代价 {np.median(z):.0f} km", xy=(-150, -40),
                fontsize=9, color="#f0f0f0")
    ax.set_xlim(-180, 180); ax.set_ylim(-60, 70)
    ax.set_xlabel("经度 (°)"); ax.set_ylabel("纬度 (°)")
    ax.set_title("(a) 全球粗搜代价地形（1° 网格）\n"
                 "130×360 个格点中仅 1 个 < 10 km", fontsize=10.5)
    ax.legend(fontsize=8.4, loc="lower left", framealpha=0.92)
    cb = fig.colorbar(im, ax=ax, pad=0.015, fraction=0.045)
    cb.set_label("代价 (km，对数)", fontsize=8.5)
    cb.ax.tick_params(labelsize=8)

    # ---- (b) 采集点附近放大（0.02° ≈ 2 km 网格）：锥形极小值 + 饱和平台
    lats2 = np.arange(50.4, 55.6, 0.02)
    lons2 = np.arange(-5.2, 3.0, 0.02)
    z2 = cost_km(lats2, lons2)
    ax = axes[1]
    im = ax.pcolormesh(lons2, lats2, np.maximum(z2, 0.01), cmap="viridis_r",
                       norm=LogNorm(vmin=0.01, vmax=100.0), shading="auto")
    cs = ax.contour(lons2, lats2, z2, levels=[1.0, 10.0, 30.0, 60.0],
                    colors="#4da3ff", linewidths=0.85)
    ax.clabel(cs, fmt=lambda v: f"{v:g} km", fontsize=7.4, inline=True)
    ax.plot(lon, lat, marker="o", ms=10, mfc="none", mec="#ff4d4d", mew=2.2,
            label="解算点（网格极小值）")
    ax.plot(tlon, tlat, marker="*", ms=18, color="#ffd166", mec="#222", mew=0.5,
            label="公开采集点")
    ax.set_xlabel("经度 (°)"); ax.set_ylabel("纬度 (°)")
    ax.set_title("(b) 采集点附近放大（0.02° ≈ 2 km 网格）\n"
                 f"解算点与采集点相距 "
                 f"{great_circle_km(sol['pos'], llh_to_ecef(tlat, tlon, 0.0)):.2f} km",
                 fontsize=10.5)
    ax.legend(fontsize=8.4, loc="upper left", framealpha=0.92)
    cb = fig.colorbar(im, ax=ax, pad=0.015, fraction=0.045)
    cb.set_label("代价 (km，对数)", fontsize=8.5)
    cb.ax.tick_params(labelsize=8)

    fig.suptitle("冷启动消模糊：不使用任何先验位置 —— 全球唯一的低代价点即采集点所在处",
                 fontsize=12.5)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    out = FIGDIR / "pvt_coldstart.png"
    fig.savefig(out); plt.close(fig)
    print(f"  已保存 {out.name}  全球最小格点={lats[jl]:.0f}°,{lons[jo]:.0f}° "
          f"({z[jl, jo]:.2f} km)")


# ---------------------------------------------------------------- 图 3：精度阶梯

def fig_accuracy(nott, sige) -> None:
    # 真实数据：Nottingham 冷启动（无先验）；SiGe 用 ~150 km 内先验
    sol_n = solve_cold_start(nott)
    ln, on, _ = ecef_to_llh(sol_n["pos"])
    _, tln, ton, thn = TRUTH["nottingham"]
    err_n = great_circle_km(sol_n["pos"], llh_to_ecef(tln, ton, thn))

    prior = llh_to_ecef(*SIGE_PRIOR)
    sol_s = solve_robust(sige, p0=prior, init=prior, earth_constraint=False)
    ls, os_, _ = ecef_to_llh(sol_s["pos"])
    _, tls, tos, ths = TRUTH["sige"]
    err_s = great_circle_km(sol_s["pos"], llh_to_ecef(tls, tos, ths))

    rows = list(SYNTHETIC) + [
        ("合成：3 星 + 地球约束（固有歧义）", 11450.0, "#bfa14a"),
        ("真实：Nottingham（1-bit，9 星）", err_n * 1e3, "#1f6feb"),
        ("真实：SiGe GN3S（8-bit，7 星）", err_s * 1e3, "#8e44ad"),
    ]
    labels = [r[0] for r in rows]
    vals = np.array([r[1] for r in rows], float)
    colors = [r[2] for r in rows]

    fig, ax = plt.subplots(figsize=(11.6, 5.6), dpi=130)
    y = np.arange(len(rows))[::-1]
    ax.barh(y, vals, color=colors, height=0.66)
    for yi, v in zip(y, vals):
        unit = f"{v:.3f} m" if v < 1000 else f"{v/1e3:.2f} km"
        ax.text(v * 1.25, yi, unit, va="center", fontsize=9.2)
    ax.axvline(1.0, color="#2e9e5b", lw=1.1, ls="--", alpha=0.8)
    ax.axvline(1e3, color="#8a94a6", lw=1.1, ls=":", alpha=0.9)
    ax.set_xscale("log")
    ax.set_xlim(3e-3, 3e5)
    ax.set_ylim(-0.60, len(rows) - 0.20)
    ax.text(1.06, len(rows) - 0.42, "1 m", fontsize=8.8, color="#2e9e5b",
            ha="left", va="center")
    ax.text(1.06e3, len(rows) - 0.42, "1 km", fontsize=8.8, color="#8a94a6",
            ha="left", va="center")
    ax.set_yticks(y); ax.set_yticklabels(labels, fontsize=9.5)
    ax.set_xlabel("位置误差 (m，对数轴)")
    ax.set_title("定位精度阶梯：合成真值（厘米级）→ 真实数据（米 / 公里级）",
                 fontsize=12)
    ax.grid(alpha=0.22, axis="x")
    ax.set_axisbelow(True)
    fig.tight_layout(rect=(0, 0.045, 1, 1))
    fig.text(0.012, 0.012,
             "合成数据取自 scripts/13_synthetic_4sat_pvt.py（可复现）；"
             "真实数据取自 scripts/14_nottingham_pvt.py、scripts/07_sige_pvt_test.py",
             fontsize=7.8, color="#6b7280")
    out = FIGDIR / "pvt_accuracy.png"
    fig.savefig(out); plt.close(fig)
    print(f"  已保存 {out.name}  Nottingham={err_n:.2f} km  SiGe={err_s:.2f} km")


def main() -> int:
    FIGDIR.mkdir(parents=True, exist_ok=True)
    nott = load("nottingham_pvt_measurements.pkl")
    sige = load("sige_pvt_measurements.pkl")
    cn0 = load_cn0()
    print(f"缓存：Nottingham {len(nott)} 星、SiGe {len(sige)} 星、"
          f"C/N₀ {len(cn0)} 通道")
    print("① 定位结果图"); fig_result(nott, cn0)
    print("② 冷启动代价地形"); fig_coldstart(nott)
    print("③ 精度阶梯"); fig_accuracy(nott, sige)
    print("\n全部完成。")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
