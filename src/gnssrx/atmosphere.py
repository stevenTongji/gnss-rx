"""大气传播延迟改正：对流层（Saastamoinen）+ 电离层（Klobuchar）。

对应 RTKLIB 的 `tropmodel()` / `ionmodel()`（rtkcmn.c）。本模块的公式与 RTKLIB
逐项对齐，便于交叉验证；差异只在接口形式（这里用角度制、返回米）。

为什么需要它：广播星历给出的是卫星在真空中的几何位置，而实测伪距里含有
    · 对流层延迟：天顶约 2.3 m，低仰角（10°）可达 15 m 以上，约 90% 为干分量
    · 电离层延迟：L1 上白天可达 5–30 m（夜间约 1–3 m），与 TEC 和仰角相关
不做改正时，最小二乘会把这部分**吸收进高程与钟差**——实测：注入 2.3 m 天顶
对流层延迟会让解算高程抬高约 11 m，而水平位置只移动约 2 m、残差几乎不变。
所以「残差很小」并不能证明大气延迟不存在，必须显式建模。

参考：
    · J. Saastamoinen, "Atmospheric Correction for the Troposphere and
      Stratosphere in Radio Ranging of Satellites", 1972.
    · J. A. Klobuchar, "Ionospheric Time-Delay Algorithm for Single-Frequency
      GPS Users", IEEE Trans. AES, 1987；IS-GPS-200 Table 20-IV。
    · RTKLIB `src/rtkcmn.c::tropmodel()` / `ionmodel()`。
"""
from __future__ import annotations

import math

import numpy as np

C_LIGHT = 299792458.0

# Klobuchar 系数（α 单位 s，β 单位 s/rad）的正常量级；真实值在 SF4 第 18 页
# （由 PRN 18 播发）。
_VALID_ALPHA = (-4.0e-8, 8.0e-8)      # 2^-30 s 量级
_VALID_BETA = (0.0, 1.5e5)            # s/rad 量级


def enu_azel(recv_ecef: np.ndarray, lat_deg: float, lon_deg: float,
             sat_ecef: np.ndarray) -> tuple[float, float, float]:
    """由接收机 ECEF 与卫星 ECEF 算 (方位角°, 仰角°, 距离 m)。"""
    la, lo = math.radians(lat_deg), math.radians(lon_deg)
    e = np.array([-math.sin(lo), math.cos(lo), 0.0])
    n = np.array([-math.sin(la) * math.cos(lo), -math.sin(la) * math.sin(lo),
                  math.cos(la)])
    u = np.array([math.cos(la) * math.cos(lo), math.cos(la) * math.sin(lo),
                  math.sin(la)])
    d = np.asarray(sat_ecef, float) - np.asarray(recv_ecef, float)
    rng = float(np.linalg.norm(d))
    de, dn, du = float(d @ e), float(d @ n), float(d @ u)
    az = math.degrees(math.atan2(de, dn)) % 360.0
    el = math.degrees(math.asin(max(-1.0, min(1.0, du / rng))))
    return az, el, rng


def troposphere_delay(el_deg: float, height_m: float = 0.0,
                      humidity: float = 0.7) -> float:
    """对流层天顶延迟映射到给定仰角的斜距延迟（米，正值）。

    与 RTKLIB `tropmodel()` 的 Saastamoinen 公式一致：

        tr = 0.002277 / cos z · (P + (1255/T + 0.05)·e − tan²z)

    其中标准大气 P(h)、T(h) 由接收机高程算出，水汽压 e 由相对湿度给出。
    """
    if el_deg <= 0.0:
        return 0.0
    h = min(10000.0, max(0.0, float(height_m)))
    P = 1013.25 * (1.0 - 2.2557e-5 * h) ** 5.2568          # 气压 (mbar)
    T = 15.0 - 6.5e-3 * h + 273.16                          # 温度 (K)
    e = 6.108 * humidity * math.exp((17.15 * T - 4684.0) / (T - 38.45))  # 水汽压
    z = math.radians(90.0 - el_deg)
    return float(0.002277 / math.cos(z)
                 * (P + (1255.0 / T + 0.05) * e - math.tan(z) ** 2))


def ionosphere_delay(alpha, beta, lat_deg: float, lon_deg: float,
                     az_deg: float, el_deg: float, gps_tow_s: float) -> float:
    """Klobuchar 电离层延迟（米，正值）。

    与 RTKLIB `ionmodel()` / IS-GPS-200 Table 20-IV 一致；角度内部换算为
    IS-GPS-200 惯用的「半周（semicircle）」单位。夜间（|X| ≥ 1.57）取 5e-9 s 底值。
    """
    a = np.asarray(alpha, float).ravel()[:4]
    b = np.asarray(beta, float).ravel()[:4]
    if el_deg <= 0.0:
        return 0.0
    # 角度 → 半周（semicircle）。IS-GPS-200 的 Klobuchar 公式全部用半周为单位：
    # 1 半周 = 180°。注意 RTKLIB 内部用弧度，故其代码写的是 `/PI`；本模块入参是
    # 角度制，必须 `/180`——写错会让 phi_i 远超 ±0.416 被钳位，算出荒唐的延迟。
    SC = 1.0 / 180.0
    phi_u = lat_deg * SC
    lam_u = lon_deg * SC
    el = el_deg * SC
    az = az_deg * SC

    psi = 0.0137 / (el + 0.11) - 0.022                       # 地心角
    phi_i = phi_u + psi * math.cos(az)
    if abs(phi_i) > 0.416:
        phi_i = 0.416 if phi_i > 0 else -0.416
    lam_i = lam_u + psi * math.sin(az) / math.cos(phi_i)
    phi_m = phi_i + 0.064 * math.cos((lam_i - 1.617) * math.pi)

    t = 4.32e4 * lam_i + gps_tow_s
    t = t % 86400.0

    amp = sum(float(a[k]) * phi_m ** k for k in range(4))     # 振幅项 (s)
    per = sum(float(b[k]) * phi_m ** k for k in range(4))     # 周期项 (s)
    if per < 72000.0:
        per = 72000.0
    x = 2.0 * math.pi * (t - 50400.0) / per
    if abs(x) < 1.57:
        F = 5.0e-9 + amp * (1.0 - x * x / 2.0 + x ** 4 / 24.0)
    else:
        F = 5.0e-9
    F *= 1.0 + 16.0 * (0.53 - el_deg * SC) ** 3               # 倾斜因子
    return float(F * C_LIGHT)


# 中纬度「标称」Klobuchar 系数（α0=2^-30·α_k 形式的值直接给出）。
# ⚠️ 这不是本数据集电文里的真值：SF4 第 18 页由 PRN 18 播发，本数据可见星中
#    不含 PRN 18，48 s 窗口也覆盖不到完整 12.5 min 的 SF4 页循环，故只能取标称值。
#    真实接收机从 SF4 第 18 页（或 RINEX 导航文件）读取该组系数。
NOMINAL_ALPHA = (7.6294e-9, 0.0, -5.9605e-8, 0.0)
NOMINAL_BETA = (1.1469e5, 0.0, -1.3107e5, 0.0)


def apply_atmospheric(measurements: list[dict], recv_ecef: np.ndarray,
                      lat_deg: float, lon_deg: float, height_m: float,
                      alpha=None, beta=None, gps_tow_s: float = 0.0,
                      use_iono: bool = True, use_tropo: bool = True,
                      humidity: float = 0.7) -> list[dict]:
    """把大气延迟从伪距里扣掉（伪距 = 几何距离 + 钟差 + 延迟）。

    返回新的测量列表，并记录 `_tropo_m` / `_iono_m` / `_el_deg` 便于审计。
    """
    if use_iono:
        alpha = NOMINAL_ALPHA if alpha is None else alpha
        beta = NOMINAL_BETA if beta is None else beta
    out = []
    for m in measurements:
        sat = np.asarray(m["sat_ecef"], float)
        az, el, _ = enu_azel(recv_ecef, lat_deg, lon_deg, sat)
        tropo = troposphere_delay(el, height_m, humidity) if use_tropo else 0.0
        iono = (ionosphere_delay(alpha, beta, lat_deg, lon_deg, az, el, gps_tow_s)
                if use_iono else 0.0)
        mm = dict(m)
        mm["pseudorange"] = m["pseudorange"] - tropo - iono
        mm["_tropo_m"] = tropo
        mm["_iono_m"] = iono
        mm["_el_deg"] = el
        mm["_az_deg"] = az
        out.append(mm)
    return out
