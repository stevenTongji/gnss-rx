"""GPS L1 单点定位（PVT）：伪距 + 迭代加权最小二乘。

================================================================================
定位模型
================================================================================
每颗卫星提供一个观测方程：

    ρ_i^meas = || sat_i(ECEF) − recv || + b + ε_i

其中：
    ρ_i^meas   伪距（由码相位 + TOW 得到，并已做卫星钟差改正，单位米）
    sat_i      卫星在信号发射时刻的 ECEF 坐标（由星历算出）
    recv       接收机位置 [x, y, z]（ECEF，米）
    b          接收机钟差等效距离（米，待估）
    ε_i        噪声/残余误差

未知量 [x, y, z, b] 共 4 个。4 颗及以上卫星可解（完整 SPS）；
3 颗卫星时欠定 1 个自由度，本模块额外加"地球表面约束"
（||recv|| ≈ R_e），无需预先知道接收机高度即可给出 3D 解。

迭代过程：以初猜展开线性化，加权最小二乘更新，直到收敛。
"""

from __future__ import annotations

import numpy as np

# ---- 物理 / WGS84 常数 ----
C_LIGHT = 299792458.0            # 光速 (m/s)
MU_EARTH = 3.986004418e14        # 地球引力常数 (m^3/s^2)
OMEGA_E = 7.2921150e-5           # 地球自转角速度 (rad/s)
R_EARTH = 6371000.0              # 平均地球半径 (m)，用于 3 星地球约束先验投影


def satellite_clock_correction(eph, t_sat: float) -> float:
    """卫星钟差等效距离（米）：δt_sat × c。

    包含：多项式（af0 + af1·Δt + af2·Δt²）与相对论项。
    t_sat 为发射时刻 GPS 秒（与电文 toc 同基准）。
    """
    dt = t_sat - eph.toc
    # 多项式
    dtc = eph.af0 + eph.af1 * dt + eph.af2 * dt * dt
    # 相对论修正：−2·(r·v)/c² 的近似形式，标准电文里用 e·sqrt(A)·sin(E)/c
    # 这里用 GPS 接口控制文档给出的等效式
    dt_rel = eph.relativistic_term(t_sat) if hasattr(eph, "relativistic_term") else 0.0
    return C_LIGHT * (dtc + dt_rel)


def _ecef_from_subframe_transmit(eph, t_sat: float) -> np.ndarray:
    """发射时刻的卫星 ECEF 坐标（含地球自转 Sagnac 修正）。

    1. 用 tk = t_sat − toe 算出"发射时"的卫星位置（惯性系视角下的轨道平面坐标）。
    2. 信号飞行 ~Δt ≈ ρ/c，期间地球转过 θ = OMEGA_E·Δt，把卫星坐标反方向
       旋转回去，得到发射瞬间卫星相对接收机处 ECEF 的真实位置。
    """
    tk = t_sat - eph.toes
    # 先猜一个大致距离做 Sagnac（首次迭代用标称 2.66e7 m）
    rho_approx = 2.66e7
    for _ in range(2):
        theta = OMEGA_E * (rho_approx / C_LIGHT)
        # 在发射时刻的位置（先按 tk 算，再绕 Z 轴转 -theta 回到接收 ECEF）
        pos_inertial = np.asarray(eph.position(tk), dtype=float)
        # 发射时卫星在 ECEF 下应再被地球自转"追上" theta，即接收机 ECEF 看过去要
        # 把卫星位置从 ECI 转到 ECEF：绕 Z 轴 +theta
        cos_t, sin_t = np.cos(theta), np.sin(theta)
        rot = np.array([[cos_t, sin_t, 0.0],
                        [-sin_t, cos_t, 0.0],
                        [0.0, 0.0, 1.0]])
        pos_ecef = rot @ pos_inertial
        rho_approx = float(np.linalg.norm(pos_ecef))
    return pos_ecef


def build_measurement(eph, prn: int, tow_seconds: float, t_user_seconds: float,
                      fs: float) -> dict:
    """由单颗卫星的（星历 + TOW + 接收机时刻）组装一个伪距观测。

    参数
    ----
    eph              Ephemeris 对象（已解出 SF1/2/3）
    prn              卫星号（Ephemeris 不存 prn，单独传）
    tow_seconds     read_tow 返回的 TOW（秒，已是 GPS 秒）
    t_user_seconds  该子帧起始处的接收机时刻（秒），必须来自"解缠后的累计码相位"，
                    否则只有毫秒级精度 → 伪距误差达数百 km。
    fs              采样率

    返回 dict：{prn, sat_ecef, pseudorange, sat_clock_m, t_sat}
    伪距的绝对时基偏置由接收机钟差 b 吸收，因此 TOW 只需同星期内一致即可。
    """
    chip_dur = 1.0 / 1.023e6
    # 发射时刻 ≈ 该子帧起始（TLM 前导）：TOW 指向下一子帧，故本子帧起始 = TOW − 6
    t_sat = tow_seconds - 6.0

    sat_ecef = _ecef_from_subframe_transmit(eph, t_sat)
    sat_clock_m = satellite_clock_correction(eph, t_sat)
    # 卫星钟差符号（关键）：c(t_user − t_sat) = ρ_true + c·b − c·δt_sat，
    # 故「已钟差改正的伪距」(= ρ_true + c·b) = c(t_user − t_sat) + c·δt_sat。
    # 符号必须取【加】——合成与真实数据（脚本 15/14）都已验证：取减号会引入
    # 随卫星变化的 ~0.5 ms 误差。绝对时基偏置 b 由最小二乘估计吸收。
    pseudorange = C_LIGHT * (t_user_seconds - t_sat) + sat_clock_m

    return {
        "prn": prn,
        "sat_ecef": sat_ecef,
        "pseudorange": float(pseudorange),
        "sat_clock_m": float(sat_clock_m),
        "t_sat": float(t_sat),
    }


def resolve_ms_ambiguity(measurements: list[dict],
                         p0: np.ndarray,
                         c: float = C_LIGHT) -> list[dict]:
    """消除单历元伪距的整数毫秒模糊度（integer-millisecond ambiguity）。

    单历元只能测到亚毫秒码相位；传播时延约 67–94 ms（≈ 2–2.8 万 km），其整数
    个 1 ms 码周期是模糊的。由于各星子帧相对 GPS 时对齐，单靠本地块计数会把
    每颗星约 8000 km 的几何距离差压缩到亚毫秒残差里，导致定位错误。

    做法：用先验位置 p0 估计 ||sat_i − p0||，给每颗星选一个整数 N_i 使伪距与先验
    距离一致：N_i = round((||sat_i−p0|| − ρ_i) / (c·1ms))，再把 N_i·c·1ms 加回伪距。
    p0 只需大致在地球表面（几万 km 误差只会让 N_i 差 ≤ 1，迭代即可收敛）。
    """
    out = []
    for m in measurements:
        sat = np.asarray(m["sat_ecef"], dtype=float)
        r_prior = float(np.linalg.norm(sat - np.asarray(p0, dtype=float)))
        n = round((r_prior - m["pseudorange"]) / (c * 1e-3))
        mm = dict(m)
        mm["pseudorange"] = m["pseudorange"] + n * c * 1e-3
        mm["_ambig_n"] = int(n)
        out.append(mm)
    return out


def solve_robust(measurements: list[dict],
                init: np.ndarray | None = None,
                earth_constraint: bool = False,
                p0: np.ndarray | None = None,
                max_iter: int = 20,
                ambig_iter: int = 4) -> dict:
    """先消除整数毫秒模糊度再迭代最小二乘定位（3 星静态定位必需）。

    先以卫星质心投影到地球表面作为先验 p0，解出位置后再用新位置回代重解模糊度，
    迭代 ambig_iter 次即可收敛到正确整数毫秒。
    """
    if p0 is None:
        cvec = np.mean([np.asarray(m["sat_ecef"], dtype=float) for m in measurements],
                       axis=0)
        nrm = float(np.linalg.norm(cvec))
        p0 = (cvec / nrm * R_EARTH) if nrm > 1.0 else np.array([0.0, 0.0, R_EARTH])
    p = np.asarray(p0, dtype=float)
    sol = None
    for _ in range(ambig_iter):
        meas = resolve_ms_ambiguity(measurements, p)
        sol = solve(meas, init=init if init is not None else p,
                    earth_constraint=earth_constraint, max_iter=max_iter)
        if not sol["converged"]:
            break
        p = sol["pos"]
    assert sol is not None
    sol["ambig_applied"] = True
    return sol


def _geo_matrix(sats: list[np.ndarray], recv: np.ndarray) -> np.ndarray:
    """几何矩阵 G（每颗星一行 [ux, uy, uz, 1]）。

    ρ_i = ||sat_i − recv|| + b，对 recv 的雅可比 = (recv − sat_i)/||·|| = u，
    对 b 的雅可比 = 1。故 G 行 = [u, 1]，线性化方程 G·Δx = 残差。
    """
    n = len(sats)
    G = np.zeros((n, 4))
    for i, s in enumerate(sats):
        rng = float(np.linalg.norm(s - recv))
        if rng < 1.0:
            rng = 1.0
        u = (recv - s) / rng                 # 单位矢量 recv→sat
        G[i, 0:3] = u
        G[i, 3] = 1.0
    return G


def _solve_3sat_spherical(sats: list[np.ndarray], rho: np.ndarray,
                          init: np.ndarray, max_iter: int = 30) -> dict:
    """3 星 + 地球表面约束：把接收机限制在球面上，未知量 (lat, lon, b) 共 3 个，
    3 个距离方程正好良态（避免 4 未知欠定导致的零空间发散）。

    高度假设为 0（平均地球半径处）；对屋顶天线等会有数十~百米水平偏差，
    属正常。返回与 solve() 同结构的 dict。
    """
    # 初猜投影到球面
    n0 = float(np.linalg.norm(init))
    if n0 < 1.0:
        lat, lon = 0.0, 0.0
    else:
        lat = np.arcsin(init[2] / n0)
        lon = np.arctan2(init[1], init[0])
    b = 0.0

    def recv_at(la, lo):
        return R_EARTH * np.array([np.cos(la) * np.cos(lo),
                                    np.cos(la) * np.sin(lo),
                                    np.sin(la)])

    converged = False
    for it in range(max_iter):
        recv = recv_at(lat, lon)
        dr_dlat = R_EARTH * np.array([-np.sin(lat) * np.cos(lon),
                                      -np.sin(lat) * np.sin(lon),
                                      np.cos(lat)])
        dr_dlon = R_EARTH * np.array([-np.cos(lat) * np.sin(lon),
                                      np.cos(lat) * np.cos(lon),
                                      0.0])
        J = np.zeros((3, 3))          # 3 方程 × 3 未知 (lat, lon, b)
        r = np.zeros(3)
        for i, s in enumerate(sats):
            rng = float(np.linalg.norm(s - recv))
            u = (recv - s) / rng                  # 单位矢量 recv→sat
            J[i, 0] = u @ dr_dlat                 # dρ/dlat
            J[i, 1] = u @ dr_dlon                 # dρ/dlon
            J[i, 2] = 1.0                         # dρ/db
            r[i] = rho[i] - (rng + b)
        try:
            d = np.linalg.solve(J.T @ J, J.T @ r)
        except np.linalg.LinAlgError:
            break
        # r = rho − rng − b，对参数的雅可比为 −J，故牛顿步 参数 += d
        lat += d[0]
        lon += d[1]
        b += d[2]
        if np.max(np.abs(d)) < 1e-9:
            converged = True
            break

    recv = recv_at(lat, lon)
    pred = np.array([float(np.linalg.norm(s - recv)) for s in sats])
    residuals = rho - (pred + b)
    return {
        "pos": recv,
        "clock_bias": b,
        "residuals": residuals,
        "n_sat": 3,
        "iter": it + 1,
        "converged": converged,
        "gdop": float("nan"),
    }


def solve(measurements: list[dict],
          init: np.ndarray | None = None,
          earth_constraint: bool = False,
          max_iter: int = 20) -> dict:
    """迭代加权最小二乘定位。

    measurements : build_measurement 返回的列表，每颗星一个。
    init          : [x,y,z] 初猜（米，ECEF）。
    earth_constraint : True 且恰好 3 颗星时，走球面参数化求解（接收机固定在地球表面）。

    返回 dict：{pos, clock_bias, residuals, n_sat, iter, converged, gdop}
    """
    sats = [np.asarray(m["sat_ecef"], dtype=float) for m in measurements]
    rho = np.array([m["pseudorange"] for m in measurements], dtype=float)
    n = len(sats)

    if n == 3 and earth_constraint:
        if init is None:
            c = np.mean(sats, axis=0)
            nrm = float(np.linalg.norm(c))
            init = (c / nrm * R_EARTH) if nrm > 1.0 else np.array([0.0, 0.0, R_EARTH])
        return _solve_3sat_spherical(sats, rho, init, max_iter=max_iter)

    if init is None:
        c = np.mean(sats, axis=0)
        nrm = float(np.linalg.norm(c))
        init = (c / nrm * R_EARTH) if nrm > 1.0 else np.array([0.0, 0.0, R_EARTH])
    recv = init.copy()
    b = 0.0

    w_earth = 1.0 if (earth_constraint and n >= 4) else 0.0

    converged = False
    for it in range(max_iter):
        G = _geo_matrix(sats, recv)
        pred = np.array([float(np.linalg.norm(s - recv)) for s in sats])
        resid = rho - (pred + b)

        rows = [G]
        rside = [resid]
        weights = [np.ones(n)]
        if w_earth > 0.0:
            rng = float(np.linalg.norm(recv))
            if rng < 1.0:
                rng = 1.0
            u_radial = recv / rng
            g_earth = np.array([u_radial[0], u_radial[1], u_radial[2], 0.0])
            rows.append(g_earth.reshape(1, 4))
            rside.append(np.array([rng - R_EARTH]))
            weights.append(np.array([w_earth]))

        Gf = np.vstack(rows)
        rf = np.concatenate(rside)
        W = np.diag(np.concatenate(weights))

        GTW = Gf.T @ W
        lhs = GTW @ Gf
        rhs = GTW @ rf
        try:
            dx = np.linalg.solve(lhs, rhs)
        except np.linalg.LinAlgError:
            break
        recv = recv + dx[0:3]
        b = b + dx[3]

        if np.max(np.abs(dx)) < 1e-4:
            converged = True
            break

    pred = np.array([float(np.linalg.norm(s - recv)) for s in sats])
    residuals = rho - (pred + b)
    G = _geo_matrix(sats, recv)
    Q = np.linalg.inv(G.T @ G) if G.shape[0] >= 4 else None
    gdop = float(np.sqrt(np.trace(Q))) if Q is not None else float("nan")

    return {
        "pos": recv,
        "clock_bias": b,
        "residuals": residuals,
        "n_sat": n,
        "iter": it + 1,
        "converged": converged,
        "gdop": gdop,
    }


def ecef_to_llh(xyz: np.ndarray) -> tuple[float, float, float]:
    """ECEF (m) → [lat_deg, lon_deg, height_m]（WGS84 近似，用平均半径）。"""
    x, y, z = float(xyz[0]), float(xyz[1]), float(xyz[2])
    lon = np.degrees(np.arctan2(y, x))
    p = np.hypot(x, y)
    lat = np.degrees(np.arctan2(z, p * (1.0 - 1.0 / 298.257223563)))
    # 几次迭代精化
    a = 6378137.0
    f = 1.0 / 298.257223563
    e2 = f * (2.0 - f)
    for _ in range(5):
        sinlat = np.sin(np.radians(lat))
        n = a / np.sqrt(1.0 - e2 * sinlat * sinlat)
        lat = np.degrees(np.arctan2(z + e2 * n * sinlat, p))
    sinlat = np.sin(np.radians(lat))
    n = a / np.sqrt(1.0 - e2 * sinlat * sinlat)
    height = p / np.cos(np.radians(lat)) - n
    return lat, lon, height


def llh_to_ecef(lat_deg: float, lon_deg: float, height_m: float = 0.0) -> np.ndarray:
    """[lat_deg, lon_deg, height_m] → ECEF (m)，WGS84。"""
    a = 6378137.0
    f = 1.0 / 298.257223563
    e2 = f * (2.0 - f)
    lat = np.radians(lat_deg)
    lon = np.radians(lon_deg)
    n = a / np.sqrt(1.0 - e2 * np.sin(lat) ** 2)
    return np.array([(n + height_m) * np.cos(lat) * np.cos(lon),
                     (n + height_m) * np.cos(lat) * np.sin(lon),
                     (n * (1.0 - e2) + height_m) * np.sin(lat)])


def _grid_ecef(lats: np.ndarray, lons: np.ndarray) -> np.ndarray:
    """经纬网格 → ECEF（h=0），返回 (nlat, nlon, 3)。矢量化。"""
    a = 6378137.0
    f = 1.0 / 298.257223563
    e2 = f * (2.0 - f)
    la = np.radians(lats)[:, None]
    lo = np.radians(lons)[None, :]
    n = a / np.sqrt(1.0 - e2 * np.sin(la) ** 2)
    x = n * np.cos(la) * np.cos(lo)
    y = n * np.cos(la) * np.sin(lo)
    z = n * (1.0 - e2) * np.sin(la) + 0.0 * lo        # 广播到 (nlat, nlon)
    return np.stack([x, y, z], axis=-1)


def solve_cold_start(measurements: list[dict], *,
                     coarse_step_deg: float = 1.0,
                     fine_step_deg: float = 0.2,
                     refine_deg: float = 1.5) -> dict:
    """冷启动定位：无需任何先验位置，全球网格搜索 + 局部精化 + 最小二乘。

    代价函数：对候选点 p 先用它做先验消整数毫秒模糊（逐星取整），再算
    「消模糊后伪距 − 几何距离」的**散布**。只有 p 落在真位置 ~150 km 内时，
    逐星取整才全部正确、散布才趋近 0；否则至少一颗星差一个 1 ms（≈300 km），散布骤增。
    因此该代价在真位置处有尖锐极小值。粗搜用 1° 网格保证不越过这个 ~150 km 的阱，
    再在最优附近细搜，最后交给 solve_robust 精化。

    适用于 ≥4 星（完整 SPS）。返回与 solve 相同结构的字典。
    """
    cms = C_LIGHT * 1e-3
    sats = np.array([np.asarray(m["sat_ecef"], dtype=float) for m in measurements])
    raw = np.array([m["pseudorange"] for m in measurements], dtype=float)

    def cost_grid(lats: np.ndarray, lons: np.ndarray) -> np.ndarray:
        g = _grid_ecef(lats, lons)                                 # (nlat,nlon,3)
        d = np.linalg.norm(g[:, :, None, :] - sats[None, None, :, :], axis=-1)  # (…,nsat)
        n = np.round((d - raw) / cms)
        return (raw[None, None, :] + n * cms - d).std(axis=-1)      # (nlat,nlon)

    def best_of(lats, lons):
        c = cost_grid(np.asarray(lats, float), np.asarray(lons, float))
        i = int(np.argmin(c))
        ilat, ilon = np.unravel_index(i, c.shape)
        return float(c[ilat, ilon]), float(np.asarray(lats)[ilat]), float(np.asarray(lons)[ilon])

    # ① 全球粗搜（有人居住的纬度带，1°）
    c0, b_lat, b_lon = best_of(np.arange(-80.0, 80.0 + 1e-9, coarse_step_deg),
                               np.arange(-180.0, 180.0, coarse_step_deg))
    # ② 最优附近细搜
    lats = np.arange(b_lat - refine_deg, b_lat + refine_deg + 1e-9, fine_step_deg)
    lons = np.arange(b_lon - refine_deg, b_lon + refine_deg + 1e-9, fine_step_deg)
    lats = lats[(-90.0 <= lats) & (lats <= 90.0)]
    c1, b_lat, b_lon = best_of(lats, lons)

    p0 = llh_to_ecef(b_lat, b_lon, 0.0)
    earth = (len(measurements) == 3)
    sol = solve_robust(measurements, p0=p0, init=p0, earth_constraint=earth)
    sol["cold_start_llh"] = (b_lat, b_lon)
    sol["cold_start_cost"] = min(c0, c1)
    return sol
