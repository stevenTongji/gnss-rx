#!/usr/bin/env python3
"""验证整毫秒模糊度消除（resolve_ms_ambiguity + solve_robust）在大整数模糊度下仍正确。

方法：用真实卫星位置（从缓存测量里取），把接收机"放"在已知 ECEF（Nottingham），
算出真实距离 R_i，构造带大整数毫秒模糊度的伪距（模拟 GPS 时基不明 + 码相位只测亚毫秒），
再用 solve_robust 从已知的粗略先验（Nottingham 自身 / 加噪声）恢复位置，比对真值。
"""
from __future__ import annotations
import pickle, sys, math
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from gnssrx.pvt import build_measurement, ecef_to_llh, resolve_ms_ambiguity, solve_robust

C = 299792458.0
CACHE = ROOT / "data" / "processed" / "sige_pvt_measurements.pkl"

def llh_str(xyz):
    lat, lon, h = ecef_to_llh(np.asarray(xyz))
    return f"lat={lat:.4f}° lon={lon:.4f}° h={h:.1f} m"

# 已知真值：Nottingham
lat0, lon0 = math.radians(52.95), math.radians(-1.15)
a = 6378137.0; f = 1/298.257223563; e2 = f*(2-f)
N0 = a/math.sqrt(1-e2*math.sin(lat0)**2)
truth = np.array([(N0)*math.cos(lat0)*math.cos(lon0),
                  (N0)*math.cos(lat0)*math.sin(lon0),
                  (N0*(1-e2))*math.sin(lat0)])
print("真值 ECEF:", truth, "\n真值", llh_str(truth))

with open(CACHE, "rb") as fh:
    meas = pickle.load(fh)
sats = [np.asarray(m["sat_ecef"], dtype=float) for m in meas]
print(f"使用 {len(sats)} 颗真实卫星位置做仿真\n")

# 构造"测量伪距"：真实距离 + 钟差 - 大整数毫秒模糊度（模拟 GPS 时基未知）
b_true = 0.0      # 钟差设为 0：真实接收机的钟差由解算的 b 连续吸收，不量化进整数毫秒
K_ambig = 434579  # 约 434579 s 的 GPS 时基偏移 → 整数毫秒模糊度（干净整数 ms）
print(f"注入模糊度 K={K_ambig} ms (= {K_ambig*C*1e-3/1e6:.1f} Mm)，钟差 b={b_true:.0f} m")
sim_meas = []
for m, s in zip(meas, sats):
    R = float(np.linalg.norm(s - truth))
    # 真实伪距（含钟差）
    rho_true = R + b_true
    # 接收机"只测到" rho_true 减去一个大整数毫秒（模糊度），剩下的由码相位亚毫秒给出
    rho_meas = rho_true - K_ambig * C * 1e-3
    mm = dict(m)
    mm["pseudorange"] = rho_meas
    mm["sat_ecef"] = s
    sim_meas.append(mm)

# 情形 1：先验 = 真值 Nottingham（理想）
print("\n=== 情形1：先验=Nottingham 真值 ===")
sol1 = solve_robust(sim_meas, earth_constraint=(len(sim_meas) == 3), p0=truth.copy())
print("  解", llh_str(sol1["pos"]), "收敛=", sol1["converged"], "残差", [f"{r:+.2f}" for r in sol1["residuals"]])
err1 = np.linalg.norm(sol1["pos"] - truth)
print(f"  与真值误差 = {err1:.2f} m")

# 情形 2：先验 = 真值 + 5 km 扰动（现实中的粗略先验）
p0b = truth + np.array([3000.0, -4000.0, 2000.0])
print("\n=== 情形2：先验=真值+~5km 扰动 ===")
sol2 = solve_robust(sim_meas, earth_constraint=(len(sim_meas) == 3), p0=p0b)
print("  解", llh_str(sol2["pos"]), "收敛=", sol2["converged"])
err2 = np.linalg.norm(sol2["pos"] - truth)
print(f"  与真值误差 = {err2:.2f} m")

# 情形 3：先验 = 卫星质心投影（最粗糙，验证是否仍可用）
cvec = np.mean(sats, axis=0); cvec = cvec/np.linalg.norm(cvec)*6371000.0
print("\n=== 情形3：先验=卫星质心投影（最粗糙）===")
sol3 = solve_robust(sim_meas, earth_constraint=(len(sim_meas) == 3), p0=cvec)
print("  解", llh_str(sol3["pos"]), "收敛=", sol3["converged"])
err3 = np.linalg.norm(sol3["pos"] - truth)
print(f"  与真值误差 = {err3:.2f} m")

print("\n结论：情形1/2 应误差 < 100 m 即证明模糊度消除正确；情形3 可能偏差大（先验太粗）。")
