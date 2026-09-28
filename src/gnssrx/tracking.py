"""GPS L1 跟踪环路：码环（DLL）+ 载波环（PLL）。

================================================================================
一、捕获和跟踪的分工
================================================================================
捕获只回答"天上有哪几颗星、粗码相位、粗多普勒"，精度是 1 个码片 / 几百赫兹，
而且只用了 1 ms 数据。跟踪要做的才是真正的解调：
    从这一毫秒到下一毫秒，持续咬住码相位和载波相位，
    每毫秒吐出一组 I/Q，后面位同步、帧同步、伪距全都基于这个输出。

================================================================================
二、信号模型（先统一记号，否则后面全是坑）
================================================================================
第 n 个采样点（中频、实数）：

    s[n] = A · d[n] · c( χ_sig(n) ) · cos( 2π(f_IF + f_d)·n/fs + φ_s )

    A     幅度
    d[n]  导航电文比特，±1，每 20 ms 可能翻转一次
    c(·)  C/A 码，周期 1023 chip
    χ_sig 信号的码相位（单位：chip），推进速率 = 码速率 f_code
    f_d   多普勒
    φ_s   载波初相

本地复现两份东西：

    本地码      c( χ_loc(n) )，χ_loc 推进速率 = 码 NCO 频率 f_code_loc
    本地载波    exp( -j(2π f_carr·n/fs + φ_loc) )，f_carr = f_IF + f_d 估计值

定义两个误差量（这是全模块的符号基准）：

    Δ  = χ_loc - χ_sig       码相位误差，单位 chip（>0 表示本地码跑在前面）
    φ  = φ_s  - φ_loc        载波相位误差，单位 rad

下变频 + 与本地码相关 + 1 ms 积分后（记 P 为即时支路）：

    I_P = (A·N/2)·d·R(Δ)·cos φ
    Q_P = (A·N/2)·d·R(Δ)·sin φ          ← 注意符号，见下面"Q 为什么要取负"

    N = 一个积分周期内的采样点数，R(·) 是码自相关（三角峰，R(0)=1）

================================================================================
三、Q 支路为什么要取负号
================================================================================
用复指数写最干净。本振取 exp(-j θ_loc)，则

    s[n]·exp(-jθ_loc) = (A d c/2)·[ e^{j(θ_sig-θ_loc)} + 高频项 ]
                      → 低通后 = (A d c/2)·[ cos φ + j sin φ ]

而 exp(-jθ) = cos θ - j sin θ，所以

    Re(基带) = s[n]·cos θ_loc   →  I ∝ cos φ
    Im(基带) = -s[n]·sin θ_loc  →  Q ∝ sin φ      ← 负号就是这么来的

不取这个负号，鉴相器会变成 -φ，环路变成正反馈，一秒之内必失锁。

================================================================================
四、鉴相器
================================================================================
【码环 DLL】归一化、非相干、超前减滞后：

    E = |IE + jQE|,  L = |IL + jQL|        （取模 → 不受数据比特翻转影响）
    e_code = -(E - L) / (2(E + L))         单位：chip

    推导：三路相关器间距 ±d/2 chip（本实现 d = 1 chip，即 E 在 +0.5、L 在 -0.5）。
    三角自相关下，|Δ| < 1 - d/2 时
        E = R(Δ + 0.5) = 0.5 - Δ
        L = R(Δ - 0.5) = 0.5 + Δ
        ⇒ (E - L)/(E + L) = -2Δ   ⇒   Δ = -(E - L)/(2(E + L))
    鉴相器增益 = -2 / chip。归一化是为了让增益不随信号强度变化，
    否则强星和弱星要分别调环路参数。

【载波环 PLL】Costas 鉴相器：

    e_carr = atan( Q_P / I_P )             单位：rad

    为什么是 atan 而不是 atan2？
    atan2 会给 (-π, π] 的完整象限信息，但导航电文每 20 ms 把 d 翻一次符号，
    I_P 和 Q_P 会同时变号 —— atan2 的结果会跳变 π，环路直接被踢飞。
    而 atan(Q/I) 在 (I,Q) → (-I,-Q) 下不变（这就是 Costas 的本质），
    代价是留下 180° 相位模糊，交给后面的帧同步去解。

================================================================================
五、环路滤波器（二阶 PI）
================================================================================
模拟域：F(s) = (1 + τ2·s) / (τ1·s)
数字化（冲激不变法，采样周期 T = PDI）：

    v[n] = v[n-1] + (T/τ1)·e[n] + (τ2/τ1)·(e[n] - e[n-1])

系数（Kaplan & Hegarty《Understanding GPS/GNSS》）：

    ω_n = 8ζ·B_n / (4ζ² + 1)
    τ2  = 2ζ / ω_n
    τ1  = K / ω_n²

K 是环路总增益（鉴相器增益 × NCO 增益）。**本实现把鉴相器输出归一化成
被控量本身**（码：chip；载波：rad），于是 K = 1、τ1 = 1/ω_n²。
好处是 B_n 就是字面意义上的噪声带宽，换数据集/换采样率不用重新标定增益。

闭环特征方程推导（检验用）：
    设 de/dt = -v，v = F(s)·e  ⇒  s + F(s) = 0
    ⇒ τ1 s² + K τ2 s + K = 0 ⇒ s² + (K τ2/τ1)s + K/τ1 = 0
    与 s² + 2ζω_n s + ω_n² 对比 ⇒ τ2 = 2ζ/ω_n，τ1 = K/ω_n²  ✓ 自洽

二阶环路（II 型）对**频率阶跃**零稳态误差 —— 这正是它必须至少二阶的原因：
多普勒和码多普勒对环路而言都是频率阶跃，一阶环会留一个恒定的相位偏差。

================================================================================
六、符号约定（最易错，逐条列出）
================================================================================
码环：  s = Δ，ds/dt = f_code_loc - f_code_sig
        反馈取 f_code_loc = CODE_RATE - v_code   ⇒ ds/dt = -v_code - 常数项
        Δ > 0（本地码超前）→ e > 0 → v > 0 → 码速率降低 → 本地码慢下来 ✓

载波环：dφ/dt = 2π(f_sig - f_loc) = -2π·u，u = f_loc - f_sig = v_carr/(2π)
        ⇒ dφ/dt = -v_carr ✓
        ⇒ f_carr = f_IF + f_d(捕获值) + v_carr/(2π)

初始码相位：本地码相位 τ（采样点）与捕获得到的码延迟 D 的关系是
        τ_0 = (-D) mod spc
    因为捕获的相关峰下标 k 就是码延迟 D（见 acquisition.py 的推导），
    而本地码要"提前" D 个采样点才能和信号对齐。

================================================================================
七、C/N0 估计（用来判断锁没锁上）
================================================================================
锁定后 I_P 承载信号、Q_P 只剩噪声。
设噪声方差 σ²，N 点积分：
    |I_P| = A·N/2        （信号）
    Var(Q_P) = σ²·N/2    （噪声，因 sin 与 cos 各分一半功率）
    C = A²/2，N0 = 2σ²/fs
消掉 A 和 σ 得：
        C/N0 = |I_P|² / (2·T·Var(Q_P))
这个式子可以直接用来校验仿真：合成时给的 C/N0 是 45 dB-Hz，
估计出来应该就在 45 附近 —— 这是验证整个跟踪链路增益是否正确的独立手段。
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .ca_code import CODE_LENGTH, CODE_RATE_HZ, L1_HZ, generate_ca

TWO_PI = 2.0 * np.pi


@dataclass
class TrackingConfig:
    """跟踪环路参数。默认值是 GPS L1 C/A、1 ms 积分下的常用取值。"""
    dll_bandwidth_hz: float = 2.0          # 码环噪声带宽。窄 → 平滑但跟不上动态
    dll_damping: float = 0.7071            # 阻尼，1/√2 是阶跃响应无超调的经典值
    pll_bandwidth_hz: float = 25.0         # 载波环噪声带宽。宽 → 跟得上动态但噪声大
    pll_damping: float = 0.7071
    correlator_spacing_chips: float = 1.0  # E 与 L 的总间距（d）。窄相关可减小多径误差
    pdi_ms: int = 1                        # 预检测积分时间（ms）
    carrier_aiding: bool = True            # 用载波多普勒辅助码环，见下面说明

    # ------------------------------------------------------------------
    # 为什么必须有载波辅助（carrier aiding）
    # ------------------------------------------------------------------
    # 码多普勒 f_code = 1.023 MHz × (1 + f_d / f_L1)。3000 Hz 的多普勒
    # 对应码速率偏差约 1.95 chip/s，即码相位每秒漂 2 个码片。
    # 而 2 Hz 码环的自然频率 ω_n = 8ζB_n/(4ζ²+1) ≈ 3.8 rad/s，
    # 对频率阶跃的瞬态峰值误差约 Δf/ω_n ≈ 0.5 chip —— 正好撞在
    # 超前减滞后鉴相器线性区（±(1 - d/2) = ±0.5 chip）的边界上，一冲就滑出去。
    #
    # 载波环测出的多普勒精度远高于码环，而且两者同源：
    #     f_code = f_chip_nominal × (1 + f_d / f_L1)
    # 所以直接用载波多普勒去推码 NCO，码环只需处理残差（近零），
    # 瞬态峰值和稳态误差都掉一个数量级。这是所有真实接收机的标准做法。
    # ------------------------------------------------------------------


class LoopFilter:
    """二阶比例-积分环路滤波器，见模块文档第五节。"""

    def __init__(self, bandwidth_hz: float, damping: float, pdi_s: float):
        self.bandwidth_hz = float(bandwidth_hz)
        self.damping = float(damping)
        self.pdi = float(pdi_s)

        self.omega_n = bandwidth_hz * 8.0 * damping / (4.0 * damping * damping + 1.0)
        self.tau1 = 1.0 / (self.omega_n ** 2)      # K = 1（鉴相器输出已归一化）
        self.tau2 = 2.0 * damping / self.omega_n

        self._prev_error = 0.0
        self.output = 0.0

    def reset(self) -> None:
        self._prev_error = 0.0
        self.output = 0.0

    def step(self, error: float) -> float:
        """喂入本次鉴相器输出，返回累积的 NCO 控制量。"""
        self.output += (self.pdi / self.tau1) * error \
                       + (self.tau2 / self.tau1) * (error - self._prev_error)
        self._prev_error = error
        return self.output

    def __repr__(self) -> str:
        return (f"LoopFilter(Bn={self.bandwidth_hz} Hz, zeta={self.damping}, "
                f"wn={self.omega_n:.2f} rad/s, tau1={self.tau1:.3e}, tau2={self.tau2:.3e})")


def _code_replica(code: np.ndarray, tau_samples: float, n_samples: int, fs: float) -> np.ndarray:
    """按本地码相位 τ（采样点，可为小数）重采样出一段本地码。"""
    idx = np.floor((np.arange(n_samples) + tau_samples) * (CODE_RATE_HZ / fs)).astype(np.int64)
    return code[idx % CODE_LENGTH]


@dataclass
class TrackingResult:
    """单个通道的跟踪输出，每个数组长度 = 积分块数。"""
    prn: int
    ie: np.ndarray
    qe: np.ndarray
    ip: np.ndarray
    qp: np.ndarray
    il: np.ndarray
    ql: np.ndarray
    dll_error_chips: np.ndarray
    pll_error_rad: np.ndarray
    code_freq_hz: np.ndarray
    doppler_hz: np.ndarray
    code_phase_samples: np.ndarray
    config: TrackingConfig = field(default_factory=TrackingConfig)

    @property
    def n_blocks(self) -> int:
        return len(self.ip)

    def prompt_amp(self) -> np.ndarray:
        return np.abs(self.ip)

    def cn0_dbhz(self, pdi_s: float, start: int = 0) -> float:
        """C/N0 估计，见模块文档第七节。

        start 用来跳过开头的牵引瞬态 —— 那段时间 |I_P| 还没起来，
        算进去会低估好几个 dB。
        """
        amp = float(np.mean(np.abs(self.ip[start:])))
        q_var = float(np.var(self.qp[start:]))
        if q_var <= 0.0 or amp <= 0.0:
            return float("nan")
        return float(10.0 * np.log10(amp ** 2 / (2.0 * pdi_s * q_var)))

    def phase_error_deg(self) -> np.ndarray:
        return np.degrees(self.pll_error_rad)


def track_channel(
    data: np.ndarray,
    fs: float,
    f_if: float,
    prn: int,
    *,
    init_code_phase_samples: int,
    init_doppler_hz: float,
    n_ms: int,
    cfg: TrackingConfig | None = None,
) -> TrackingResult:
    """对一个卫星通道做闭环跟踪。

    参数
    ----
    data                    中频采样（实数）
    fs, f_if                采样率与中频，必须和捕获时用的一致
    prn                     卫星号
    init_code_phase_samples 捕获给出的码延迟（采样点）
    init_doppler_hz         捕获给出的多普勒（Hz）
    n_ms                    跟踪时长
    cfg                     环路参数
    """
    cfg = cfg or TrackingConfig()

    spc = int(round(fs / 1e3))                  # 1 ms 的采样点数
    block = spc * cfg.pdi_ms                    # 一个积分块的采样点数
    n_blocks = int(n_ms // cfg.pdi_ms)
    pdi_s = cfg.pdi_ms * 1e-3

    if data.size < block * n_blocks:
        raise ValueError(f"数据不足：需要 {block * n_blocks} 点，实际 {data.size}")

    code = generate_ca(prn).astype(np.float64)
    # 关键：E 在 +d/2、L 在 -d/2，d 是两者的总间距。
    # 若误用 ±d（即间距 2d），d=1 时 E/L 落在 ±1 码片 —— 码自相关在那里恰好为 0，
    # 鉴相器在锁定点附近只剩噪声、没有增益，环路必然跑飞。
    spacing_samples = 0.5 * cfg.correlator_spacing_chips * fs / CODE_RATE_HZ
    samples_per_chip = fs / CODE_RATE_HZ

    dll = LoopFilter(cfg.dll_bandwidth_hz, cfg.dll_damping, pdi_s)
    pll = LoopFilter(cfg.pll_bandwidth_hz, cfg.pll_damping, pdi_s)

    # 本地码相位（采样点）。符号约定见模块文档第六节：τ0 = (-D) mod spc
    tau = float((-init_code_phase_samples) % spc)
    code_freq = CODE_RATE_HZ
    carr_freq = f_if + init_doppler_hz
    carr_phase = 0.0

    k = np.arange(block)
    ts = 1.0 / fs

    ie = np.zeros(n_blocks); qe = np.zeros(n_blocks)
    ip = np.zeros(n_blocks); qp = np.zeros(n_blocks)
    il = np.zeros(n_blocks); ql = np.zeros(n_blocks)
    dll_err = np.zeros(n_blocks); pll_err = np.zeros(n_blocks)
    code_freqs = np.zeros(n_blocks); dopplers = np.zeros(n_blocks)
    taus = np.zeros(n_blocks)

    for m in range(n_blocks):
        seg = data[m * block:(m + 1) * block]

        # ---- 载波剥离：生成本地 I/Q 本振并混频 ----
        ph = carr_phase + TWO_PI * carr_freq * ts * k
        mixed_i = seg * np.cos(ph)
        mixed_q = seg * np.sin(ph)

        # ---- 三路相关器：超前 / 即时 / 滞后 ----
        c_e = _code_replica(code, tau + spacing_samples, block, fs)
        c_p = _code_replica(code, tau, block, fs)
        c_l = _code_replica(code, tau - spacing_samples, block, fs)

        ie[m] = mixed_i @ c_e
        qe[m] = -(mixed_q @ c_e)          # 负号见文档第三节
        ip[m] = mixed_i @ c_p
        qp[m] = -(mixed_q @ c_p)
        il[m] = mixed_i @ c_l
        ql[m] = -(mixed_q @ c_l)

        # ---- 鉴相器 ----
        e_mag = float(np.hypot(ie[m], qe[m]))
        l_mag = float(np.hypot(il[m], ql[m]))
        denom = e_mag + l_mag
        dll_err[m] = -(e_mag - l_mag) / (2.0 * denom) if denom > 0.0 else 0.0

        # Costas：atan(Q/I)，对 180° 翻转不敏感。I 接近 0 时做数值保护
        i_val = ip[m]
        if abs(i_val) < 1e-9:
            i_val = 1e-9 if i_val >= 0.0 else -1e-9
        pll_err[m] = float(np.arctan(qp[m] / i_val))

        # ---- 环路滤波 + NCO ----
        v_code = dll.step(dll_err[m])
        v_carr = pll.step(pll_err[m])

        # 载波环先更新：它的多普勒估计要拿去辅助码环
        carr_freq = f_if + init_doppler_hz + v_carr / TWO_PI    # 文档第六节

        # 码 NCO：标称码速率（可选地被载波多普勒牵引）减去码环修正量
        if cfg.carrier_aiding:
            code_rate_ref = CODE_RATE_HZ * (1.0 + (carr_freq - f_if) / L1_HZ)
        else:
            code_rate_ref = CODE_RATE_HZ
        code_freq = code_rate_ref - v_code

        # ---- 推进本地相位到下一个积分块 ----
        # 标称部分（每块整数个码周期）在 mod spc 下自动抵消，只留偏差量
        tau = (tau + (code_freq - CODE_RATE_HZ) * pdi_s * samples_per_chip) % spc
        carr_phase = (carr_phase + TWO_PI * carr_freq * block * ts) % TWO_PI

        code_freqs[m] = code_freq
        dopplers[m] = carr_freq - f_if
        taus[m] = tau

    return TrackingResult(
        prn=prn, ie=ie, qe=qe, ip=ip, qp=qp, il=il, ql=ql,
        dll_error_chips=dll_err, pll_error_rad=pll_err,
        code_freq_hz=code_freqs, doppler_hz=dopplers,
        code_phase_samples=taus, config=cfg,
    )


def track_all(
    data: np.ndarray,
    fs: float,
    f_if: float,
    acquisitions: list[dict],
    n_ms: int,
    cfg: TrackingConfig | None = None,
) -> list[TrackingResult]:
    """对捕获给出的每颗星各开一个通道跟踪。

    acquisitions 用 acquisition.acquire_all / detect 的输出即可，
    需要包含 prn、code_phase_samples、doppler_hz 三个字段。
    若做过精频估计（refine_doppler），把结果放在 doppler_refined_hz 字段，
    这里会优先用它 —— 直接用粗捕获值会因频差过大而牵不进去。
    """
    results = []
    for acq in acquisitions:
        doppler = acq.get("doppler_refined_hz", acq["doppler_hz"])
        results.append(
            track_channel(data, fs, f_if, acq["prn"],
                          init_code_phase_samples=int(acq["code_phase_samples"]),
                          init_doppler_hz=float(doppler),
                          n_ms=n_ms, cfg=cfg)
        )
    return results
